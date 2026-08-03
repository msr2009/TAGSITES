from shiny import reactive, ui, render, module
import base64, json, os, re
from modules.json_card import json_upload_card as _json_upload_card

from utils.results import (
    load_data_from_json,
    load_run_metadata,
    load_reagents_df,
    build_plot_payload,
    assign_task_colors,
    residue_colors_for_track,
    residue_colors_for_annotations,
    residue_colors_for_phobius,
    residue_colors_for_isoforms,
    residue_colors_for_patches,
    residue_colors_jet,
    _guess_analysis_type,
    _pick_gradient_cmap,
    _build_isoform_pane,
    _VIRIDIS, _PLASMA, _COOL, _BWR,
)
from utils.scoring import (
    score_tag_sites, failed_flags, mask_reasons, score_cmap_hex, score_max,
    load_scoring_config, pick_suggested_sites, position_mask_criterion, with_extra_criteria,
)
from utils.tag_filters import isoform_key, isoform_allowed_positions, topology_tag_allowed_positions
from config import RESULTS_TYPE_DICT, DOMAIN_SOURCE_COLORS

# color for residues masked outright (issue #38) — distinct from the white->green
# score gradient; must match the client-side MASKED_HEX in www/tagsites.js
_MASKED_HEX = "#d9b3b3"


def _safe_id(s):
    """Strip non-alphanumeric characters for use in a Shiny input ID."""
    return re.sub(r'[^A-Za-z0-9]', '_', str(s))


def _job_label(col, analysis):
    """Strip the trailing '_<analysis>' suffix from a task column name, e.g. 'BROAD_blast' -> 'BROAD'."""
    suffix = "_" + analysis
    return col[:-len(suffix)] if col.lower().endswith(suffix) else col

# pLDDT 4-band legend items (values stored 0-1 after /100 at load time)
_PLDDT_LEGEND = [
    {"color": "#0053D6", "label": "Very high (≥0.90)"},
    {"color": "#65CBF3", "label": "Confident (0.70–0.90)"},
    {"color": "#FFDB13", "label": "Low (0.50–0.70)"},
    {"color": "#FF7D45", "label": "Very low (<0.50)"},
]


def _build_colors_and_legend(track, task_name, scheme, aa_df, range_df, seq_len, hex_color,
                             iso_result=None):
    """Single entry point: compute per-residue structure colors AND the matching legend.

    Returns (colors, legend) where legend is a dict ready for JSON serialization,
    or (None, None) when the track cannot be resolved.
    """
    # ── Annotation coloring (domains / modifications / Phobius) ─────────────────
    if track == "__domains__":
        if range_df is None or not seq_len:
            return None, None
        colors, items = residue_colors_for_annotations(range_df, seq_len)
        return colors, {"type": "categorical", "items": items}

    if track == "__phobius__":
        if range_df is None or not seq_len:
            return None, None
        colors, items = residue_colors_for_phobius(range_df, seq_len)
        return colors, {"type": "categorical", "items": items}

    if track == "__isoforms__":
        if iso_result is None or not seq_len:
            return None, None
        colors, items = residue_colors_for_isoforms(iso_result, seq_len)
        return colors, {"type": "categorical", "items": items}

    if track == "__hydrophobic_patch__":
        if range_df is None or not seq_len:
            return None, None
        colors, items = residue_colors_for_patches(range_df, seq_len)
        return colors, {"type": "categorical", "items": items}

    # ── Continuous track coloring ────────────────────────────────────────────────
    if aa_df is None or task_name not in aa_df.columns:
        return None, None

    colors = residue_colors_for_track(aa_df, task_name, hex_color)

    # pLDDT categorical (4-band) — excludes SASA
    if _guess_analysis_type(task_name) == "plddt" and not task_name.endswith("_sasa"):
        legend = {"type": "categorical", "items": _PLDDT_LEGEND}
    else:
        _CMAP_NAMES = {
            id(_VIRIDIS): "viridis", id(_PLASMA): "plasma",
            id(_COOL): "cool",       id(_BWR): "bwr",
        }
        cmap = _pick_gradient_cmap(task_name)
        cmap_name = _CMAP_NAMES.get(id(cmap), "viridis")
        # hydrophobic exposure uses named endpoints instead of 0/1
        vmin, vmax = ("polar", "hydrophobic") if cmap is _BWR else (0, 1)
        legend = {"type": "gradient", "label": task_name, "vmin": vmin, "vmax": vmax,
                  "cmap": cmap_name}

    return colors, legend


@module.server
def results_server(input, output, session, shared_json, shared_sites, shared_results_trigger=None):

    # ── Per-run data ────────────────────────────────────────────────────────────
    aa_data       = reactive.Value()     # continuous scores DataFrame
    range_data    = reactive.Value()     # range annotations DataFrame
    iso_data      = reactive.Value(None) # isoforms JSON dict (see derive_isoforms.py), or None
    reagents_data = reactive.Value(None) # reagents-TSV subset used for scoring (see load_reagents_df)
    aln_meta      = reactive.Value([])   # list of (path, task_name, params)
    run_name      = reactive.Value(None)
    run_meta      = reactive.Value({})   # {query_seq, pdb_path, seq_len}
    task_colors   = reactive.Value({})   # task_name → hex color (stable across renders)
    # bumped once per _do_load; lets _push_rescore skip the redundant incremental send
    # that would otherwise immediately follow the full tagsites_set_plot on load. A
    # plain mutable counter, NOT a reactive.Value: _do_load runs inside the
    # auto_load_results effect, and a reactive.Value read-then-written by the same
    # effect that produced it (x.set(x.get() + 1)) is a classic Shiny self-invalidation
    # loop — the effect becomes a dependent of the very value it just changed, so it
    # reruns forever. _push_rescore only needs this as a plain, non-reactive marker.
    _plot_epoch = [0]

    # ── Color-by choices (drives the button row above the structure) ────────────
    color_by_choices = reactive.Value({})  # val → label dict, populated on load

    # canonical on-disk path (not the Shiny upload temp path) — used for saves
    json_path = reactive.Value(None)

    # ── Advanced / Rescoring: isoform + topology restrictions (session-only) ────

    def _checked(cid, default=True):
        """Value of a dynamically rendered checkbox, or `default` before it exists."""
        try:
            v = input[cid]()
        except Exception:
            return default
        return default if v is None else bool(v)

    def _phobius_labels():
        """Distinct Phobius topology labels present in the current run, in first-seen order."""
        range_df = range_data.get()
        if range_df is None or range_df.empty:
            return []
        phobius = range_df[range_df["source"] == "Phobius"]
        return list(dict.fromkeys(phobius["description"]))

    @reactive.calc
    def user_filter_spec():
        """(visible_keys, tagged_keys, constitutive, checked_topology) from the advanced controls."""
        isoforms = (iso_data.get() or {}).get("isoforms", [])
        visible_keys, tagged_keys = set(), set()
        for i, iso in enumerate(isoforms):
            k = isoform_key(iso, i)
            if _checked(f"iso_show_{_safe_id(k)}", True):
                visible_keys.add(k)
            if _checked(f"iso_tag_{_safe_id(k)}", True):
                tagged_keys.add(k)
        constitutive = _checked("iso_constitutive_only", False)
        checked_topology = {lbl for lbl in _phobius_labels()
                            if _checked(f"topology_tag_{_safe_id(lbl)}", True)}
        return visible_keys, tagged_keys, constitutive, checked_topology

    def _visible_iso_result():
        """iso_data with hidden isoforms filtered out (the query isoform is never hidden)."""
        iso_result = iso_data.get()
        if not iso_result:
            return iso_result
        visible_keys, _, _, _ = user_filter_spec()
        isoforms = [iso for i, iso in enumerate(iso_result.get("isoforms", []))
                   if iso.get("is_query") or isoform_key(iso, i) in visible_keys]
        return {**iso_result, "isoforms": isoforms}

    @reactive.calc
    def site_scores():
        """Tag-site scores for the loaded run under the current advanced-panel filters.

        Returns (scores_df, merged_config) — the merged config (base + any synthesized
        user-restriction criteria) must be reused by every mask_reasons()/failed_flags()
        call downstream, or the tooltip silently loses the user-restriction labels.
        """
        aa_df, range_df = aa_data.get(), range_data.get()
        meta = run_meta.get() or {}
        seq_len = meta.get("seq_len", 0)
        cfg = load_scoring_config()  # reloaded each time so config edits take effect

        isoforms = (iso_data.get() or {}).get("isoforms", [])
        visible_keys, tagged_keys, constitutive, checked_topology = user_filter_spec()
        all_topology = set(_phobius_labels())
        all_positions = set(range(1, seq_len + 1))

        # isoform and topology restrictions are computed independently and each becomes
        # its own mask criterion; since mask criteria OR together (any firing masks the
        # position), the net effect is exactly their intersection — same as ANDing the
        # two restrictions directly
        extra = []
        iso_allowed = isoform_allowed_positions(isoforms, visible_keys, tagged_keys,
                                                constitutive, seq_len)
        if iso_allowed is not None:
            label = "not constitutive" if constitutive else "not in selected isoforms"
            extra.append(position_mask_criterion(all_positions - iso_allowed,
                                                  "user_isoform_mask", label))
        topo_allowed = topology_tag_allowed_positions(range_df, all_topology, checked_topology, seq_len)
        if topo_allowed is not None:
            label = ("no topology tagged" if not checked_topology
                     else "outside " + ", ".join(sorted(checked_topology)))
            extra.append(position_mask_criterion(all_positions - topo_allowed,
                                                  "user_topology_mask", label))
        cfg = with_extra_criteria(cfg, extra)

        scores_df = score_tag_sites(aa_df, range_df, meta.get("query_seq", ""),
                                    reagents_data.get(), cfg)
        return scores_df, cfg

    # ── Load results ────────────────────────────────────────────────────────────

    async def _do_load(json_content):
        """Parse a loaded JSON dict and push all results to the UI."""
        aa_df, range_df, alns, iso_result = load_data_from_json(json_content, RESULTS_TYPE_DICT)
        meta = load_run_metadata(json_content)
        aa_data.set(aa_df)
        range_data.set(range_df)
        aln_meta.set(alns)
        iso_data.set(iso_result)
        run_name.set(json_content["global"]["run_name"])
        run_meta.set(meta)
        task_colors.set(assign_task_colors(aa_df) if aa_df is not None else {})

        # tag-site suggestion scoring (see utils/scoring.py, issue #32) — reagents
        # data is optional and only present once the Reagents task has run; site_scores()
        # itself reloads scores.config.json each call so config edits take effect without
        # restarting the app
        reagents_data.set(load_reagents_df(json_content))

        # build color-by button choices, grouped by analysis type (BLAST jobs together,
        # then pLDDT jobs together) rather than left in raw column order
        # pLDDT tracks use the standard AlphaFold 4-band categorical scheme only
        choices = {"(none)": "N→C"}
        if aa_df is not None:
            blast_cols, plddt_cols, other_cols = [], [], []
            for col in aa_df.columns[1:]:
                if _guess_analysis_type(col) == "blast":
                    blast_cols.append(col)
                elif _guess_analysis_type(col) == "plddt" and not col.endswith("_sasa"):
                    plddt_cols.append(col)
                elif col.endswith("_sasa"):
                    other_cols.append((col, "Solv Access"))
                elif col.endswith("_hydro"):
                    other_cols.append((col, "Hydrophobic exposure"))
                else:
                    other_cols.append((col, col))
            for label, cols in (("BLAST", blast_cols), ("pLDDT", plddt_cols)):
                for col in cols:
                    job = _job_label(col, _guess_analysis_type(col))
                    choices[col] = label if len(cols) == 1 else f"{label} ({job})"
            for col, name in other_cols:
                choices[col] = name
        if range_df is not None and not range_df.empty:
            if not range_df[range_df["source"].isin({"Pfam", "modification", "UniProt", "UniProt_site"})].empty:
                choices["__domains__"] = "Domains"
            if not range_df[range_df["source"] == "Phobius"].empty:
                choices["__phobius__"] = "Phobius"
            if not range_df[range_df["source"] == "hydrophobic_patch"].empty:
                choices["__hydrophobic_patch__"] = "Hydrophobic patches"
        if iso_result and len(iso_result.get("isoforms", [])) > 1:
            choices["__isoforms__"] = "Isoforms"
        # NOTE: deliberately not calling site_scores() here — this effect (via _do_load)
        # is what sets aa_data/range_data/run_meta/iso_data/reagents_data, which
        # site_scores() depends on; reading it from the same execution that just wrote
        # those values made this effect a dependent of its own writes, causing an
        # infinite invalidate/rerun loop. seq_len > 0 is exactly score_tag_sites'
        # condition for producing a non-empty frame (see utils.scoring.score_tag_sites).
        if meta.get("seq_len", 0) > 0:
            choices["__score__"] = "Score"
        color_by_choices.set(choices)
        ui.update_select("color_by", choices=choices, selected="(none)")

        await _send_plot(aa_df, range_df, meta, json_content["global"]["run_name"], iso_result)
        _plot_epoch[0] += 1  # marks the full payload just sent above as current
        # always reset structure colors on load — don't rely on color_by reactive firing
        # (it won't fire if color_by was already "(none)" before this load)
        jet_colors = residue_colors_jet(meta.get("seq_len", 0))
        await session.send_custom_message("tagsites_set_colors",
                                          {"colors": jet_colors, "legend": {"type": "rainbow"}})
        if meta.get("pdb_path"):
            await _send_struct(meta["pdb_path"])
        else:
            # no AFDB structure for this run — clear out any structure left over
            # from a previously loaded run rather than leaving it stale on screen
            await session.send_custom_message("tagsites_clear_struct", {})
        # record the canonical on-disk path so _save_sites can write back
        wd = json_content["global"].get("working_dir", "")
        rn_str = json_content["global"].get("run_name", "")
        canon = os.path.join(wd, rn_str + ".run.json") if wd and rn_str else None
        json_path.set(canon)

        # restore previously saved site selection (or start fresh)
        saved = json_content["global"].get("selected_sites", [])
        shared_sites.set(sorted(int(s) for s in saved))

    @reactive.effect
    @reactive.event(input.json_file_input)
    def load_results():
        """Propagate a Results-tab JSON upload to shared state.

        Validates the file first, then sets shared_json so that auto_load_results
        (and the Progress / Reagents tabs) all react to the same path — same path
        as when JSON is loaded through any other tab.
        """
        info = input.json_file_input()
        if not info:
            return
        path = info[0]["datapath"]
        # bundle ZIPs are valid — the interceptor in server.py handles extraction
        from utils.bundle import validate_run_json
        error = validate_run_json(path)
        if error:
            ui.notification_show(error, type="error", duration=8)
            return
        shared_json.set(path)

    @reactive.effect
    async def auto_load_results():
        """Auto-plot results whenever the JSON path changes or a run completes."""
        path = shared_json.get()
        # take a dependency on the completion trigger so we reload after tasks finish
        # even when the JSON path hasn't changed
        if shared_results_trigger is not None:
            shared_results_trigger.get()
        if not path:
            return
        try:
            with open(path) as f:
                json_content = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, UnicodeDecodeError):
            return
        await _do_load(json_content)

    @reactive.effect
    def _save_sites():
        """Write selected_sites back to the run JSON whenever they change."""
        sites = shared_sites.get()
        path  = json_path.get()
        if not path or not os.path.exists(path):
            return
        try:
            with open(path) as f:
                j = json.load(f)
            j.setdefault("global", {})["selected_sites"] = sites
            with open(path, "w") as f:
                json.dump(j, f, indent=2)
        except Exception:
            pass

    def _input_names():
        """Build the dict of namespaced Shiny input IDs to send to JS."""
        return {
            "residue_click": session.ns("residue_click"),
            "struct_click":  session.ns("struct_click"),
            "remove_site":   session.ns("remove_site"),
        }

    def _score_payload(scores, cfg, seq_len):
        """Build the score* payload keys shared by the full plot and incremental rescore.

        `cfg` must be the merged config that produced `scores` (base criteria plus any
        synthesized user-restriction ones) or mask_reasons() will silently drop the
        user-restriction tooltip labels.
        """
        if scores is None or scores.empty or not seq_len:
            return {"scoreTrack": [], "scoreFlags": [], "scoreMasked": [],
                    "scoreMaskReasons": [], "scoreMax": 1, "suggestedSites": []}
        by_pos = scores["score"]
        flags_by_pos = dict(zip(scores.index, failed_flags(scores, cfg)))
        masked_by_pos = scores["masked"]
        mask_reasons_by_pos = dict(zip(scores.index, mask_reasons(scores, cfg)))
        return {
            "scoreTrack": [int(by_pos[p]) if p in by_pos.index else None
                          for p in range(1, seq_len + 1)],
            "scoreFlags": [flags_by_pos.get(p, []) for p in range(1, seq_len + 1)],
            "scoreMasked": [bool(masked_by_pos[p]) if p in masked_by_pos.index else False
                           for p in range(1, seq_len + 1)],
            "scoreMaskReasons": [mask_reasons_by_pos.get(p, []) for p in range(1, seq_len + 1)],
            "scoreMax": max(1, int(score_max(cfg))),
            "suggestedSites": pick_suggested_sites(scores),
        }

    async def _send_plot(aa_df, range_df, meta, title, iso_result=None):
        """Build and send all plot data to the native canvas renderer.

        Called from _do_load, so this computes scores directly with score_tag_sites()
        rather than through the site_scores() calc — site_scores() depends (via
        user_filter_spec) on the Advanced/Rescoring checkboxes, and _do_load is what
        (re)renders those checkboxes each load; calling site_scores() from here would
        make the load effect a listener of its own dynamically-recreated inputs and
        loop forever. No user restriction has been set yet at load time anyway (the
        panel is session-only and starts unrestricted), so this matches site_scores()'s
        unfiltered case exactly. See _push_rescore for the filtered incremental path.
        """
        payload = build_plot_payload(aa_df, range_df, iso_result=iso_result, title=title)
        payload["seq"] = meta.get("query_seq", "")
        payload["inputs"] = _input_names()

        # tag-site score heatmap row (see utils/scoring.py, issue #32)
        cfg = load_scoring_config()
        scores = score_tag_sites(aa_df, range_df, meta.get("query_seq", ""),
                                 reagents_data.get(), cfg)
        payload.update(_score_payload(scores, cfg, meta.get("seq_len", 0)))

        await session.send_custom_message("tagsites_set_plot", payload)

    _last_full_plot_epoch = {"v": None}

    @reactive.effect
    async def _push_rescore():
        """Push an incremental score/isoform-pane update when the advanced controls change.

        Skips the first invalidation after each _do_load — that rescore was already
        covered by the full tagsites_set_plot send — so only later advanced-panel edits
        trigger this (issue: don't wipe zoom/committed sites with a full plot re-send).
        """
        scores, cfg = site_scores()
        epoch = _plot_epoch[0]  # plain read, not reactive — see _plot_epoch's definition
        if _last_full_plot_epoch["v"] != epoch:
            _last_full_plot_epoch["v"] = epoch
            return
        meta = run_meta.get() or {}
        payload = _score_payload(scores, cfg, meta.get("seq_len", 0))
        payload["isoformPane"] = _build_isoform_pane(_visible_iso_result(), range_data.get())
        await session.send_custom_message("tagsites_set_scores", payload)

    async def _send_struct(pdb_path):
        """Read PDB file and send to 3Dmol viewer; clears the viewer if unreadable."""
        try:
            with open(pdb_path) as f:
                pdb_str = f.read()
        except OSError:
            await session.send_custom_message("tagsites_clear_struct", {})
            return
        await session.send_custom_message("tagsites_init_struct", {
            "pdb": pdb_str,
            "inputs": _input_names(),
        })

    # ── Click/selection handlers ─────────────────────────────────────────────────

    def _click_pos(raw):
        """Unwrap a [pos] click payload; JS wraps positions in a fresh array each
        click so re-clicking the same residue still invalidates server-side
        (py-shiny skips updates when the new value `is` the cached small int).
        """
        if raw is None:
            return None
        if isinstance(raw, (list, tuple)):
            raw = raw[0]
        return int(raw)

    @reactive.effect
    @reactive.event(input.residue_click)
    def on_residue_click():
        """Add a residue to the committed (green) site list."""
        pos = _click_pos(input.residue_click())
        if pos is None:
            return
        sites = sorted(set(shared_sites.get()) | {pos})
        shared_sites.set(sites)

    @reactive.effect
    @reactive.event(input.struct_click)
    def on_struct_click():
        """Add a residue to the committed (green) site list from a 3D structure click."""
        pos = _click_pos(input.struct_click())
        if pos is None:
            return
        sites = sorted(set(shared_sites.get()) | {pos})
        shared_sites.set(sites)

    @reactive.effect
    @reactive.event(input.add_nterm_button)
    def on_add_nterm():
        """Add position 1 (N-terminus) to committed sites."""
        sites = sorted(set(shared_sites.get()) | {1})
        shared_sites.set(sites)

    @reactive.effect
    @reactive.event(input.add_cterm_button)
    def on_add_cterm():
        """Add the last residue (C-terminus) to committed sites."""
        meta = run_meta.get() or {}
        seq_len = meta.get("seq_len", 0)
        if not seq_len:
            return
        sites = sorted(set(shared_sites.get()) | {seq_len})
        shared_sites.set(sites)

    @reactive.effect
    @reactive.event(input.add_suggested_button)
    def on_add_suggested():
        """Add the top-scoring, well-spaced tag sites (see utils.scoring.score_tag_sites)."""
        scores, _ = site_scores()
        picked = pick_suggested_sites(scores) if scores is not None and not scores.empty else []
        if not picked:
            return
        sites = sorted(set(shared_sites.get()) | set(picked))
        shared_sites.set(sites)

    @reactive.effect
    @reactive.event(input.clear_highlights_button)
    def on_clear():
        """Clear all committed tag sites."""
        shared_sites.set([])

    @reactive.effect
    @reactive.event(input.remove_site)
    def on_remove_site():
        """Remove a committed site via its chip × button."""
        pos = input.remove_site()
        if pos is None:
            return
        pos = int(pos)
        sites = [s for s in shared_sites.get() if s != pos]
        shared_sites.set(sites)

    # ── Sync JS highlight state whenever Python state changes ───────────────────

    @reactive.effect
    async def sync_js_states():
        """Push committed state to both JS surfaces."""
        await session.send_custom_message("tagsites_set_states", {
            "committed": sorted(shared_sites.get()),
        })

    # ── Color-by handler ─────────────────────────────────────────────────────────

    @reactive.effect
    async def sync_js_colors():
        """Compute per-residue colors for the selected track and send to 3D viewer.

        Not gated with @reactive.event(input.color_by) — the __score__ branch below reads
        site_scores(), so a live rescore (advanced-panel change) also refreshes the
        structure coloring, not just an explicit color_by change.
        """
        track = input.color_by()
        if not track or track == "(none)":
            meta = run_meta.get()
            seq_len = meta.get("seq_len", 0) if meta else 0
            colors = residue_colors_jet(seq_len) if seq_len else []
            await session.send_custom_message("tagsites_set_colors",
                                              {"colors": colors,
                                               "legend": {"type": "rainbow"}})
            return

        if track == "__score__":
            scores, _ = site_scores()
            meta = run_meta.get()
            seq_len = meta.get("seq_len", 0) if meta else 0
            vmax = max(1, int(scores["score"].max())) if scores is not None and not scores.empty else 1
            colors = []
            if scores is not None and not scores.empty and seq_len:
                by_pos = scores["score"]
                masked_by_pos = scores["masked"]
                for p in range(1, seq_len + 1):
                    if p not in by_pos.index:
                        colors.append("#e8e8e8")
                    elif masked_by_pos[p]:
                        colors.append(_MASKED_HEX)
                    else:
                        colors.append(score_cmap_hex(by_pos[p] / vmax))
            legend = {"type": "gradient", "cmap": "viridis", "label": "Score",
                     "vmin": 0, "vmax": vmax, "maskedColor": _MASKED_HEX, "maskedLabel": "masked"}
            await session.send_custom_message("tagsites_set_colors",
                                              {"colors": colors, "legend": legend})
            return

        task_name, scheme = track.rsplit(":", 1) if ":" in track else (track, "categorical")
        meta = run_meta.get()
        seq_len = meta.get("seq_len", 0) if meta else 0
        hex_color = task_colors.get().get(task_name, "#888888")

        colors, legend = _build_colors_and_legend(
            track, task_name, scheme,
            aa_data.get(), range_data.get(), seq_len, hex_color,
            iso_result=_visible_iso_result(),
        )
        if colors is None:
            return

        await session.send_custom_message("tagsites_set_colors",
                                          {"colors": colors, "legend": legend})

    # ── Outputs ──────────────────────────────────────────────────────────────────

    @render.ui
    def json_upload_card():
        """Render the JSON upload card; collapsed and warning-free when data is loaded."""
        return _json_upload_card(
            "json_file_input", "ts-json-body", "Upload results JSON/ZIP",
            run_name.get() is not None,
        )

    @render.ui
    def color_buttons_ui():
        """Render color-by buttons above the structure viewer."""
        choices = color_by_choices.get()
        if not choices:
            return ui.div()
        color_by_input_id = session.ns("color_by")
        _isoforms_tooltip = (
            "Isoform coverage from UniProt annotation where available, otherwise inferred from "
            "same-organism BLAST hits. Segments reflect protein sequence coverage, not genomic "
            "exon boundaries."
        )
        buttons = [
            ui.tags.button(
                label,
                onclick=f"tsSetColorBy(this, '{val}', '{color_by_input_id}')",
                class_="btn btn-sm ts-colorby-btn" + (" ts-colorby-active" if val == "(none)" else ""),
                title=_isoforms_tooltip if val == "__isoforms__" else None,
            )
            for val, label in choices.items()
        ]
        return ui.div(*buttons, class_="ts-colorby-row")

    @render.ui
    def chosen_sites_display():
        """Render chosen-site chips with × remove buttons."""
        sites = shared_sites.get()
        if not sites:
            return ui.div(
                ui.span("No sites chosen yet — click residues in the plot or structure to select.",
                        class_="ts-empty-hint"),
                id="ts-sites-box",
            )
        meta = run_meta.get() or {}
        seq = meta.get("query_seq", "")
        remove_id = session.ns("remove_site")
        chips = []
        for pos in sites:
            aa = seq[pos - 1] if seq and 0 < pos <= len(seq) else "?"
            label = f"{aa}{pos}"
            chip = ui.span(
                label,
                ui.tags.button(
                    "×",
                    class_="ts-chip-remove",
                    onclick=f"tsRemoveSite({pos}, '{remove_id}')",
                    title=f"Remove {label}",
                ),
                class_="ts-site-chip",
            )
            chips.append(chip)
        return ui.div(*chips, id="ts-sites-box")

    @render.ui
    def advanced_panel():
        """Advanced / Rescoring accordion — isoform + topology tag restrictions.

        Restricts scoring/suggested-sites/structure coloring only — never which sites get
        reagents designed (Reagents tab annotates isoform specificity instead; see
        modules/reagents_server.py). Hidden entirely when there's nothing to restrict.

        Deliberately depends only on iso_data/range_data (i.e. re-renders on a new run
        load), NOT on site_scores() — site_scores() itself depends on the checkboxes this
        function creates, so reading it here would re-render the checkboxes every time
        their own value changes, which re-fires their "set" event and never settles. The
        live position count instead lives in the separate filter_summary output below.
        """
        iso_result = iso_data.get() or {}
        isoforms = sorted(iso_result.get("isoforms", []), key=lambda x: x.get("length", 0),
                          reverse=True)
        topology_labels = _phobius_labels()
        if len(isoforms) <= 1 and not topology_labels:
            return ui.div()

        body = []

        if len(isoforms) > 1:
            # Show before Tag so Tag ends up the rightmost column in both matrices —
            # matches the Phobius matrix below (Name, Tag) at the same 340px width, so
            # the Tag checkboxes line up vertically between the two matrices.
            header = ui.div(
                ui.span("Isoform (Uniprot ID)", class_="ts-matrix-name-label"),
                ui.span("Show", class_="ts-matrix-col-label"),
                ui.span("Tag", class_="ts-matrix-col-label"),
                class_="ts-iso-matrix-row ts-matrix-header",
            )
            iso_rows = []
            for i, iso in enumerate(isoforms):
                k = _safe_id(isoform_key(iso, i))
                # UniProt accession is the preferred row label — it's the stable,
                # unambiguous identifier; fall back to name/index for BLAST-inferred
                # isoforms that have no accession (see derive_isoforms.py)
                label = iso.get("accession") or iso.get("name") or f"isoform {i + 1}"
                length = iso.get("length")
                if length:
                    label += f" ({length} aa)"
                is_query = bool(iso.get("is_query"))
                # the query row's Show is a real (Bootstrap-styled) checkbox, same as the
                # other rows, just disabled client-side below — a raw unstyled <input>
                # looked visually different (no Bootstrap checkbox skin) from the rest
                show_class = "ts-query-show-cell" if is_query else ""
                iso_rows.append(ui.div(
                    ui.span(label + (" (query)" if is_query else ""), class_="ts-matrix-name"),
                    ui.div(ui.input_checkbox(f"iso_show_{k}", "", value=True, width="auto"),
                          class_=show_class),
                    ui.div(ui.input_checkbox(f"iso_tag_{k}", "", value=True, width="auto"),
                          class_="ts-iso-tag-cell"),
                    class_="ts-iso-matrix-row",
                ))
            body.append(ui.div(
                header, *iso_rows,
                id="ts-iso-matrix",
                class_="ts-matrix",
            ))
            body.append(ui.div(
                ui.input_switch("iso_constitutive_only", "Constitutive only", value=False),
                ui.span("— overrides the Tag column above (disables those checkboxes); "
                        "restricts tagging to sequence present in every shown isoform.",
                        class_="ts-filter-hint"),
                class_="ts-constitutive-toggle",
            ))
            # purely client-side: disables the isoform matrix's Tag checkboxes while
            # Constitutive-only is on, and permanently disables the query row's Show
            # checkbox (it's always visible; the input is never read server-side).
            # Not done by re-rendering this output on the switch's value (see the
            # docstring above — that would recreate the switch itself and reset it,
            # since this whole function would become a dependent of its own
            # dynamically-created input).
            toggle_id = session.ns("iso_constitutive_only")
            body.append(ui.tags.script(ui.HTML(f"""
                (function() {{
                  var toggle = document.getElementById("{toggle_id}");
                  var tagBoxes = document.querySelectorAll("#ts-iso-matrix .ts-iso-tag-cell input[type=checkbox]");
                  function sync() {{
                    tagBoxes.forEach(function(cb) {{ cb.disabled = toggle.checked; }});
                  }}
                  if (toggle) {{ toggle.addEventListener("change", sync); sync(); }}
                  document.querySelectorAll("#ts-iso-matrix .ts-query-show-cell input[type=checkbox]")
                    .forEach(function(cb) {{ cb.disabled = true; cb.title = "The query isoform can't be hidden"; }});
                }})();
            """)))

        if len(isoforms) > 1 and topology_labels:
            body.append(ui.tags.hr(class_="ts-filter-divider"))

        if topology_labels:
            topo_header = ui.div(
                ui.span("Topology (Phobius)", class_="ts-matrix-name-label"),
                ui.span("Tag", class_="ts-matrix-col-label"),
                class_="ts-topo-matrix-row ts-matrix-header",
            )
            topo_rows = [
                ui.div(
                    ui.span(lbl, class_="ts-matrix-name"),
                    ui.div(ui.input_checkbox(f"topology_tag_{_safe_id(lbl)}", "", value=True,
                                             width="auto")),
                    class_="ts-topo-matrix-row",
                )
                for lbl in topology_labels
            ]
            body.append(ui.div(
                topo_header, *topo_rows,
                class_="ts-matrix",
            ))

        body.append(ui.output_text("filter_summary"))

        return ui.accordion(
            ui.accordion_panel("Rescoring (Advanced)", *body),
            open=False, id="advanced_accordion", class_="ts-advanced-accordion",
        )

    @render.text
    def filter_summary():
        """Live '412 / 854 positions available' line — separate from advanced_panel

        so a rescore updates the count without re-rendering (and re-triggering) the
        Show/Tag checkboxes themselves.
        """
        scores, _ = site_scores()
        meta = run_meta.get() or {}
        seq_len = meta.get("seq_len", 0)
        if scores is None or scores.empty or not seq_len:
            return ""
        available = int((~scores["masked"]).sum())
        return f"{available} / {seq_len} positions available"

    @render.ui
    def alignments_container():
        """Render blast alignment PDFs as collapsible accordion panes."""
        alns = aln_meta.get()
        if not alns:
            return ui.div()

        panels = []
        for aln_path, task_name, params in alns:
            pdf_path = aln_path.removesuffix(".aln") + ".pdf"

            param_rows = [
                ui.tags.tr(ui.tags.td(k), ui.tags.td(str(v)))
                for k, v in params.items()
            ]
            param_table = ui.tags.table(
                ui.tags.tbody(*param_rows),
                class_="ts-aln-params",
            ) if param_rows else ui.div()

            if os.path.exists(pdf_path):
                with open(pdf_path, "rb") as f:
                    b64 = base64.b64encode(f.read()).decode("ascii")
                embed = ui.tags.iframe(
                    src=f"data:application/pdf;base64,{b64}",
                    style="width:100%; height:400px; border:none;",
                )
                body = ui.div(
                    param_table,
                    ui.div(embed, class_="ts-aln-svg-wrap"),
                )
            else:
                body = ui.div(
                    param_table,
                    ui.p(f"Alignment image not found: {os.path.basename(pdf_path)}",
                         style="color:#c00;"),
                )

            panels.append(ui.accordion_panel(task_name, body))

        if not panels:
            return ui.div()
        return ui.div(
            ui.h4("Alignments"),
            ui.accordion(*panels, id="ts_aln_accordion", open=False),
        )
