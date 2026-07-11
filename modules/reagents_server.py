"""reagents_server.py — Shiny module server for the Design Reagents tab.

Loads the run JSON (via shared_json) and the precomputed .reagents.tsv, then
presents one accordion card per chosen tag site (shared_sites) with ASCII
diagrams, guide-RNA selection checkboxes, and output downloads.
"""

import io
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd
from shiny import module, reactive, render, ui

# ensure scripts/ is importable from the app context
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from reagent_sequences import (
    ascii_diagram,
    build_ssodn,
    calc_tm,
    check_recut_risk,
    design_pcr_primers,
    load_tags,
    truncate_arms,
)
from plasmid_assembly import (
    OLIGO_MAX,
    SAPI_SITE,
    SAPI_RC,
    SAPTRAP_OVERHANGS,
    build_saptrap_site,
    find_internal_sapi,
)

from modules.json_card import json_upload_card as _json_upload_card
from modules.ui_helpers import compact_file_input

_TAGS_PATH = Path(__file__).parent.parent / "tables" / "tags.tsv"
_TAGS = load_tags(_TAGS_PATH) if _TAGS_PATH.exists() else {}

# Arm-length defaults by repair type
_ARM_DEFAULTS = {"ha": 1000, "ssodn": 50, "amplicon": 50, "plasmid": 57}


def _safe_id(s):
    """Strip non-alphanumeric characters for use in a Shiny input ID."""
    return re.sub(r'[^A-Za-z0-9]', '_', str(s))


def _parse_seq_file(path, filename=''):
    """Return an uppercase DNA string from a FASTA or GenBank file."""
    from Bio import SeqIO
    ext = os.path.splitext(filename or path)[1].lower()
    fmts = (['genbank', 'fasta'] if ext in ('.gb', '.gbk', '.genbank', '.gbff')
            else ['fasta', 'genbank'])
    for fmt in fmts:
        try:
            recs = list(SeqIO.parse(path, fmt))
            if recs:
                return str(recs[0].seq).upper()
        except Exception:
            continue
    raise ValueError("Not a valid FASTA or GenBank file")


# ── module server ──────────────────────────────────────────────────────────────

@module.server
def reagents_server(input, output, session, shared_json, shared_sites):

    # ── reactive state ────────────────────────────────────────────────────────
    run_name         = reactive.Value(None)   # str or None
    working_dir      = reactive.Value("")
    stored_arm_len   = reactive.Value(1000)   # arm_length used when TSV was generated
    reagents_df      = reactive.Value(None)   # full TSV as DataFrame or None
    selected_guides  = reactive.Value({})     # {residue_index (int): set(guide_id)}
    _template_seq    = reactive.Value("")     # PCR template, persists across repair-type switches

    # reactive path cell: drives _poll_for_tsv even when reagents_df stays None
    _tsv_path_cache = reactive.Value("")

    # ── populate tag dropdown once at session start ───────────────────────────

    @reactive.effect
    def _init_tag_choices():
        """Populate the tag dropdown on session start; no reactive deps → runs once."""
        choices = {n: n for n in _TAGS} if _TAGS else {}
        choices["Custom"] = "Custom…"
        ui.update_select("tag_name", choices=choices,
                         selected=next(iter(_TAGS), "Custom"))

    @reactive.effect
    def _load_template_file():
        """Parse an uploaded PCR template file and store the sequence."""
        try:
            info = input.pcr_template_file()
        except Exception:
            return
        if not info:
            return
        path = info[0]["datapath"]
        name = info[0]["name"]
        try:
            seq = _parse_seq_file(path, name)
        except Exception as e:
            ui.notification_show("Could not parse template file: {}".format(e),
                                 type="error", duration=6)
            return
        _template_seq.set(seq)

    # ── load run JSON whenever shared path changes ────────────────────────────

    @reactive.effect
    @reactive.event(shared_json)
    def _load_json():
        """Parse the run JSON and load the reagents TSV."""
        path = shared_json.get()
        if not path:
            return
        try:
            with open(path) as fh:
                j = json.load(fh)
        except Exception:
            return

        gb  = j.get("global", {})
        rn  = gb.get("run_name", "")
        wd  = gb.get("working_dir", "")
        run_name.set(rn if rn else None)
        working_dir.set(wd)

        # read reagent args to discover stored arm length and TSV path
        reagent_args = j.get("tasks", {}).get("REAGENTS_reagents", {}).get("args", {})
        s_arm = int(reagent_args.get("arm_length", 1000))
        stored_arm_len.set(s_arm)

        tsv_path = reagent_args.get("output", "")
        _tsv_path_cache.set(tsv_path)   # reactive set so _poll_for_tsv re-fires
        if tsv_path and os.path.exists(tsv_path):
            df = pd.read_csv(tsv_path, sep='\t')
            df["guide_id"] = (
                df["guide_strand"].astype(str)
                + "_" + df["pam_fwd_start"].astype(str)
                + "_" + df["spacer"].astype(str)
            )
            # stored arm length = minimum arm in the TSV (some sites may be shorter
            # than the requested arm_length when near sequence boundaries)
            try:
                actual_arm = int(min(
                    df["left_arm"].str.len().min(),
                    df["right_arm"].str.len().min(),
                ))
                stored_arm_len.set(actual_arm)
            except Exception:
                stored_arm_len.set(s_arm)
            reagents_df.set(df)
        else:
            reagents_df.set(None)

        # restore previously saved site selection if present
        saved_sites = gb.get("selected_sites", [])
        if saved_sites:
            shared_sites.set(sorted(int(s) for s in saved_sites))

        selected_guides.set({})

    # ── poll for reagents TSV when it wasn't ready at load time ─────────────

    @reactive.effect
    def _poll_for_tsv():
        """Retry loading the reagents TSV every 2 s until it appears on disk.

        Fires when reagents_df is None (TSV missing at load time) and we have a
        known path.  Stops automatically once the file is found and loaded.
        """
        if reagents_df.get() is not None:
            return
        path = _tsv_path_cache.get()
        if not path:
            return
        reactive.invalidate_later(2.0)   # schedule next check before file test
        if not os.path.exists(path):
            return
        try:
            df = pd.read_csv(path, sep='\t')
            df["guide_id"] = (
                df["guide_strand"].astype(str)
                + "_" + df["pam_fwd_start"].astype(str)
                + "_" + df["spacer"].astype(str)
            )
            try:
                actual_arm = int(min(
                    df["left_arm"].str.len().min(),
                    df["right_arm"].str.len().min(),
                ))
                stored_arm_len.set(actual_arm)
            except Exception:
                pass
            reagents_df.set(df)
        except Exception:
            pass

    # ── upload from the reagents tab itself ───────────────────────────────────

    @reactive.effect
    @reactive.event(input.reagents_json_file)
    def _handle_upload():
        """Push an uploaded JSON path into shared state."""
        info = input.reagents_json_file()
        if info:
            shared_json.set(info[0]["datapath"])

    # ── update arm-length default when repair type changes ────────────────────

    @reactive.effect
    def _arm_length_default():
        """Set a sensible default arm length for the chosen repair type."""
        repair = input.repair_type()
        default = _ARM_DEFAULTS.get(repair, 500)
        ui.update_numeric("arm_length", value=default)

    # ── filtered per-site dataframe ───────────────────────────────────────────

    @reactive.calc
    def sites_df():
        """DataFrame rows for the currently chosen tag sites only."""
        df    = reagents_df.get()
        sites = shared_sites.get()
        if df is None or not sites:
            return None
        mask = df["residue_index"].isin(sites)
        return df[mask].copy() if mask.any() else None

    # ── selection sync: scan per-guide checkboxes into selected_guides ─────────

    @reactive.effect
    def _sync_selection():
        """Collect checkbox states for every rendered guide into selected_guides."""
        df = sites_df()
        if df is None:
            return
        sel = {}
        for _, row in df.iterrows():
            rid = int(row["residue_index"])
            gid = str(row["guide_id"])
            cid = "sel_{}_{}".format(rid, _safe_id(gid))
            try:
                checked = input[cid]()
            except Exception:
                checked = False
            if checked:
                sel.setdefault(rid, set()).add(gid)
        selected_guides.set(sel)

    # ── insert sequence (from dropdown or custom text) ─────────────────────────

    @reactive.calc
    def _insert_seq():
        """Resolve the insert/tag DNA sequence from the dropdown or custom input."""
        name = input.tag_name() if hasattr(input, "tag_name") else "Custom"
        if name == "Custom":
            return input.tag_custom().strip().upper() if hasattr(input, "tag_custom") else ""
        return _TAGS.get(name, "").upper()

    # ── arm length validation ─────────────────────────────────────────────────

    @reactive.calc
    def _arm_ok():
        """True when the requested arm_length <= stored arm length."""
        try:
            return int(input.arm_length()) <= stored_arm_len.get()
        except Exception:
            return False

    # ── pre-compute guide content (diagrams + truncated arms) ─────────────────

    @reactive.calc
    def _guide_content():
        """Truncate arms and build ASCII diagrams for all selected guides.

        Cached per (sites_df, arm_length) — does not re-run on checkbox changes.
        Returns dict {guide_id: {left, right, diagram}} or {} on error.
        """
        df = sites_df()
        if df is None or not _arm_ok():
            return {}
        try:
            arm_len = int(input.arm_length())
        except Exception:
            return {}
        result = {}
        wt_extra = 10   # extra bp beyond arm_length shown in the WT row
        for _, row in df.iterrows():
            gid = str(row["guide_id"])
            try:
                left, right = truncate_arms(
                    str(row["left_arm"]), str(row["right_arm"]), arm_len
                )
                left_wt  = str(row["left_arm"])[-(arm_len + wt_extra):]
                right_wt = str(row["right_arm"])[:(arm_len + wt_extra)]
                diag = ascii_diagram(row, left, right,
                                     left_wt=left_wt, right_wt=right_wt)
            except Exception as e:
                left, right, diag = "", "", str(e)
            result[gid] = {
                "left": left, "right": right, "diagram": diag,
                "left_bp": len(left), "right_bp": len(right), "arm_len": arm_len,
            }
        return result

    # ── recut risk: does the insert reconstitute the guide target? ────────────

    @reactive.calc
    def _recut_risk():
        """Dict {guide_id: bool} — True if the insert seq reconstitutes the guide target at the junction."""
        df = sites_df()
        insert = _insert_seq()
        content = _guide_content()
        if df is None or not insert or not content:
            return {}
        result = {}
        for _, row in df.iterrows():
            gid = str(row["guide_id"])
            c = content.get(gid, {})
            left = c.get("left", "")
            right = c.get("right", "")
            if not left or not right:
                result[gid] = False
                continue
            result[gid] = check_recut_risk(
                left, right, insert,
                str(row["spacer"]), str(row["pam_seq"]), str(row["guide_strand"])
            )
        return result

    # ── Outputs ───────────────────────────────────────────────────────────────

    @render.ui
    def json_upload_card():
        """Collapsible JSON upload card."""
        return _json_upload_card(
            "reagents_json_file", "ts-reagents-json-body", "Upload run JSON",
            run_name.get() is not None,
        )

    @render.ui
    def run_header():
        """One-line run summary."""
        rn = run_name.get()
        if rn is None:
            return ui.span()
        df = reagents_df.get()
        n_sites = len(shared_sites.get())
        n_rows  = len(df) if df is not None else 0
        tsv_ok  = df is not None
        return ui.div(
            ui.span("Run: "),
            ui.strong(rn),
            ui.span("  |  {}{} site(s) selected".format(
                n_sites, "  |  reagents TSV loaded ({} rows)".format(n_rows) if tsv_ok
                else "  |  ⚠ reagents TSV not found — run the Reagents task first"
            )),
            class_="run-header",
        )

    @render.ui
    def tag_custom_area():
        """Conditional custom text area or sequence preview below the static tag select."""
        try:
            name = input.tag_name()
        except Exception:
            name = "Custom"
        if name == "Custom":
            return ui.input_text_area(
                "tag_custom", "Custom sequence (DNA, 5'→3', + strand)",
                value="", rows=2,
            )
        seq = _TAGS.get(name, "")
        if seq:
            return ui.p(
                "Sequence: ", ui.tags.code(seq[:60] + ("…" if len(seq) > 60 else "")),
                class_="section-hint",
            )
        return ui.span()

    @render.ui
    def repair_options():
        """Repair-type-specific options block rendered server-side."""
        repair = input.repair_type()

        if repair == "ssodn":
            return ui.div(
                ui.output_ui("ssodn_len_warning"),
            )

        if repair == "amplicon":
            tmpl = _template_seq.get()
            children = [
                compact_file_input(
                    "pcr_template_file",
                    "Template (FASTA or GenBank)",
                    accept=[".fa", ".fasta", ".fsa", ".fna", ".gb", ".gbk", ".genbank"],
                    multiple=False,
                ),
                ui.input_text_area(
                    "pcr_template", "Or paste template sequence",
                    value=tmpl, rows=3,
                ),
            ]
            if not tmpl:
                children.append(ui.div(
                    "⚠ A template sequence is required for PCR primer design.",
                    class_="ts-warn",
                ))
            children += [
                ui.input_numeric("pcr_tm", "Primer binding Tm target (°C)",
                                 value=60, min=40, max=80),
                ui.input_checkbox(
                    "pcr_phos",
                    "5′ phosphorylate one primer (lambda-exo single-strand degradation)",
                    value=False,
                ),
                ui.hr(style="margin: 0.5rem 0;"),
                ui.p("Optional 5′ extensions (e.g. Gibson assembly overhangs, "
                     "restriction sites, SapI tails). Prepended to the arm overhang.",
                     class_="section-hint"),
                ui.input_text("pcr_fwd_ext", "Forward primer 5′ extension", value="",
                              placeholder="e.g. AAACCCGGG"),
                ui.input_text("pcr_rev_ext", "Reverse primer 5′ extension", value="",
                              placeholder="e.g. TTTGGGCCC"),
            ]
            return ui.div(*children)

        if repair == "plasmid":
            df = sites_df()
            children = []

            # ── arm-length mode hint ─────────────────────────────────────────
            try:
                arm_len = int(input.arm_length())
            except Exception:
                arm_len = None
            if arm_len is not None:
                if arm_len <= OLIGO_MAX:
                    mode_msg = (
                        "Arms ≤ {} bp → ordered as annealed oligo pairs (no PCR).".format(OLIGO_MAX)
                    )
                else:
                    mode_msg = (
                        "Arms {} bp > {} bp → PCR primers with SapI sites will be designed.".format(
                            arm_len, OLIGO_MAX)
                    )
                children.append(ui.p(mode_msg, class_="section-hint"))

            # ── Tm target (shown only for PCR-arm mode) ───────────────────────
            if arm_len is not None and arm_len > OLIGO_MAX:
                children.append(
                    ui.input_numeric("saptrap_tm", "Primer binding Tm target (°C)",
                                     value=60, min=40, max=80)
                )

            # ── optional GenBank output ──────────────────────────────────────
            children.append(
                ui.input_checkbox(
                    "saptrap_genbank",
                    "Generate ApE-compatible GenBank files",
                    value=False,
                )
            )

            return ui.div(*children)

        return ui.span()   # ha: no extra options

    @render.ui
    def ssodn_len_warning():
        """Warn if the total ssODN length exceeds the 200 nt IDT synthesis limit."""
        try:
            arm_len = int(input.arm_length())
            insert  = _insert_seq()
            total   = 2 * arm_len + len(insert)
        except Exception:
            return ui.span()
        if total > 200:
            return ui.div(
                "⚠ Estimated ssODN length is ~{} nt (limit ≈ 200 nt for IDT ultramer). "
                "Reduce arm length or use a shorter tag.".format(total),
                class_="ts-warn",
            )
        return ui.span()

    @render.ui
    def arm_length_warning():
        """Error shown next to the arm length input when it exceeds stored arm length."""
        if _arm_ok():
            return ui.span()
        try:
            arm_len = int(input.arm_length())
        except Exception:
            arm_len = "?"
        return ui.div(
            "⚠ Requested arm length ({} bp) exceeds the stored arm length "
            "({} bp). Use a smaller value or rerun the Reagents task with "
            "a larger arm_length.".format(arm_len, stored_arm_len.get()),
            class_="ts-error",
        )

    @render.ui
    def site_cards():
        """One accordion panel per chosen tag site."""
        sites = shared_sites.get()
        df    = sites_df()

        if not sites:
            return ui.div(
                ui.p("No tag sites selected yet. Go to the Results tab and "
                     "click residues to choose insertion sites.",
                     class_="section-hint"),
            )
        if df is None:
            return ui.div(
                ui.p("Reagents TSV not found. Run the Reagents task on the "
                     "Progress tab first.", class_="ts-error"),
            )

        content = _guide_content()   # cached; does not re-run on checkbox changes
        recut   = _recut_risk()

        panels = []
        for rid in sorted(sites):
            site_rows = df[df["residue_index"] == rid].copy()
            if site_rows.empty:
                continue
            aa = str(site_rows.iloc[0]["amino_acid"])
            site_label = "{}{}".format(aa, rid)

            guide_divs = _build_guide_divs(site_rows, rid, content, recut)

            panels.append(
                ui.accordion_panel(
                    site_label,
                    *guide_divs,
                    value=str(rid),
                )
            )

        if not panels:
            return ui.div(ui.p("No reagent data for the chosen sites.", class_="section-hint"))

        return ui.accordion(
            *panels,
            multiple=True,
            open=False,
            class_="reagents-accordion",
            id="reagents_accordion",
        )

    @render.ui
    def download_row_top():
        """Download buttons above the site cards."""
        return _download_buttons("top")

    @render.ui
    def download_row_bottom():
        """Download buttons below the site cards."""
        return _download_buttons("bot")

    # ── download handlers ─────────────────────────────────────────────────────

    @render.download(filename=lambda: "{}_reagents.zip".format(run_name.get() or "reagents"))
    def download_sequences_top():
        """Download ZIP of all reagent files (top button)."""
        yield _build_zip()

    @render.download(filename=lambda: "{}_reagents.zip".format(run_name.get() or "reagents"))
    def download_sequences_bot():
        """Download ZIP of all reagent files (bottom button)."""
        yield _build_zip()

    # ── private helpers ───────────────────────────────────────────────────────

    def _arm_sapi_html(seq, label):
        """HTML for an arm sequence with SapI recognition sites highlighted in red."""
        from html import escape as _esc
        seq_up = seq.upper()
        sites = set()
        for site_str in [SAPI_SITE, SAPI_RC]:
            i = 0
            while True:
                p = seq_up.find(site_str, i)
                if p == -1:
                    break
                for k in range(p, p + 7):
                    sites.add(k)
                i = p + 1

        parts = [_esc(label) + ': <code>']
        for i, ch in enumerate(seq):
            if i in sites:
                parts.append('<span style="color:#c0392b;font-weight:bold">' + _esc(ch) + '</span>')
            else:
                parts.append(_esc(ch))
        parts.append('</code>')
        return ''.join(parts)

    def _build_guide_divs(site_rows, rid, content, recut_risk):
        """Build guide-panel divs for one tag site.

        Best guide (lowest distance) is always visible.  Guides 2-N are placed
        inside a Bootstrap collapse toggled by a button — no server round-trip.
        Diagram content comes from the pre-computed _guide_content() cache.
        """
        site_rows = site_rows.sort_values("distance")
        rows_list = list(site_rows.iterrows())
        n_extra   = len(rows_list) - 1

        def _one_panel(row, i):
            gid    = str(row["guide_id"])
            cid    = "sel_{}_{}".format(rid, _safe_id(gid))
            dist   = int(row["distance"])
            spacer = str(row["spacer"]).upper()
            best   = i == 0

            cached = content.get(gid, {})
            diag_text = cached.get("diagram", "")
            if diag_text:
                diagram_content = ui.div(
                    ui.HTML("<pre>" + diag_text + "</pre>"), class_="guide-diagram"
                )
            else:
                diagram_content = ui.div("Diagram unavailable — check arm length.", class_="ts-error")

            # warn when an arm was clipped at the sequence boundary
            _req = cached.get("arm_len", 0)
            _short = []
            if _req and cached.get("left_bp", _req) < _req:
                _short.append("5′ ({} bp)".format(cached["left_bp"]))
            if _req and cached.get("right_bp", _req) < _req:
                _short.append("3′ ({} bp)".format(cached["right_bp"]))
            arm_clip_warning = (
                ui.div(
                    "⚠ Short arm at sequence boundary: {} arm(s) truncated "
                    "(requested {} bp).".format(" and ".join(_short), _req),
                    class_="ts-warn",
                ) if _short else ui.span()
            )

            recut_warning = (
                ui.div(
                    "⚠ Recut risk: the insert sequence reconstitutes the guide target "
                    "at the arm/insert junction. Consider a different guide or modifying "
                    "the insert sequence.",
                    class_="ts-warn",
                ) if recut_risk.get(gid, False) else ui.span()
            )

            # per-guide plasmid-specific warnings (shown inline in the guide card)
            plasmid_warnings = []
            try:
                repair = input.repair_type()
            except Exception:
                repair = ""
            if repair == "plasmid":
                # 3'C guide expression warning
                if spacer.endswith("C"):
                    plasmid_warnings.append(
                        "guide contains 3′C and may be inefficient when expressed from a plasmid"
                    )
                # SapI domestication warning for spacer and arms
                for seq_label, seq_val in [
                    ("spacer", spacer),
                    ("5′ HA", cached.get("left", "")),
                    ("3′ HA", cached.get("right", "")),
                ]:
                    if seq_val and find_internal_sapi(seq_val):
                        plasmid_warnings.append(
                            "internal SapI site in {} — must domesticate before assembly".format(seq_label)
                        )
            plasmid_warning_div = (
                ui.div(
                    *[ui.div("⚠ " + msg, class_="ts-warn") for msg in plasmid_warnings]
                ) if plasmid_warnings else ui.span()
            )

            # arm sequence display with SapI sites highlighted (SapTrap mode only)
            arm_sapi_display = ui.span()
            if repair == "plasmid":
                arm_divs = []
                for arm_label, arm_seq, is_left in [
                    ("5′ HA", cached.get("left", ""), True),
                    ("3′ HA", cached.get("right", ""), False),
                ]:
                    if arm_seq and find_internal_sapi(arm_seq):
                        short = len(arm_seq) <= OLIGO_MAX
                        note = "(auto-domesticated in oligo)" if short else "(include mutation in inner primer)"
                        arm_divs.append(ui.div(
                            ui.HTML(_arm_sapi_html(arm_seq, arm_label)),
                            ui.span(" " + note, style="color:#888;font-size:0.85em"),
                        ))
                if arm_divs:
                    arm_sapi_display = ui.div(
                        *arm_divs,
                        style="font-size:0.85em;margin:0.25rem 0;",
                    )

            meta_cells = []
            for label, val in [
                ("Spacer", str(row["spacer"])),
                ("PAM", str(row["pam_seq"])),
                ("Strand", str(row["guide_strand"])),
                ("Distance (bp)", str(dist)),
                ("Recut block", str(row["recut_block_method"])),
                ("Mutation", str(row["mutation_desc"]) or "—"),
            ]:
                meta_cells.append(ui.div(label, class_="param-label"))
                meta_cells.append(ui.div(val, class_="param-value"))

            return ui.div(
                ui.div(
                    ui.input_checkbox(cid, "Use this guide", value=best),
                    ui.span("Guide {}".format(i + 1), class_="fw-semibold"),
                    ui.span("{} bp from cut to insert".format(dist), class_="dist-badge"),
                    plasmid_warning_div,
                    class_="guide-header",
                ),
                arm_clip_warning,
                recut_warning,
                arm_sapi_display,
                diagram_content,
                ui.div(*meta_cells, class_="param-grid mt-2"),
                class_="guide-panel{}".format(" guide-best" if best else ""),
            )

        # Best guide always visible
        _, best_row = rows_list[0]
        best_panel  = _one_panel(best_row, 0)

        divs = [best_panel]

        # Extra guides in a Bootstrap collapse
        if n_extra > 0:
            collapse_id = "guide-alt-{}".format(rid)
            extra_panels = [
                _one_panel(row, i + 1)
                for i, (_, row) in enumerate(rows_list[1:])
            ]
            toggle_btn = ui.tags.button(
                "{} other guide{} ▾".format(n_extra, "s" if n_extra > 1 else ""),
                class_="btn btn-link btn-sm guide-toggle-btn p-0",
                **{"data-bs-toggle": "collapse",
                   "data-bs-target": "#{}".format(collapse_id),
                   "aria-expanded": "false",
                   "aria-controls": collapse_id},
            )
            collapse_div = ui.div(
                *extra_panels,
                id=collapse_id,
                class_="collapse",
            )
            divs.extend([toggle_btn, collapse_div])

        return divs

    def _download_buttons(loc):
        """Single download button per location; produces a ZIP of all reagent files."""
        return ui.div(
            ui.download_button(
                "download_sequences_{}".format(loc),
                "⬇ Download reagents (.zip)",
                class_="btn-primary btn-sm",
            ),
        )

    def _selected_rows():
        """Yield (row, left_arm, right_arm) for each selected guide."""
        df  = reagents_df.get()
        sel = selected_guides.get()
        if df is None or not sel:
            return
        try:
            arm_len = int(input.arm_length())
        except Exception:
            return
        for rid, gid_set in sel.items():
            for gid in gid_set:
                matches = df[(df["residue_index"] == rid) & (df["guide_id"] == gid)]
                if matches.empty:
                    continue
                row = matches.iloc[0]
                try:
                    left, right = truncate_arms(
                        str(row["left_arm"]), str(row["right_arm"]), arm_len
                    )
                except ValueError:
                    continue
                yield row, left, right

    def _guide_label(row):
        """Compact guide label: {site}{strand}{distance} e.g. K45+3."""
        return '{}{}{}{}'.format(
            row['amino_acid'], int(row['residue_index']),
            str(row['guide_strand']), int(row['distance']),
        )

    def _insert_tag_name():
        """Tag name for use in oligo filenames; falls back to 'custom'."""
        try:
            name = input.tag_name()
        except Exception:
            name = 'Custom'
        return 'custom' if name == 'Custom' else name.replace(' ', '-')

    def _assemble_ha_tsv():
        """[runname]_homology-arms.tsv: name, 5p_HA, 3p_HA, guide."""
        rn = run_name.get() or 'run'
        lines = ["name\t5p_HA\t3p_HA\tguide"]
        for row, left, right in _selected_rows():
            site = '{}{}'.format(row['amino_acid'], int(row['residue_index']))
            lines.append('{}_{}\t{}\t{}\t{}'.format(
                rn, site, left, right, str(row['spacer'])))
        return '\n'.join(lines) + '\n'

    def _assemble_guides_tsv():
        """[runname]_guides.tsv: name (run+site+strand+dist), guide sequence."""
        rn = run_name.get() or 'run'
        lines = ["name\tsequence"]
        for row, left, right in _selected_rows():
            lines.append('{}_{}\t{}'.format(rn, _guide_label(row), str(row['spacer'])))
        return '\n'.join(lines) + '\n'

    def _assemble_oligos_tsv():
        """[runname]_oligos.tsv for ssodn/amplicon; None when not applicable."""
        repair = input.repair_type()
        if repair not in ('ssodn', 'amplicon'):
            return None
        rn     = run_name.get() or 'run'
        insert = _insert_seq()
        tag    = _insert_tag_name()
        lines  = ["name\tsequence"]
        for row, left, right in _selected_rows():
            base = '{}_{}_{}'.format(rn, _guide_label(row), tag)
            cut  = int(row['cut_pos'])
            ins  = int(row['insert_pos'])
            if repair == 'ssodn':
                seq, _, _ = build_ssodn(left, right, insert, cut, ins, True)
                lines.append('{}_repair\t{}'.format(base, seq))
            elif repair == 'amplicon':
                # prefer text area (may have manual edits), fall back to file-loaded value
                template = ''
                try:
                    template = input.pcr_template().strip().upper()
                except Exception:
                    pass
                if not template:
                    template = _template_seq.get().strip().upper()
                tm = 60.0
                try:
                    tm = float(input.pcr_tm())
                except Exception:
                    pass
                phos = False
                try:
                    phos = input.pcr_phos()
                except Exception:
                    pass
                if template:
                    (_, fwd), (_, rev) = design_pcr_primers(
                        left, right, template, tm, phos, cut, ins
                    )
                    fwd_ext = ''
                    rev_ext = ''
                    try:
                        fwd_ext = input.pcr_fwd_ext().strip().upper()
                    except Exception:
                        pass
                    try:
                        rev_ext = input.pcr_rev_ext().strip().upper()
                    except Exception:
                        pass
                    if fwd_ext:
                        fwd = fwd_ext + fwd
                    if rev_ext:
                        rev = rev_ext + rev
                    lines.append('{}_F\t{}'.format(base, fwd))
                    lines.append('{}_R\t{}'.format(base, rev))
        return '\n'.join(lines) + '\n' if len(lines) > 1 else None

    def _assemble_saptrap_tsv():
        """[runname]_saptrap-oligos.tsv: sgRNA oligos + 5'HA + 3'HA oligos/primers.

        Returns (tsv_text, genbank_records) where genbank_records is a list of
        (filename, SeqRecord) pairs (empty list when GenBank option is off).
        """
        rn = run_name.get() or 'run'
        tm = 60.0
        try:
            tm = float(input.saptrap_tm())
        except Exception:
            pass
        do_gb = False
        try:
            do_gb = bool(input.saptrap_genbank())
        except Exception:
            pass

        lines = ["name\tsequence\tkind\tnotes"]
        gb_pairs = []   # list of (filename, SeqRecord)
        for row, left, right in _selected_rows():
            oligo_rows, records = build_saptrap_site(
                row, left, right, tm_target=tm, genbank=do_gb
            )
            for orow in oligo_rows:
                lines.append("{}\t{}\t{}\t{}".format(
                    orow["name"], orow["sequence"],
                    orow.get("kind", ""), orow.get("notes", "")))
            for rec in records:
                gb_pairs.append((
                    "{}_{}_saptrap.gb".format(rn, rec.id), rec
                ))

        tsv_text = '\n'.join(lines) + '\n' if len(lines) > 1 else None
        return tsv_text, gb_pairs

    def _build_zip():
        """Produce ZIP bytes: always HA + guides TSVs; oligos when applicable."""
        import zipfile as _zf
        from Bio import SeqIO
        repair = input.repair_type()
        rn  = run_name.get() or 'run'
        buf = io.BytesIO()
        with _zf.ZipFile(buf, 'w', _zf.ZIP_DEFLATED) as zf:
            zf.writestr('{}_homology-arms.tsv'.format(rn), _assemble_ha_tsv())
            zf.writestr('{}_guides.tsv'.format(rn),        _assemble_guides_tsv())
            if repair == 'plasmid':
                tsv_text, gb_pairs = _assemble_saptrap_tsv()
                if tsv_text:
                    zf.writestr('{}_saptrap-oligos.tsv'.format(rn), tsv_text)
                for fname, rec in gb_pairs:
                    gb_buf = io.StringIO()
                    SeqIO.write(rec, gb_buf, "genbank")
                    zf.writestr(fname, gb_buf.getvalue())
            else:
                oligos = _assemble_oligos_tsv()
                if oligos:
                    zf.writestr('{}_oligos.tsv'.format(rn), oligos)
        buf.seek(0)
        return buf.read()
