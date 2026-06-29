"""setup_server.py — Shiny module server for the Setup tab.

Thin reactive layer: reads user inputs, delegates data construction to
setup_logic.py, writes the run JSON, and publishes its path via shared_json.
"""

import json
import os
import requests
import shutil
from pathlib import Path

from shiny import module, reactive, render, ui

from config import (
    DEFAULT_SPECIES, GLOBAL_TOOLTIPS, TASK_PARAMETERS,
    task_hidden, task_defaults, task_choices, task_tooltips, task_output_suffix,
    GLOBAL_DEFAULTS,
)
from modules.setup_ui import label_with_tip, TASK_DESCRIPTIONS
from modules.setup_logic import (
    build_defaults_entry,
    build_global_block,
    build_reagents_entry,
    build_run_json,
    build_task_entry,
    make_task,
)
from scripts.site_selection_util import get_sequence, save_fasta, extract_bfactors_from_pdb

_ROOT = Path(__file__).parent.parent
_PARAMS_DIR = _ROOT / "params"
_TABLES_DIR = _ROOT / "tables"


def _table_choices(ext="*.tsv"):
    """Return {absolute_path_str: display_name} for files in tables/ matching ext."""
    return {str(f): f.stem for f in sorted(_TABLES_DIR.glob(ext))}


# ── widget builders ───────────────────────────────────────────────────────────

def _make_param_widget(widget_id, param, value, tip, choices=None):
    """Return one input widget for a task parameter."""
    lbl = label_with_tip(param.replace("_", " "), tip)
    if param == "scores_file":
        file_choices = _table_choices("*.tsv")
        selected = value if value in file_choices else next(iter(file_choices), "")
        return ui.input_select(widget_id, label=lbl, choices=file_choices, selected=selected)
    if choices:
        selected = value if value in choices else choices[0]
        return ui.input_select(widget_id, label=lbl, choices=choices, selected=selected, size=1)
    return ui.input_text(widget_id, label=lbl, value=str(value) if value != "" else "")


def _build_param_inputs(task, task_values_snap):
    """Build the list of param widgets for one task card."""
    tips       = task.get("tooltips", {})
    all_choices = task.get("choices", {})
    hidden     = task_hidden(task["name"])
    widgets = []
    for p, default in task["params"].items():
        if p in hidden:
            continue
        widget_id = f"{task['id']}_{p}"
        # use the snapshotted value if present (preserves edits across re-renders)
        current = task_values_snap.get(task["id"], {}).get(p, default)
        widgets.append(_make_param_widget(widget_id, p, current, tips.get(p, ""), all_choices.get(p)))
    return widgets


# ── module server ─────────────────────────────────────────────────────────────

@module.server
def setup_server(input, output, session, shared_json):

    tasks      = reactive.Value([])   # list of in-memory task dicts
    task_snap  = reactive.Value({})   # {task_id: {param: value}} — preserved across re-renders

    organism_taxid   = reactive.Value("")   # resolved taxid string from organism selection
    _pdb_plddt_valid = reactive.Value(True) # False when uploaded PDB lacks valid pLDDT B-factors
    _search_hits     = reactive.Value([])   # [(name, taxid_str)] from UniProt taxonomy search

    # ── organism selection ────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.organism)
    def _organism_changed():
        """Resolve taxid from preset dropdown; pre-fill existing blast task widgets."""
        name = input.organism()
        if not name or name not in DEFAULT_SPECIES:
            return
        taxid = DEFAULT_SPECIES[name]
        if taxid is None:   # "Other (search...)" sentinel — wait for search result
            _search_hits.set([])
            return
        organism_taxid.set(str(taxid))

    @render.ui
    def organism_search_ui():
        """Show a search row only when 'Other (search…)' is selected."""
        if input.organism() != "Other (search...)":
            return None
        hits = _search_hits()
        result_widget = None
        if hits:
            choices = {"": "— select —",
                       **{str(taxid): f"{name} ({taxid})" for name, taxid in hits}}
            result_widget = ui.div(
                ui.input_select("organism_search_result", "", choices=choices),
                style="margin-top: 0.2rem;",
            )
        return ui.div(
            ui.div(
                ui.input_text("organism_search_text", "",
                              placeholder="Type species name…"),
                ui.input_action_button("organism_search_btn", "Search",
                                       class_="btn-sm btn-outline-secondary"),
                style="display:flex; gap:0.4rem; align-items:flex-end;",
            ),
            *([result_widget] if result_widget else []),
            style="margin-top: 0.2rem;",
        )

    @reactive.effect
    @reactive.event(input.organism_search_btn)
    def _organism_search():
        """Query UniProt taxonomy API and populate the result select."""
        query = input.organism_search_text().strip() if "organism_search_text" in input else ""
        if not query:
            return
        try:
            resp = requests.get(
                "https://rest.uniprot.org/taxonomy/search",
                params={"query": query, "format": "json", "size": 5},
                timeout=10,
            )
            resp.raise_for_status()
            results = resp.json().get("results", [])
            hits = [(r["scientificName"], r["taxonId"])
                    for r in results if "taxonId" in r]
            _search_hits.set(hits)
        except Exception as e:
            ui.notification_show(f"Taxonomy search failed: {e}", type="error", duration=5)

    @reactive.effect
    def _organism_search_result_changed():
        """Set taxid when user picks from search results."""
        if "organism_search_result" not in input:
            return
        val = input.organism_search_result()
        if not val:
            return
        organism_taxid.set(val)

    # ── taxonomic lineage ─────────────────────────────────────────────────────

    # ranks worth showing as BLAST scope options (clades and "no rank" filtered out)
    _BLAST_RANKS = {
        "domain", "superkingdom", "kingdom", "subkingdom",
        "phylum", "subphylum", "class", "subclass",
        "order", "suborder", "infraorder", "superfamily",
        "family", "subfamily", "genus", "species",
    }

    _lineage = reactive.Value([])  # [{"rank": str, "name": str, "taxid": str}, ...]

    @reactive.effect
    def _fetch_lineage():
        """Fetch taxonomic lineage from UniProt whenever the resolved taxid changes."""
        tid = organism_taxid()
        if not tid:
            _lineage.set([])
            return
        try:
            resp = requests.get(
                f"https://rest.uniprot.org/taxonomy/{tid}",
                params={"format": "json"},
                timeout=10,
            )
            resp.raise_for_status()
            data = resp.json()
            entries = [
                {"rank": e["rank"], "name": e["scientificName"], "taxid": str(e["taxonId"])}
                for e in data.get("lineage", [])
                if e.get("rank") in _BLAST_RANKS
            ]
            # append the organism itself
            entries.append({
                "rank": data.get("rank", "species"),
                "name": data.get("scientificName", ""),
                "taxid": tid,
            })
            _lineage.set(entries)
        except Exception as e:
            ui.notification_show(f"Lineage fetch failed: {e}", type="warning", duration=4)
            _lineage.set([])

    @render.ui
    def lineage_ui():
        """Display taxonomic lineage as a reference table (rank, name, taxid)."""
        entries = _lineage()
        if not entries:
            return None
        rows = [
            ui.tags.tr(
                ui.tags.td(e["rank"],   style="padding:1px 6px; color:#6c757d;"),
                ui.tags.td(ui.tags.em(e["name"]), style="padding:1px 6px;"),
                ui.tags.td(e["taxid"],  style="padding:1px 8px; font-family:monospace;"),
            )
            for e in entries
        ]
        return ui.div(
            ui.tags.p("Use taxids to limit BLAST searches to specific taxa:",
                      style="font-size:0.78rem; color:#495057; margin:0.4rem 0 0.15rem;"),
            ui.tags.table(
                ui.tags.tbody(*rows),
                style="font-size:0.75rem; border-collapse:collapse;",
            ),
        )

    # ── PDB pLDDT validity check ──────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.input_file)
    def _check_input_file():
        """Warn if the uploaded PDB has no valid AlphaFold pLDDT B-factors."""
        files = input.input_file()
        if not files:
            return
        if not files[0]["name"].lower().endswith(".pdb"):
            _pdb_plddt_valid.set(True)
            return
        try:
            bfs = extract_bfactors_from_pdb(files[0]["datapath"])
            valid = bool(bfs) and all(0 <= v <= 100 for v in bfs) and (max(bfs) - min(bfs)) > 5
        except Exception:
            valid = False
        _pdb_plddt_valid.set(valid)
        if not valid:
            ui.notification_show(
                "The uploaded PDB does not appear to contain AlphaFold pLDDT scores "
                "(B-factors outside 0–100 or uniform). "
                "The AlphaFold Database will be searched for a matching structure.",
                type="warning",
                duration=10,
            )

    # ── preset refresh (called after saving a new preset) ────────────────────

    def _populate_presets():
        """Refresh the Load preset dropdown after a new preset is saved."""
        names = [p.stem for p in sorted(_PARAMS_DIR.glob("*.json"))]
        choices = {"": "— select —", **{n: n for n in names}}
        ui.update_select("load_preset", choices=choices, selected="")

    # ── button enable / disable ───────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.task_label)
    def _toggle_add():
        """Enable Add only when a task label has been typed."""
        ui.update_action_button("add_task", disabled=not bool(input.task_label()))

    @reactive.effect
    def _toggle_load():
        """Enable Load only when a preset is selected."""
        ui.update_action_button("load_preset_btn",
                                disabled=not bool(input.load_preset()))

    @reactive.effect
    def _toggle_save_preset():
        """Enable Save preset only when a name is typed and tasks exist."""
        ui.update_action_button(
            "save_preset_btn",
            disabled=not (bool(input.preset_name()) and len(tasks()) > 0),
        )

    @reactive.calc
    def _ready():
        """True when all required fields for Save Analysis are present."""
        return (
            bool(input.email())
            and bool(input.input_file())
            and bool(input.working_dir())
            and len(tasks()) > 0
        )

    @reactive.effect
    @reactive.event(_ready)
    def _toggle_save():
        ui.update_action_button("save_analysis", disabled=not _ready())

    # ── auto-fill working directory ───────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.run_name)
    def _autofill_dir():
        """Derive working_dir from run_name when the user types a name."""
        if not input.run_name():
            return
        base = GLOBAL_DEFAULTS.get("working_dir") or "data"
        ui.update_text("working_dir", value=f"{base}/{input.run_name()}/")

    # ── save-status hint ──────────────────────────────────────────────────────

    @render.text
    def save_status():
        """One-line hint next to the Save button."""
        if not _ready():
            return "Requires email, a sequence file, an analysis name, and at least one task."
        if input.save_analysis() == 0:
            return "Ready — click Save Analysis to write the run configuration."
        return f"Saved → {input.working_dir()}{input.run_name()}.json"

    @render.text
    def task_type_desc():
        """One-line description of the currently selected analysis type."""
        return TASK_DESCRIPTIONS.get(input.task_type(), "")

    # ── task state helpers ────────────────────────────────────────────────────

    def _snapshot():
        """Read all current param inputs into task_snap (copy-and-set for reactivity)."""
        snap = {}
        for t in tasks():
            hidden = task_hidden(t["name"])
            snap[t["id"]] = {}
            for p in t["params"]:
                if p in hidden:
                    continue
                wid = f"{t['id']}_{p}"
                if wid in input:
                    snap[t["id"]][p] = input[wid]()
        task_snap.set(snap)

    def _collect_args(task):
        """Read input values for all params in a task, falling back to registry defaults."""
        args = {}
        for p, pcfg in TASK_PARAMETERS[task["name"]]["params"].items():
            wid = f"{task['id']}_{p}"
            # is_set() is False for hidden params (no widget) — fall back to default
            if input[wid].is_set():
                args[p] = input[wid]()
            else:
                args[p] = pcfg["default"]
        return args

    # ── task CRUD ─────────────────────────────────────────────────────────────

    def _register_remove(task_id):
        """Register a one-off reactive effect that removes a task when its button fires."""
        @reactive.effect
        @reactive.event(lambda: input[f"{task_id}_remove"]())
        def _handler():
            _snapshot()
            tasks.set([t for t in tasks() if t["id"] != task_id])

    def _register_table_upload(task_id):
        """Copy an uploaded TSV to /tables/, update snap, and re-render the card."""
        @reactive.effect
        @reactive.event(lambda: input[f"{task_id}_add_table"]())
        def _handler():
            files = input[f"{task_id}_add_table"]()
            if not files:
                return
            f = files[0]
            dest = _TABLES_DIR / f["name"]
            shutil.copy(f["datapath"], dest)
            snap = dict(task_snap())
            snap.setdefault(task_id, {})["scores_file"] = str(dest)
            task_snap.set(snap)
            tasks.set(list(tasks()))
            ui.notification_show(f"Table '{f['name']}' added.", type="message", duration=3)

    @reactive.effect
    @reactive.event(input.add_task)
    def _add_task():
        """Append a new task card for the selected analysis type."""
        ttype = input.task_type()
        label = input.task_label().strip()
        if not ttype or not label:
            return
        _snapshot()
        hidden = task_hidden(ttype)
        params  = {p: v for p, v in task_defaults(ttype).items() if p not in hidden}
        tips    = task_tooltips(ttype)
        choices = task_choices(ttype)
        task = make_task(ttype, label, params, tips, choices)
        task["start_open"] = True   # newly added tasks open by default
        tasks.set(tasks() + [task])
        _register_remove(task["id"])
        if ttype == "scores":
            _register_table_upload(task["id"])
        ui.update_text("task_label", value="")

    # ── preset load / save ────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.load_preset_btn)
    def _load_preset():
        """Load a saved task set from params/ and append to the current list."""
        name = input.load_preset()
        if not name:
            return
        path = _PARAMS_DIR / f"{name}.json"
        try:
            data = json.load(open(path))
        except FileNotFoundError:
            ui.notification_show(f"Preset '{name}' not found.", type="error")
            return
        _snapshot()
        new_tasks = []
        for tid, entry in data.items():
            ttype  = entry.get("type") or entry.get("analysis", "")  # accept old "analysis" key
            hidden = task_hidden(ttype)
            # use saved args as params, supplement tooltips + choices from registry
            params  = {p: v for p, v in entry["args"].items() if p not in hidden and p != "output"}
            tips    = task_tooltips(ttype)
            choices = task_choices(ttype)
            task = {"id": tid, "name": ttype, "params": params, "tooltips": tips,
                    "choices": choices, "start_open": False}   # preloaded tasks start collapsed
            new_tasks.append(task)
            _register_remove(tid)
            if ttype == "scores":
                _register_table_upload(tid)
        tasks.set(tasks() + new_tasks)

    @reactive.effect
    @reactive.event(input.save_preset_btn)
    def _save_preset():
        """Write the current task set (no global block) to a new params/ file."""
        _snapshot()
        name = input.preset_name().strip()
        if not name:
            return
        result = {}
        wd, rn = input.working_dir(), input.run_name()
        for task in tasks():
            args = _collect_args(task)
            result.update(build_defaults_entry(task["id"], task["name"], args,
                                               task_output_suffix(task["name"]), wd, rn))
        path = _PARAMS_DIR / f"{name}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(result, f, indent=4)
        ui.notification_show(f"Saved preset '{name}'.", type="message", duration=3)
        _populate_presets()

    # ── task card rendering ───────────────────────────────────────────────────

    @render.ui
    def task_cards():
        """Render tasks as a collapsible accordion; new tasks open, preloaded closed."""
        current_tasks = tasks()
        snap = task_snap()

        if not current_tasks:
            return ui.p(
                "No tasks yet — choose a type and label above, then click Add.",
                style="color:#6c757d; margin-top:0.4rem; font-size:0.8rem;",
            )

        open_ids = [t["id"] for t in current_tasks if t.get("start_open", False)]

        panels = []
        for task in current_tasks:
            label = task["id"].rsplit("_", 1)[0]
            widgets = _build_param_inputs(task, snap)
            add_table = (
                ui.div(
                    ui.input_file(f"{task['id']}_add_table", "＋ Add table (.tsv)",
                        accept=[".tsv"]),
                    style="margin-top:0.2rem;",
                )
                if task["name"] == "scores" else None
            )
            panels.append(
                ui.accordion_panel(
                    ui.span(
                        label,
                        ui.tags.span(task["name"], class_="task-type-badge"),
                    ),
                    # params grid
                    ui.layout_column_wrap(*widgets, width=1/2) if widgets else
                    ui.p("No configurable parameters.", style="color:#6c757d; font-size:0.8rem;"),
                    # scores-only: upload a new property table
                    *([add_table] if add_table else []),
                    # remove button at bottom of body
                    ui.div(
                        ui.input_action_button(
                            f"{task['id']}_remove", "✕ Remove task",
                            class_="btn-sm btn-outline-danger",
                        ),
                        style="margin-top:0.6rem;",
                    ),
                    value=task["id"],
                )
            )

        return ui.accordion(
            *panels,
            open=open_ids,
            multiple=True,
            class_="task-accordion",
        )

    # ── save analysis ─────────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.save_analysis)
    def _save_analysis():
        """Write the run JSON and publish its path via shared_json."""
        wd  = input.working_dir()
        rn  = input.run_name()
        email = input.email()

        Path(wd).mkdir(parents=True, exist_ok=True)

        # copy uploaded protein file to the working directory
        tmp_protein = Path(input.input_file()[0]["datapath"])
        ext = ".pdb" if tmp_protein.suffix.lower() == ".pdb" else ".fa"
        protein_dest = wd + rn + ext
        shutil.copy(tmp_protein, protein_dest)

        # if PDB: extract FASTA and set pdb field; otherwise pdb stays empty
        pdb_path = ""
        input_file = protein_dest
        if ext == ".pdb":
            pdb_path = protein_dest
            fasta_dest = protein_dest.replace(".pdb", ".fa")
            save_fasta(Path(pdb_path).stem, get_sequence(pdb_path), fasta_dest)
            input_file = fasta_dest

        # copy optional genomic FASTA
        genomic_path = ""
        if input.input_genomic():
            tmp_gen = Path(input.input_genomic()[0]["datapath"])
            genomic_dest = wd + rn + ".genomic" + tmp_gen.suffix
            shutil.copy(tmp_gen, genomic_dest)
            genomic_path = genomic_dest

        # build the global block
        global_block = build_global_block(
            email        = email,
            run_name     = rn,
            working_dir  = wd,
            input_file   = input_file,
            pdb          = pdb_path,
            genomic_file = genomic_path,
        )

        # build per-task entries
        task_entries = []
        for task in tasks():
            args = _collect_args(task)
            if task["name"] == "plddt":
                # use the supplied PDB only when it contains valid pLDDT B-factors
                use_supplied_pdb = bool(pdb_path) and _pdb_plddt_valid()
                args["existing_AF2"] = 0 if use_supplied_pdb else 1
                # write pdb into task args so parse_run merge doesn't lose it to
                # global_block being overridden by the empty-string task default
                args["pdb"] = pdb_path if use_supplied_pdb else ""
                # always inject the resolved taxid so AFDB lookup can filter by organism
                args["taxid"] = organism_taxid()
            task_entries.append(
                build_task_entry(task["id"], task["name"], args,
                                 task_output_suffix(task["name"]), wd, rn)
            )

        # build reagents entry or warn
        reagents_entry = None
        if genomic_path:
            reagents_entry = build_reagents_entry(
                defaults     = task_defaults("reagents"),
                genomic_path = genomic_path,
                out_suffix   = task_output_suffix("reagents"),
                working_dir  = wd,
                run_name     = rn,
            )
        else:
            ui.notification_show(
                "No genomic FASTA uploaded — CRISPR reagent design will not run. "
                "Upload a genomic region FASTA and re-save to enable it.",
                type="warning",
                duration=8,
            )

        run_json = build_run_json(
            global_block   = global_block,
            task_entries   = task_entries,
            reagents_entry = reagents_entry,
        )

        out_path = f"{wd}{rn}.json"
        with open(out_path, "w") as f:
            json.dump(run_json, f, indent=4)

        shared_json.set(out_path)
