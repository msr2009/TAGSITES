"""setup_server.py — Shiny module server for the Setup tab.

Thin reactive layer: reads user inputs, delegates data construction to
setup_logic.py, writes the run JSON, and publishes its path via shared_json.
"""

import json
import os
import shutil
from pathlib import Path

from shiny import module, reactive, render, ui

from config import DEFAULT_JSON, EXCLUDE_ARGS, GLOBAL_TOOLTIPS, TASK_PARAMETERS
from modules.setup_ui import label_with_tip, TASK_DESCRIPTIONS
from modules.setup_logic import (
    build_defaults_entry,
    build_global_block,
    build_reagents_entry,
    build_run_json,
    build_task_entry,
    make_task,
)
from scripts.site_selection_util import get_sequence, save_fasta

_ROOT = Path(__file__).parent.parent
_PARAMS_DIR = _ROOT / "params"


# ── widget builders ───────────────────────────────────────────────────────────

def _make_param_widget(widget_id, param, value, tip):
    """Return one input widget for a task parameter."""
    lbl = label_with_tip(param.replace("_", " "), tip)
    if isinstance(value, list):
        return ui.input_select(widget_id, label=lbl,
                               choices=value, selected=value[0], size=1)
    return ui.input_text(widget_id, label=lbl, value=str(value) if value != "" else "")


def _build_param_inputs(task, task_values_snap):
    """Build the list of param widgets for one task card."""
    tips = task.get("tooltips", {})
    widgets = []
    for p, default in task["params"].items():
        if p in EXCLUDE_ARGS:
            continue
        widget_id = f"{task['id']}_{p}"
        # use the snapshotted value if present (preserves edits across re-renders)
        current = task_values_snap.get(task["id"], {}).get(p, default)
        widgets.append(_make_param_widget(widget_id, p, current, tips.get(p, "")))
    return widgets


# ── module server ─────────────────────────────────────────────────────────────

@module.server
def setup_server(input, output, session, shared_json):

    # reload defaults fresh per session — avoids cross-session mutation
    INPUT_JSON = json.load(open(DEFAULT_JSON))

    tasks      = reactive.Value([])   # list of in-memory task dicts
    task_snap  = reactive.Value({})   # {task_id: {param: value}} — preserved across re-renders

    # ── preset refresh (called after saving a new preset) ────────────────────

    def _populate_presets():
        """Refresh the Load preset dropdown after a new preset is saved."""
        names = [p.stem for p in sorted(_PARAMS_DIR.glob("*.json"))]
        ui.update_selectize("load_preset", choices=names)
        ui.update_selectize("load_preset", selected="")

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
        base = INPUT_JSON["global"].get("working_dir") or os.getcwd()
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
            snap[t["id"]] = {}
            for p in t["params"]:
                if p in EXCLUDE_ARGS:
                    continue
                wid = f"{t['id']}_{p}"
                if wid in input:
                    snap[t["id"]][p] = input[wid]()
        task_snap.set(snap)

    def _collect_args(task):
        """Read input values for all params in a task, falling back to config defaults."""
        base = TASK_PARAMETERS[task["name"]]["args"]
        args = {}
        for p, default in base.items():
            wid = f"{task['id']}_{p}"
            # is_set() returns False for params excluded from the UI (no widget exists)
            if input[wid].is_set():
                args[p] = input[wid]()
            else:
                args[p] = default
        return args

    # ── task CRUD ─────────────────────────────────────────────────────────────

    def _register_remove(task_id):
        """Register a one-off reactive effect that removes a task when its button fires."""
        @reactive.effect
        @reactive.event(lambda: input[f"{task_id}_remove"]())
        def _handler():
            _snapshot()
            tasks.set([t for t in tasks() if t["id"] != task_id])

    @reactive.effect
    @reactive.event(input.add_task)
    def _add_task():
        """Append a new task card for the selected analysis type."""
        ttype = input.task_type()
        label = input.task_label().strip()
        if not ttype or not label:
            return
        _snapshot()
        cfg = TASK_PARAMETERS[ttype]
        # filter excluded args from the in-memory params (they'll still be written at save time)
        params = {p: v for p, v in cfg["args"].items() if p not in EXCLUDE_ARGS}
        tips   = {p: v for p, v in cfg.get("tooltips", {}).items() if p not in EXCLUDE_ARGS}
        task = make_task(ttype, label, params, tips)
        task["start_open"] = True   # newly added tasks open by default
        tasks.set(tasks() + [task])
        _register_remove(task["id"])
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
            ttype = entry["analysis"]
            cfg = TASK_PARAMETERS.get(ttype, {})
            # use saved args as params, supplement tooltips from config
            params = {p: v for p, v in entry["args"].items() if p not in EXCLUDE_ARGS}
            tips   = {p: v for p, v in cfg.get("tooltips", {}).items() if p not in EXCLUDE_ARGS}
            task = {"id": tid, "name": ttype, "params": params, "tooltips": tips,
                    "start_open": False}   # preloaded tasks start collapsed
            new_tasks.append(task)
            _register_remove(tid)
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
            out_suffix = TASK_PARAMETERS[task["name"]]["args"].get("output", "")
            result.update(build_defaults_entry(task["id"], task["name"], args, out_suffix, wd, rn))
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
            panels.append(
                ui.accordion_panel(
                    ui.span(
                        label,
                        ui.tags.span(task["name"], class_="task-type-badge"),
                    ),
                    # params grid
                    ui.layout_column_wrap(*widgets, width=1/2) if widgets else
                    ui.p("No configurable parameters.", style="color:#6c757d; font-size:0.8rem;"),
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
            base_global  = INPUT_JSON["global"],
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
            # plddt: auto-set existing_AF2 based on whether a PDB was supplied
            if task["name"] == "plddt":
                args["existing_AF2"] = 0 if pdb_path else 1
            out_suffix = TASK_PARAMETERS[task["name"]]["args"].get("output", "")
            task_entries.append(
                build_task_entry(task["id"], task["name"], args, global_block,
                                 out_suffix, wd, rn)
            )

        # build reagents entry or warn
        reagents_entry = None
        if genomic_path:
            reagents_entry = build_reagents_entry(
                defaults_cfg = TASK_PARAMETERS["reagents"]["args"],
                genomic_path = genomic_path,
                email        = email,
                working_dir  = wd,
                run_name     = rn,
                input_file   = input_file,
            )
        else:
            ui.notification_show(
                "No genomic FASTA uploaded — CRISPR reagent design will not run. "
                "Upload a genomic region FASTA and re-save to enable it.",
                type="warning",
                duration=8,
            )

        run_json = build_run_json(
            scripts        = INPUT_JSON["scripts"],
            global_block   = global_block,
            task_entries   = task_entries,
            reagents_entry = reagents_entry,
        )

        out_path = f"{wd}{rn}.json"
        with open(out_path, "w") as f:
            json.dump(run_json, f, indent=4)

        shared_json.set(out_path)
