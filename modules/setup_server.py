from shiny import ui, reactive, render, module
from pathlib import Path
import json, os, shutil

from config import DEFAULT_JSON, TASK_PARAMETERS, SELECTABLE_TASKS, EXCLUDE_ARGS, GLOBAL_TOOLTIPS
from scripts.site_selection_util import get_sequence, save_fasta

# paths to bundled data dirs — relative to this file, not cwd
_ROOT = Path(__file__).parent.parent
_PARAMS_DIR = _ROOT / "params"
_TABLES_DIR = _ROOT / "tables"


def label_with_tip(text, tip=""):
    """Return a label element: plain text when no tip, ⓘ icon + tooltip when tip is provided."""
    if not tip:
        return text
    return ui.span(
        text,
        ui.tooltip(
            ui.span("ⓘ", class_="tip-icon"),
            tip,
            placement="top",
        ),
    )


@module.server
def setup_server(input, output, session, shared_json):

    selected_tasks = reactive.Value([])  # list of task dicts added by the user
    task_values = reactive.Value({})     # snapshot of user-entered param values, keyed by task id

    # reload INPUT_JSON fresh per session to avoid cross-session mutation of defaults
    INPUT_JSON = json.load(open(DEFAULT_JSON, "r"))

    def populate_default_params():
        """Populate the Load Defaults dropdown from params/ json files."""
        params_files = [j.name.removesuffix(".json") for j in _PARAMS_DIR.glob("*.json")]
        ui.update_selectize("load_default_tasks", choices=params_files)
        ui.update_selectize("load_default_tasks", selected="")

    session.on_flush(populate_default_params, once=True)

    # ── sidebar global inputs (rendered so tooltips can be dynamic) ───────────

    @render.ui
    def sidebar_global_inputs():
        """Render sidebar global inputs with ⓘ tooltip icons."""
        t = GLOBAL_TOOLTIPS
        return ui.div(
            ui.input_text("email", label_with_tip("Email Address", t.get("email", "")),
                          value=INPUT_JSON["global"]["email"]),
            ui.input_file("input_file",
                          label_with_tip("Protein sequence (.fa, .fasta, or .pdb)",
                                         t.get("input_file", "")),
                          accept=[".fa", ".fasta", ".pdb"]),
            ui.input_file("input_genomic",
                          label_with_tip("Genomic region FASTA (optional — needed for CRISPR reagents)",
                                         t.get("input_genomic", "")),
                          accept=[".fasta", ".fa"]),
            ui.tags.small(
                ui.tags.em("⚠ Upload a genomic FASTA to enable CRISPR reagent design."),
                style="color:#856404; display:block; margin-top:-8px; margin-bottom:8px;",
            ),
            ui.input_text("run_name",
                          label_with_tip("Analysis name", t.get("run_name", ""))),
            ui.input_text("working_dir",
                          label_with_tip("Output directory", t.get("working_dir", "")),
                          width="100%"),
        )

    # ── button state management ───────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.run_name)
    def update_working_dir():
        """Auto-fill working_dir whenever the run name changes."""
        base = INPUT_JSON["global"]["working_dir"] or os.getcwd()
        if input.run_name():
            ui.update_text("working_dir", value=f"{base}/{input.run_name()}/")

    @reactive.effect
    def update_load_defaults_button_state():
        """Enable the Load button only when a default set is selected."""
        ui.update_action_button("load_defaults_button",
                                disabled=not bool(input.load_default_tasks()))

    @reactive.effect
    @reactive.event(input.task_desc_name)
    def update_task_name_button_state():
        """Enable the Add Task button only when a task name is typed."""
        ui.update_action_button("add_task", disabled=not bool(input.task_desc_name()))

    @reactive.calc
    def requirements_filled():
        """True when all fields required to save an analysis are present."""
        return (
            bool(input.email())
            and bool(input.input_file())
            and bool(input.working_dir())
            and len(selected_tasks()) > 0
        )

    @reactive.effect
    @reactive.event(requirements_filled)
    def update_save_analysis_button_state():
        """Enable the Save Analysis button when all requirements are met."""
        ui.update_action_button("save_analysis", disabled=not requirements_filled())

    @render.text
    def tip_save_analysis():
        """Status hint displayed next to the Save Analysis button."""
        if not requirements_filled():
            return "Requires email, sequence, analysis name, and at least one task."
        if input.save_analysis() == 0:
            return "Click to save analysis"
        return f"Saved to {input.working_dir()}{input.run_name()}.json"

    # ── task param rendering ──────────────────────────────────────────────────

    def build_task_params_input(task):
        """Build the list of UI inputs for one task's configurable parameters."""
        param_list = []
        tips = task.get("tooltips", {})
        for p in task["params"]:
            if p in EXCLUDE_ARGS:
                continue
            lbl = label_with_tip(p.replace("_", " "), tips.get(p, ""))
            if isinstance(task["params"][p], list):
                # list-valued param → select dropdown
                current = task_values().get(task["id"], {}).get(p)
                selected = current if current is not None else task["params"][p][0]
                param_list.append(ui.input_select(
                    f"{task['id']}_{p}", label=lbl,
                    selected=selected, choices=task["params"][p], size=1,
                ))
            else:
                # scalar param → free text
                param_list.append(ui.input_text(
                    f"{task['id']}_{p}",
                    value=task_values().get(task["id"], {}).get(p, task["params"][p]),
                    label=lbl,
                ))
        return param_list

    @render.ui
    def task_params_container():
        """Render all task configuration cards dynamically."""
        tasks = selected_tasks()
        if not tasks:
            return ui.div(
                ui.tags.em("No analysis tasks added yet. Use the form above to add one.",
                           style="color:#6c757d;"),
            )

        cards = []
        for task in tasks:
            label = task["id"].rsplit("_", 1)[0]  # strip trailing analysis type
            remove_btn = ui.input_action_button(
                f"{task['id']}_remove", "✕ Remove",
                class_="btn-sm btn-outline-danger",
            )
            cards.append(
                ui.card(
                    ui.card_header(
                        ui.span(f"{label}"),
                        ui.tags.small(f" — {task['name']}", style="color:#6c757d; font-weight:400;"),
                        remove_btn,
                    ),
                    ui.card_body(
                        ui.layout_column_wrap(
                            *build_task_params_input(task),
                            width=1/2,
                        ),
                    ),
                    class_="task-card",
                )
            )
        return ui.div(*cards)

    # ── task remove handler (dynamic, one observer per task) ──────────────────

    def _make_remove_handler(task_id):
        """Register a reactive effect that removes this task when its button is clicked."""
        @reactive.effect
        @reactive.event(lambda: input[f"{task_id}_remove"]())
        def _handler():
            save_task_values()
            remaining = [t for t in selected_tasks() if t["id"] != task_id]
            selected_tasks.set(remaining)

    # ── task loading ──────────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.load_defaults_button)
    def add_tasks_from_default():
        """Load a saved task list from a params/ JSON file."""
        default_name = input.load_default_tasks()
        if not default_name:
            return
        path = _PARAMS_DIR / f"{default_name}.json"
        default_task_json = json.load(open(path, "r"))
        new_tasks = [
            {"id": tid,
             "name": default_task_json[tid]["analysis"],
             "params": default_task_json[tid]["args"],
             "tooltips": TASK_PARAMETERS.get(default_task_json[tid]["analysis"], {}).get("tooltips", {})}
            for tid in default_task_json
        ]
        save_task_values()
        selected_tasks.set(selected_tasks() + new_tasks)
        for t in new_tasks:
            _make_remove_handler(t["id"])

    @reactive.effect
    def update_save_default_tasks_button():
        """Enable the Save as Default button when a name and tasks both exist."""
        enabled = bool(input.new_default_name()) and len(selected_tasks()) > 0
        ui.update_action_button("save_default_tasks_button", disabled=not enabled)

    @reactive.effect
    @reactive.event(input.save_default_tasks_button)
    def make_new_default():
        """Serialize the current task list to a new defaults JSON file."""
        save_task_values()
        path = _PARAMS_DIR / f"{input.new_default_name()}.json"
        def_json = {}
        for task in selected_tasks():
            def_json = {**def_json, **write_task_json(task)}
        with open(path, "w") as f:
            json.dump(def_json, f, indent=4)

    def save_task_values():
        """Snapshot current input values for all tasks into task_values reactive."""
        snapshot = {}
        for task in selected_tasks():
            task_params = {}
            for param in task["params"]:
                input_id = f"{task['id']}_{param}"
                if input_id in input:
                    task_params[param] = input[input_id]()
            snapshot[task["id"]] = task_params
        task_values.set(snapshot)

    @reactive.effect
    @reactive.event(input.add_task)
    def add_task():
        """Add a new task row for the selected analysis type."""
        task_name = input.task_selector()
        if not task_name:
            return
        save_task_values()
        task_id = f"{input.task_desc_name()}_{task_name}"
        default_params = INPUT_JSON["analyses"][task_name]["args"]
        task_params = {p: default_params[p] for p in default_params if p not in EXCLUDE_ARGS}
        default_tips = INPUT_JSON["analyses"][task_name]["tooltips"]
        task_tooltips = {p: default_tips[p] for p in default_tips if p not in EXCLUDE_ARGS}
        new_task = {
            "id": task_id, "name": task_name,
            "params": task_params, "tooltips": task_tooltips,
        }
        selected_tasks.set(selected_tasks() + [new_task])
        _make_remove_handler(task_id)
        ui.update_text("task_desc_name", value="")

    # ── serialisation helpers ─────────────────────────────────────────────────

    def write_task_json(task, global_params=None):
        """Serialize one task's current parameter inputs to a dict."""
        task_json = {task["id"]: {"analysis": task["name"], "args": {}}}
        for p in TASK_PARAMETERS[task["name"]]["args"]:
            input_id = f"{task['id']}_{p}"
            if input[input_id].is_set():
                value = input[input_id]()
            else:
                value = TASK_PARAMETERS[task["name"]]["args"][p]
            task_json[task["id"]]["args"][p] = value
        if global_params is not None:
            task_json[task["id"]]["args"] = {**task_json[task["id"]]["args"], **global_params}
        # build output path from working dir + run name + task id + default extension
        out_suffix = TASK_PARAMETERS[task["name"]]["args"]["output"]
        task_json[task["id"]]["args"]["output"] = (
            f"{input.working_dir()}{input.run_name()}.{task['id']}{out_suffix}"
        )
        return task_json

    # ── save analysis ─────────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.save_analysis)
    def save_analysis():
        """Write all parameters to a run JSON and publish the path to shared state."""
        # create working directory
        Path(input.working_dir()).mkdir(parents=True, exist_ok=True)

        # update global parameters in this session's copy of INPUT_JSON
        INPUT_JSON["global"]["email"] = input.email()
        INPUT_JSON["global"]["run_name"] = input.run_name()
        INPUT_JSON["global"]["working_dir"] = input.working_dir()

        # copy uploaded input file into working directory with a clean name
        tmp_input = Path(input.input_file()[0]["datapath"])
        new_input_filename = input.working_dir() + input.run_name()
        new_input_filename += ".pdb" if tmp_input.suffix == ".pdb" else ".fa"
        INPUT_JSON["global"]["input_file"] = new_input_filename
        shutil.copy(tmp_input, new_input_filename)

        # copy optional genomic FASTA if provided; track its path for reagents auto-inject
        genomic_path = ""
        if input.input_genomic():
            tmp_genomic = Path(input.input_genomic()[0]["datapath"])
            new_genomic = (
                input.working_dir() + input.run_name() + ".genomic" + tmp_genomic.suffix
            )
            INPUT_JSON["global"]["genomic_file"] = new_genomic
            shutil.copy(tmp_genomic, new_genomic)
            genomic_path = new_genomic

        # for PDB input, extract FASTA and point input_file at it
        if tmp_input.suffix == ".pdb":
            INPUT_JSON["global"]["pdb"] = INPUT_JSON["global"]["input_file"]
            new_fasta = INPUT_JSON["global"]["pdb"].replace(".pdb", ".fa")
            save_fasta(
                Path(INPUT_JSON["global"]["pdb"]).stem,
                get_sequence(INPUT_JSON["global"]["pdb"]),
                new_fasta,
            )
            INPUT_JSON["global"]["input_file"] = new_fasta

        # assemble and write the run JSON
        out_json_name = f"{input.working_dir()}{input.run_name()}.json"
        out_json = {
            "scripts": INPUT_JSON["scripts"],
            "global": INPUT_JSON["global"],
        }
        for task in selected_tasks():
            # plddt needs to know whether we have a PDB vs AFDB lookup
            if task["name"] == "plddt":
                task["params"]["existing_AF2"] = 0 if INPUT_JSON["global"].get("pdb") else 1
            out_json = {**out_json, **write_task_json(task, INPUT_JSON["global"])}

        # auto-inject reagents task when a genomic FASTA was uploaded
        if genomic_path:
            rea_args = dict(TASK_PARAMETERS["reagents"]["args"])
            rea_args["genomic_fasta"] = genomic_path
            rea_args["genewise"] = ""   # task_runner will run Genewise automatically
            rea_args["email"] = input.email()
            rea_args["working_dir"] = input.working_dir()
            rea_args["run_name"] = input.run_name()
            rea_args["input_file"] = INPUT_JSON["global"]["input_file"]
            out_suffix = TASK_PARAMETERS["reagents"]["args"]["output"]
            rea_args["output"] = f"{input.working_dir()}{input.run_name()}.reagents{out_suffix}"
            out_json["REAGENTS_reagents"] = {"analysis": "reagents", "args": rea_args}
        else:
            ui.notification_show(
                "No genomic FASTA uploaded — CRISPR reagent design will not run. "
                "Upload a genomic FASTA in the sidebar and re-save to enable it.",
                type="warning",
                duration=8,
            )

        with open(out_json_name, "w") as f:
            json.dump(out_json, f, indent=4)

        shared_json.set(out_json_name)
