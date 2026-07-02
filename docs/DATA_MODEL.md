# TAGSITES data model

This document explains the three file types that carry data in TAGSITES, how they
are built, and how to add a new analysis task type.

---

## File types

### 1. `task_definitions.json` — static registry (never written by the app)

The single source of truth for every task type. Edit this by hand when adding or
changing task types. The app reads it at startup via `scripts/task_registry.py`.

**Shape:**

```json
{
  "global_defaults":  { "scripts_folder": "./scripts/" },
  "global_tooltips":  { "email": "...", "run_name": "...", ... },
  "tasks": {
    "<task_type>": {
      "script":        "<filename>.py",
      "result_type":   "continuous" | "range" | "none",
      "output_suffix": ".ext",
      "colors":        ["#hex1", "#hex2", "#hex3"],   (optional; continuous tasks only)
      "selectable":    true | false,                  (false = auto-injected, not user-selectable)
      "companions":    { "<key>": ".ext", ... },       (optional sidecar files)
      "params": {
        "<param>": {
          "default":  <value>,
          "choices":  ["opt1", "opt2"],  (optional → renders a dropdown in the UI)
          "tooltip":  "Help text.",      (optional → shown next to the label)
          "hidden":   true               (optional → not shown in the form; set programmatically)
        }
      }
    }
  }
}
```

**What `config.py` derives from it:**
`AVAILABLE_TASKS`, `SELECTABLE_TASKS`, `RESULTS_TYPE_DICT`, `ANALYSIS_COLORS`,
`GLOBAL_KEYS`, `EXCLUDE_ARGS`, and the accessor helpers (`task_defaults`,
`task_choices`, `task_tooltips`, `task_hidden`, `task_script`,
`task_output_suffix`, `task_companions`, `companion_path`).

---

### 2. `{working_dir}{run_name}.run.json` — per-run state (written by the app / CLI)

Stores only what changes per run. No tooltips, no choices, no script names,
no global keys copied per-task.

**Shape:**

```json
{
  "global": {
    "email":         "you@example.com",
    "working_dir":   "data/my-gene/",
    "run_name":      "my-gene",
    "input_file":    "data/my-gene/my-gene.fa",
    "pdb":           "",
    "scripts_folder": "./scripts/",
    "selected_sites": [],
    "genomic_file":  "data/my-gene/my-gene.genomic.fa"  (optional)
  },
  "tasks": {
    "WORM_blast": {
      "type": "blast",
      "args": {
        "taxid":    "6239",
        "evalue":   "1e-10",
        "max_hits": "100",
        "output":   "data/my-gene/my-gene.WORM_blast.jsd"
      }
    },
    "HUMAN_plddt": {
      "type": "plddt",
      "args": {
        "pdb":         "",
        "existing_AF2": 1,
        "evalue":      "1e-10",
        "output":      "data/my-gene/my-gene.HUMAN_plddt.txt"
      }
    }
  }
}
```

**Written by:** `modules/setup_logic.py` → `build_run_json()`. The only
post-run write is `global.selected_sites` (set by the Results tab).

**Read by:** `modules/progress_logic.py:parse_run()` (merges global into each
task's args at read time so runners receive a complete args dict),
`utils/results.py:load_data_from_json()`, and the CLI orchestrator
`scripts/run_tag_sites_from_json.py`.

---

### 3. `{working_dir}{run_name}.status.json` — job status (written by the app / CLI)

One entry per task, updated live as tasks run. Never written by `setup_logic`.

**Shape:**

```json
{
  "WORM_blast": {
    "status":  "success" | "running" | "pending" | "failed",
    "stage":   "blast",
    "log":     "...",
    "output":  "data/my-gene/my-gene.WORM_blast.jsd",
    "job_id":  "",
    "message": ""
  }
}
```

---

### 4. `params/*.json` — saved task presets (optional)

Same shape as the `tasks` section of the run JSON, minus `global`. Loaded via
the Setup tab's "Load preset" button.

```json
{
  "WORM_blast": { "type": "blast", "args": { "taxid": "6239", "evalue": "1e-10", ... } }
}
```

---

## Adding a new task type

Three steps; nothing else needs to change.

### Step 1 — add a block to `task_definitions.json`

```json
"my_analysis": {
  "script":        "my_analysis.py",
  "result_type":   "continuous",
  "output_suffix": ".tsv",
  "colors":        ["#1a237e", "#283593", "#3949ab"],
  "selectable":    true,
  "params": {
    "threshold": { "default": 0.5, "tooltip": "Score threshold (0–1)." },
    "mode":      { "default": "fast", "choices": ["fast", "sensitive"], "tooltip": "Search speed vs. sensitivity." }
  }
}
```

The Setup form, progress display, CLI arg-string builder, and results plotting
all read this block — no other config files need changing.

### Step 2 — add the analysis script to `scripts/`

```
scripts/my_analysis.py
```

The script must have a `main(args, report=None)` function (for in-app use) **and**
an `if __name__ == "__main__":` block with `argparse` so it is also runnable as
a standalone CLI command. See any existing script (`blast_orthologs.py`, etc.) for
the pattern.

> **Dual-use rule (from CLAUDE.md):** never change the argparse interface of an
> existing script; only internal implementation may change.

### Step 3 — register a runner in `scripts/task_runners.py`

```python
def run_my_analysis(args, report=None):
    """Run my_analysis.py in-process."""
    from my_analysis import main as _main
    _main(args, report=report)

TASK_RUNNERS["my_analysis"] = run_my_analysis
```

That's it. The task now:
- appears in the Setup form's task-type dropdown
- shows parameters from its registry block (dropdowns for `choices` params)
- runs via the Progress tab and the CLI orchestrator
- has its output loaded and plotted by the Results tab (for `continuous`/`range` types)

### Verifying

1. Launch the app: `conda activate tagsites && python app_modular.py`
2. Add a task of type `my_analysis`, configure it, and click **Save Analysis**.
3. Confirm the written JSON has `{"type": "my_analysis", "args": {...}}` — no tooltips, no choices, no global keys copied in.
4. Run it from the Progress tab and confirm output loads in Results.
5. Run the same JSON through the CLI: `python scripts/run_tag_sites_from_json.py -i path/to/run.json` and confirm identical output.
