"""setup_ui.py — Shiny UI module for the Setup tab."""

from pathlib import Path

from shiny import ui, module

from config import SELECTABLE_TASKS, GLOBAL_TOOLTIPS, DEFAULT_SPECIES

_PARAMS_DIR = Path(__file__).parent.parent / "params"


def _preset_choices():
    """Return {name: name} dict for selectize, with a leading empty option."""
    names = [p.stem for p in sorted(_PARAMS_DIR.glob("*.json"))]
    return {"": "— select —", **{n: n for n in names}}


def label_with_tip(text, tip=""):
    """Return a plain label or label + hoverable ⓘ icon when tip is provided."""
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


# one-line description shown when an analysis type is selected
TASK_DESCRIPTIONS = {
    "blast":         "NCBI BLAST search for orthologs across species; computes per-residue conservation (Jensen-Shannon divergence).",
    "plddt":         "Extracts per-residue AlphaFold confidence (pLDDT) and solvent accessibility (SASA) from a PDB structure.",
    "modifications": "Identifies post-translational modification sites (phosphorylation, ubiquitination, etc.) by regex pattern matching.",
    "domains":       "Annotates protein domains and functional families via the EBI InterPro API.",
    "scores":        "Computes a sliding-window amino acid property score (e.g. hydrophobicity) along the sequence.",
}

_STYLE = """
    /* base font */
    .tab-pane { font-size: 0.875rem; }
    .form-label { font-size: 0.8rem; font-weight: 500; margin-bottom: 0.15rem; color: #343a40; }
    .form-control, .form-select { font-size: 0.8rem; padding: 0.2rem 0.45rem; }
    .form-control-sm, .form-select-sm { font-size: 0.8rem; }

    /* collapse Bootstrap's per-input bottom margin */
    .shiny-input-container { margin-bottom: 0.35rem !important; }
    .mb-3 { margin-bottom: 0.35rem !important; }

    /* accordion */
    .accordion-button { font-size: 0.875rem; font-weight: 600; padding: 0.45rem 0.75rem; }
    .accordion-body   { padding: 0.5rem 0.75rem; }

    /* card */
    .card-header { padding: 0.3rem 0.75rem; font-size: 0.8rem; font-weight: 600; }
    .card-body   { padding: 0.45rem 0.75rem; }
    .card        { margin-bottom: 0.4rem; }

    /* task accordion */
    .task-accordion .accordion-button { background: #f8f9fa; color: #212529; }
    .task-accordion .accordion-button:not(.collapsed) { background: #e9ecef; }
    .task-type-badge {
        font-size: 0.7rem; font-weight: 400; color: #6c757d;
        background: #e9ecef; border-radius: 3px; padding: 1px 5px;
        margin-left: 0.4rem; vertical-align: middle;
    }
    /* tighten the param grid inside task cards */
    .task-accordion .layout-column-wrap { gap: 0.4rem 0.75rem !important; }

    /* tooltip icon */
    .tip-icon { cursor: help; color: #adb5bd; font-size: 0.75em; margin-left: 3px; }
    .tooltip-inner { max-width: 320px; }

    /* task description line */
    .task-desc { font-size: 0.78rem; color: #495057; margin-top: 0.15rem; min-height: 1.1em; }

    /* save bar */
    .save-bar {
        position: sticky; bottom: 0; background: #fff;
        border-top: 2px solid #dee2e6; padding: 0.5rem 1rem;
        margin-top: 0.5rem; z-index: 100;
        display: flex; align-items: center; gap: 0.75rem;
        justify-content: space-between;
    }
    .save-bar #save_status { font-size: 0.8rem; color: #6c757d; margin: 0; }
    .save-bar-right { display: flex; align-items: center; gap: 0.4rem; }
    .save-bar-right .form-label { margin: 0; white-space: nowrap; }
    .save-bar-right .shiny-input-container { margin-bottom: 0 !important; }

    /* inline task builder rows */
    .task-type-row, .task-label-row {
        display: flex; align-items: center; gap: 0.5rem; flex-wrap: wrap;
    }
    .task-type-row { margin-bottom: 0.3rem; }
    .task-label-row { margin-bottom: 0.3rem; }
    .task-type-row .shiny-input-container > label,
    .task-label-row .shiny-input-container > label { display: none; }
    .task-type-row .shiny-input-container,
    .task-label-row .shiny-input-container { margin-bottom: 0 !important; }
    .task-type-row .form-label { white-space: nowrap; margin: 0; }
    .task-label-row .form-label { white-space: nowrap; margin: 0; }
    .task-type-row .task-desc { margin: 0; }

    /* misc */
    .section-hint { font-size: 0.8rem; color: #6c757d; margin-bottom: 0.3rem; margin-top: 0; }
    .layout-column-wrap > * { min-width: 0; }

    /* preset row: selectize fills remaining space, button stays natural width */
    .preset-pair { display: flex; gap: 0.4rem; align-items: flex-end; }
    .preset-pair .selectize-wrapper { flex: 1; min-width: 0; max-width: 200px; }
    .preset-pair .selectize-wrapper .shiny-input-container { margin-bottom: 0 !important; }

    /* file inputs: shrink the Browse button to match other controls */
    .shiny-input-file .btn { font-size: 0.8rem !important; padding: 0.2rem 0.5rem !important; }
    .shiny-input-file .form-control { font-size: 0.8rem; }

    /* inputs panel: keep everything in a constrained block */
    .inputs-inner { max-width: 600px; }
    .inputs-inner .shiny-input-container { margin-bottom: 0.35rem !important; }
"""

_t = GLOBAL_TOOLTIPS


def _organism_choices():
    """Build {value: label} for the organism preset dropdown."""
    choices = {"": "— select organism —"}
    for name, taxid in DEFAULT_SPECIES.items():
        choices[name] = f"{name} ({taxid})" if taxid else name
    return choices


@module.ui
def setup_ui():
    return ui.page_fluid(

        ui.tags.style(_STYLE),

        ui.accordion(

            # ── Panel 1: Inputs ───────────────────────────────────────────────
            ui.accordion_panel(
                "1 · Inputs",
                ui.p("Configure your run and upload your sequence(s).", class_="section-hint"),

                ui.div(
                    ui.input_text("email",
                        label_with_tip("Email", _t.get("email", "")),
                        value="",
                        placeholder="required for EBI submissions", width="100%"),
                    ui.input_text("run_name",
                        label_with_tip("Analysis name", _t.get("run_name", "")),
                        placeholder="e.g., your-favorite-gene-tag", width="100%"),
                    # Row 2: output directory (full width of the container)
                    ui.input_text("working_dir",
                        label_with_tip("Output directory", _t.get("working_dir", "")),
                        placeholder="auto-filled", width="100%"),
                    # Row 3: file uploads (stacked full-width)
                    ui.div(
                        ui.input_file("input_file",
                            label_with_tip("Protein sequence (.fa, .fasta, or .pdb)",
                                           _t.get("input_file", "")),
                            accept=[".fa", ".fasta", ".pdb"], width="100%"),
                        ui.input_file("input_genomic",
                            label_with_tip("Genomic region FASTA (required for reagent design)",
                                           _t.get("input_genomic", "")),
                            accept=[".fasta", ".fa"], width="100%"),
                        style="margin-top: 0.75rem;",
                    ),
                    ui.div(
                        ui.input_select(
                            "organism",
                            label_with_tip("Organism",
                                "Sets the taxid for the AlphaFold DB lookup and pre-fills "
                                "conservation BLAST. Use 'Other (search…)' to find any species."),
                            choices=_organism_choices(),
                            selected="",
                            width="100%",
                        ),
                        ui.output_ui("organism_search_ui"),
                        ui.output_ui("lineage_ui"),
                        style="margin-top: 0.5rem;",
                    ),
                    class_="inputs-inner",
                ),
                value="inputs",
            ),

            # ── Panel 2: Analyses ─────────────────────────────────────────────
            ui.accordion_panel(
                "2 · Analyses",
                ui.p("Build the set of analyses to run.", class_="section-hint"),

                # preset load row
                ui.card(
                    ui.card_body(
                        ui.div(
                            ui.tags.label("Load saved analyses set", class_="form-label"),
                            ui.input_select("load_preset", label="",
                                choices=_preset_choices(), selected=""),
                            ui.input_action_button("load_preset_btn", "Load",
                                class_="btn-sm btn-outline-secondary"),
                            class_="task-type-row",
                        ),
                    ),
                ),

                # task builder
                ui.card(
                    ui.card_body(
                        # Row 1: Analysis type | dropdown | description
                        ui.div(
                            ui.tags.label("Analysis type", class_="form-label"),
                            ui.input_select("task_type", label="", choices=SELECTABLE_TASKS),
                            ui.span(ui.output_text("task_type_desc"), class_="task-desc"),
                            class_="task-type-row",
                        ),
                        # Row 2: Task label | text input
                        ui.div(
                            ui.tags.label("Task label", class_="form-label"),
                            ui.input_text("task_label", label="",
                                placeholder="short ID, e.g. WORM or BROAD"),
                            ui.tags.small(
                                "Names output files, e.g. WORM_blast.jsd.",
                                class_="text-muted",
                            ),
                            class_="task-label-row",
                        ),
                        ui.div(
                            ui.input_action_button("add_task", "＋ Add task",
                                disabled=True, class_="btn-primary btn-sm"),
                            style="margin-top: 0.5rem;",
                        ),
                    ),
                ),

                # dynamically rendered task accordion
                ui.output_ui("task_cards"),

                value="analyses",
            ),

            open=["inputs", "analyses"],
            multiple=True,
        ),

        # ── Sticky save bar ───────────────────────────────────────────────────
        ui.div(
            # left: save + status
            ui.div(
                ui.input_action_button("save_analysis", "💾 Save Analysis",
                    disabled=True, class_="btn-success btn-sm"),
                ui.output_text("save_status"),
                style="display:flex; align-items:center; gap:0.75rem;",
            ),
            # right: save preset
            ui.div(
                ui.tags.label("Save current analyses as default", class_="form-label"),
                ui.input_text("preset_name", label="", placeholder="name"),
                ui.input_action_button("save_preset_btn", "Save",
                    disabled=True, class_="btn-sm btn-outline-secondary"),
                class_="save-bar-right",
            ),
            class_="save-bar",
        ),
    )
