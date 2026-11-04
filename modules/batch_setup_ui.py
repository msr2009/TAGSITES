"""batch_setup_ui.py — Setup tab of the batch app: open a batch, or build one from a list.

Styled like the single-protein Setup tab (modules/setup_ui.py): numbered accordion panels with
a grey hint in each header, compact inputs, and a sticky save bar. The stylesheet and helpers
are imported from there so the two apps stay in step.
"""

from shiny import module, ui

from modules.setup_ui import _STYLE, label_with_tip
from modules.ui_helpers import COMPACT_FILE_CSS, compact_file_input

_LIST_PLACEHOLDER = ("One protein per line: gene name, UniProt accession or isoform id,\n"
                     "WBGene id or WormBase sequence / transcript name. An insert name may\n"
                     "follow, e.g.\n"
                     "trxr-1\tFLAG\nunc-54\nQ17745-2\n\n"
                     "FASTA records are matched by exact sequence.")

_HINT = "font-size:0.8rem; font-weight:400; color:#6c757d; margin-left:0.75rem;"

# batch-only additions to the single app's Setup stylesheet
_BATCH_STYLE = """
    .ts-inline-row { display: flex; align-items: center; gap: 0.5rem; flex-wrap: wrap;
                     margin-bottom: 0.3rem; }
    .ts-inline-row .shiny-input-container { margin-bottom: 0 !important; width: auto !important; }
    .ts-inline-row .form-check { margin-bottom: 0; }
    .ts-inline-row .form-label { white-space: nowrap; margin: 0; }
    .ts-inline-row input[type=number] { width: 5rem; }
    .ts-filter-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(230px, 1fr));
                      gap: 0.1rem 1rem; max-width: 760px; }
    .ts-filter-grid .shiny-input-container { margin-bottom: 0.2rem !important; }
    .ts-filter-grid .form-check-label, .ts-mask-group .form-check-label { font-size: 0.8rem; }
    .ts-masks-title { font-weight: 700; font-size: 0.875rem; color: #212529;
                      margin: 0.6rem 0 0.2rem; }
    .ts-mask-group { margin: 0 0 0.35rem 0.25rem; }
    .ts-mask-group-label { font-size: 0.78rem; font-weight: 600; color: #495057;
                           margin-bottom: 0.1rem; }
    .ts-mask-group .ts-inline-row { flex-wrap: nowrap; }
    .ts-repair .form-check-inline { margin-right: 1rem; }
    .ts-repair .form-check-label { font-size: 0.8rem; white-space: nowrap; }
    .ts-db-info { font-size: 0.8rem; color: #343a40; margin: 0 0 0.4rem; }
    .ts-db-info td { padding: 0.05rem 0.75rem 0.05rem 0; vertical-align: top; }
    .ts-db-info td:first-child { color: #6c757d; white-space: nowrap; }
    .ts-db-info code { font-size: 0.78rem; word-break: break-all; }
    .ts-resolved { font-size: 0.8rem; margin: 0.4rem 0 0; }
    .ts-resolved th { font-weight: 600; color: #343a40; background: #f1f5fb !important;
                      position: sticky; top: 0; z-index: 1; }
    .ts-resolved td { vertical-align: middle; padding: 0.15rem 0.4rem; }
    .ts-resolved .shiny-input-container { margin-bottom: 0 !important; }
    .ts-resolved-wrap { max-height: 420px; overflow-y: auto; border: 1px solid #dee2e6;
                        border-radius: 4px; }
    .ts-resolved-summary { font-size: 0.8rem; font-weight: 600; color: #343a40;
                           margin-top: 0.6rem; }
    .ts-resolved-summary .ts-miss { color: #9b1c1c; }
    .save-bar .btn { white-space: nowrap; }
"""


def _hint(title, text):
    """Accordion panel title with the grey one-line hint the single app uses."""
    return ui.span(title, ui.span(text, style=_HINT))


def _inline(*children):
    """A label + input row (the inputs' own labels are omitted)."""
    return ui.div(*children, class_="ts-inline-row")


@module.ui
def batch_setup_ui():
    """Open a saved batch, or build one: proteins, reagents, site filters, EBI; then save."""
    return ui.page_fluid(
        ui.tags.style(_STYLE + COMPACT_FILE_CSS + _BATCH_STYLE),
        ui.accordion(
            # ── Panel 0: the proteome database ────────────────────────────────
            ui.accordion_panel(
                _hint("0 · Proteome database",
                      "The precomputed analyses every protein is looked up in."),
                ui.output_ui("db_card"),
                ui.div(
                    _inline(
                        ui.div(ui.input_text("db_path", None, value="", width="100%",
                                             placeholder="/path/to/proteome.sqlite3"),
                               style="flex:1; min-width:0;"),
                        ui.input_action_button("use_db", "Use database",
                                               class_="btn-sm btn-outline-primary")),
                    ui.p("Applies to this app until it restarts; a new batch records its "
                         "database and reopening the batch switches back to it.",
                         class_="section-hint"),
                    class_="inputs-inner",
                ),
                value="db",
            ),

            # ── Panel 1: open a saved batch ───────────────────────────────────
            ui.accordion_panel(
                _hint("1 · Open a saved batch",
                      "Resume a batch on this machine, or restore a batch ZIP."),
                ui.div(ui.output_ui("open_batch_card"), class_="inputs-inner"),
                value="open",
            ),

            # ── Panel 2: proteins ─────────────────────────────────────────────
            ui.accordion_panel(
                _hint("2 · Proteins",
                      "Paste or upload a list; each name is looked up in the proteome "
                      "database."),
                ui.div(
                    ui.input_text("batch_name",
                                  label_with_tip("Batch name",
                                                 "Names the batch folder and the output "
                                                 "files. Defaults to batch_<date>."),
                                  value="", placeholder="e.g. alz_genes", width="100%"),
                    ui.input_text_area("protein_list", "Proteins", rows=8, width="100%",
                                       placeholder=_LIST_PLACEHOLDER),
                    ui.div(ui.span("— or upload —"), class_="or-divider"),
                    compact_file_input("list_file", "List file (txt / tsv / csv / FASTA)",
                                       accept=[".txt", ".tsv", ".csv", ".fa", ".fasta"]),
                    ui.input_action_button("resolve", "Look up proteins",
                                           class_="btn-primary btn-sm mt-2"),
                    class_="inputs-inner",
                ),
                ui.output_ui("resolved_card"),
                value="proteins",
            ),

            # ── Panel 3: reagents ─────────────────────────────────────────────
            ui.accordion_panel(
                _hint("3 · Reagents",
                      "Repair template and insert(s) for every accepted site."),
                ui.div(
                    ui.div(
                        ui.input_radio_buttons(
                            "repair_type", "Repair template",
                            {"ha": "Homology arms", "ssodn": "ssODN", "amplicon": "PCR amplicon",
                             "plasmid": "SapTrap plasmid"},
                            selected="ssodn", inline=True, width="100%"),
                        class_="ts-repair"),
                    _inline(ui.tags.label("Arm length (bp)", class_="form-label"),
                            ui.input_numeric("arm_length", None, value=45, min=20, max=1000)),
                    ui.output_ui("insert_select"),
                    ui.panel_conditional(
                        "input.inserts && input.inserts.indexOf('Custom') >= 0",
                        ui.input_text("custom_insert", "Custom insert DNA", value="",
                                      width="100%")
                    ),
                    ui.p("A protein's own insert (second column of the list) replaces these "
                         "for that protein. Guides default to the one nearest each site and "
                         "can be changed per site in Review.", class_="section-hint mt-1"),
                    class_="inputs-inner",
                ),
                value="reagents",
            ),

            # ── Panel 4: suggested-site filters ───────────────────────────────
            ui.accordion_panel(
                _hint("4 · Suggested sites",
                      "How many sites to suggest per protein, and where not to."),
                _inline(ui.tags.label("Suggested sites per protein", class_="form-label"),
                        ui.input_numeric("n_suggested", None, value=2, min=1, max=20)),
                ui.div("Masks", class_="ts-masks-title"),
                ui.div(
                    ui.div("Position", class_="ts-mask-group-label"),
                    ui.div(
                        ui.input_checkbox("filt_no_nterm", "No N-terminus (residue 1)",
                                          value=False),
                        ui.input_checkbox("filt_no_cterm", "No C-terminus (last residue)",
                                          value=False),
                        class_="ts-filter-grid"),
                    _inline(ui.input_checkbox("filt_guide", "Guide cut within", value=False),
                            ui.input_numeric("filt_guide_bp", None, value=5, min=0, max=100),
                            ui.span("bp", class_="form-label")),
                    class_="ts-mask-group",
                ),
                ui.div(
                    ui.div("Annotation", class_="ts-mask-group-label"),
                    ui.div(
                        ui.input_checkbox("filt_constitutive", "Constitutive sites only",
                                          value=False),
                        ui.input_checkbox("filt_mask_signal", "Signal / transit peptide",
                                          value=True),
                        ui.input_checkbox("filt_mask_tm", "Transmembrane", value=True),
                        ui.input_checkbox("filt_mask_inside", "Inside (cytosolic)",
                                          value=False),
                        ui.input_checkbox("filt_mask_outside",
                                          "Outside (lumenal / extracellular)", value=False),
                        ui.input_checkbox("filt_mask_domain", "Domains (Pfam)", value=False),
                        ui.input_checkbox("filt_mask_uniprot", "All UniProt sites",
                                          value=False),
                        class_="ts-filter-grid"),
                    class_="ts-mask-group",
                ),
                ui.p("Masks remove positions from every protein's suggestions; each can be "
                     "switched off per protein in Review (Rescoring), and any site can still "
                     "be added by hand.", class_="section-hint mt-1"),
                value="filters",
            ),

            # ── Panel 5: proteins not in the database ─────────────────────────
            ui.accordion_panel(
                _hint("5 · Proteins not in the database",
                      "Optional: analyse them on EBI (slow)."),
                ui.div(
                    ui.input_checkbox("ebi_for_missing",
                                      "Run EBI analyses for proteins not in the database",
                                      value=False, width="100%"),
                    ui.panel_conditional(
                        "input.ebi_for_missing",
                        ui.input_text("email",
                                      label_with_tip("Email", "Required by EBI for job "
                                                     "submission."),
                                      value="", placeholder="required for EBI submissions",
                                      width="100%"),
                        ui.p("Only UniProt accessions and FASTA records can be analysed; "
                             "EBI runs take minutes per protein, at most 3 at once.",
                             class_="section-hint")
                    ),
                    class_="inputs-inner",
                ),
                value="ebi",
            ),

            open=["db", "open", "proteins", "reagents", "filters"],
            multiple=True,
        ),

        # ── Sticky save bar ───────────────────────────────────────────────────
        ui.div(
            ui.input_action_button("save_continue", "6 · 💾 Save batch & continue →",
                                   class_="btn-success btn-sm"),
            class_="save-bar",
        ),
    )
