"""batch_setup_ui.py — Setup tab of the batch app: open a batch, or build one from a list."""

from shiny import module, ui

_LIST_PLACEHOLDER = ("One protein per line: gene name, UniProt accession or isoform id,\n"
                     "WBGene id or WormBase transcript. An insert name may follow, e.g.\n"
                     "trxr-1\tFLAG\nunc-54\nQ17745-2\n\n"
                     "FASTA records are matched by exact sequence.")


@module.ui
def batch_setup_ui():
    """Open-batch card, protein list + resolution table, reagent settings, Save & continue."""
    return ui.page_fluid(
        ui.layout_columns(
            ui.card(
                ui.card_header("Open a saved batch"),
                ui.output_ui("open_batch_card"),
            ),
            ui.card(
                ui.card_header("New batch"),
                ui.input_text("batch_name", "Batch name", value=""),
                ui.input_text_area("protein_list", "Proteins", rows=8, width="100%",
                                   placeholder=_LIST_PLACEHOLDER),
                ui.input_file("list_file", "…or upload a list (txt / tsv / csv / FASTA)",
                              accept=[".txt", ".tsv", ".csv", ".fa", ".fasta"]),
                ui.input_action_button("resolve", "Look up proteins", class_="btn-primary"),
            ),
            col_widths=(4, 8),
        ),
        ui.output_ui("resolved_card"),
        ui.card(
            ui.card_header("Reagents for every accepted site"),
            ui.layout_columns(
                ui.input_radio_buttons(
                    "repair_type", "Repair template",
                    {"ha": "Homology arms", "ssodn": "ssODN", "amplicon": "PCR amplicon",
                     "plasmid": "SapTrap plasmid"},
                    selected="ssodn"),
                ui.div(
                    ui.input_numeric("arm_length", "Arm length (bp)", value=45, min=20,
                                     max=1000),
                    ui.output_ui("insert_select"),
                    ui.input_text("custom_insert", "Custom insert DNA (when 'Custom' is chosen)",
                                  value="", width="100%"),
                ),
                ui.div(
                    ui.input_checkbox("ebi_for_missing",
                                      "Run EBI analyses for proteins not in the database",
                                      value=False),
                    ui.panel_conditional(
                        "input.ebi_for_missing",
                        ui.input_text("email", "Email (required by EBI)", value=""),
                        ui.p("Only UniProt accessions and FASTA records can be analysed; "
                             "EBI runs take minutes per protein, at most 3 at once.",
                             class_="text-muted small"),
                    ),
                ),
                col_widths=(3, 5, 4),
            ),
            ui.div(
                ui.tags.label("Suggested-site filters", class_="form-label fw-semibold"),
                ui.layout_columns(
                    ui.div(
                        ui.tags.span("Suggested sites per protein"),
                        ui.input_numeric("n_suggested", None, value=2, min=1, max=20,
                                         width="80px"),
                        style="display:flex; gap:6px; align-items:baseline;"),
                    ui.input_checkbox("filt_no_nterm", "No N-terminus (residue 1)", value=False),
                    ui.input_checkbox("filt_no_cterm", "No C-terminus (last residue)",
                                      value=False),
                    ui.div(
                        ui.input_checkbox("filt_guide", "Guide cut within", value=False),
                        ui.input_numeric("filt_guide_bp", None, value=5, min=0, max=100,
                                         width="90px"),
                        style="display:flex; gap:6px; align-items:baseline;"),
                    ui.input_checkbox("filt_constitutive", "Constitutive sites only",
                                      value=False),
                    col_widths=(4, 2, 2, 2, 2),
                ),
                ui.p("These mask positions from the suggested sites of every protein; each "
                     "can be switched off per protein in Review (Rescoring, Advanced), and any "
                     "site can still be added by hand.", class_="text-muted small"),
            ),
            ui.p("A protein's own insert (second column of the list) replaces the batch "
                 "inserts for that protein. Guides default to the one nearest each site and "
                 "can be changed per site while reviewing.", class_="text-muted small"),
        ),
        ui.input_action_button("save_continue", "Save batch & continue →",
                               class_="btn-success btn-lg"),
    )
