# TAGSITES

A Shiny-based web application for identifying functional sites in protein sequences to guide protein tagging and CRISPR reagent design.

## What it does

Given a protein sequence (FASTA or PDB), the app runs a configurable analysis pipeline:
- **BLAST**: finds orthologs across selected organisms to compute conservation
- **pLDDT**: extracts AlphaFold2 confidence scores and solvent-accessible surface area (SASA) from PDB structures
- **Modifications**: identifies PTM sites (phosphorylation, ubiquitination, etc.) using regex patterns
- **Domains**: calls the EBI InterPro API to annotate protein domains
- **Scores**: computes sliding-window amino acid property scores (hydrophobicity, etc.)

Results are visualized as interactive Plotly plots. Users can then design CRISPR guides targeting selected sites.

## Running the app

```bash
conda activate tagsites
python app_modular.py
```

## Project layout

```
app_modular.py          # main Shiny entry point
server.py / ui.py       # top-level Shiny orchestrators
config.py               # species taxonomy, result type config, JSON defaults
default_json.json       # default task parameters and script-to-task mappings

modules/                # Shiny UI + server components
  setup_ui.py / setup_server.py       # analysis configuration, file upload, task params
  progress_ui.py / progress_server.py # job progress display
  results_ui.py / results_server.py   # plot rendering, alignment visualization
  reagents_ui.py / reagents_server.py # CRISPR reagent design

scripts/                # core analysis executables
  run_tag_sites_from_json.py  # async pipeline orchestrator (main workhorse)
  run_tag_sites.py            # CLI pipeline runner
  blast_orthologs.py          # NCBI BLAST homolog search
  extract_from_pdb.py         # pLDDT + SASA from AlphaFold PDB
  regex_sites.py              # PTM site identification
  call_interpro.py            # EBI InterPro domain annotation
  calculate_protein_scores.py # sliding-window property scoring
  site_selection_util.py      # shared library: FASTA/PDB I/O, BLAST API, sequence utils
  existing_AF_model.py        # search AFDB for existing predictions
  design_guides_across_region.py  # CRISPR guide design

utils/
  results.py    # load JSON output → DataFrames; Plotly + matplotlib visualization
  helpers.py    # taxonomy loading, Shiny reactive state helpers

tables/
  modification_sites.txt           # regex patterns for PTM sites
  hydrophobicity_kyle-doolittle.tsv # amino acid property scores

params/
  worm_default.json   # saved parameter preset for C. elegans
  *.json              # user-saved parameter presets
```

## Dual-use constraint: Shiny app + standalone CLI

**All scripts in `scripts/` must remain runnable from the command line independently of the Shiny app.** When refactoring script internals, never change argparse interfaces, CLI flag names, or `__main__` entry points. Only internal implementation may change (e.g. swapping a subprocess call for an in-process function call). The Shiny app calls script `main()` functions directly via `scripts/task_runners.py`; the CLI calls the same scripts as subprocesses via `scripts/run_tag_sites_from_json.py`. Both paths must produce identical outputs.

## Key design patterns

**Configuration-driven pipeline**: `default_json.json` defines which scripts map to which analyses and their default parameters. `config.py` defines available species (with NCBI taxonomy IDs) and which result types are continuous vs. range-based.
**Shiny reactive state**: analysis parameters are stored in a shared reactive dict (`shared_dict`) passed between modules. `utils/helpers.py:update_shared_dict()` handles updates.
**Async job submission**: `run_tag_sites_from_json.py` spawns analysis scripts as subprocesses and polls for completion, enabling parallel execution of independent analyses.
**Alignment rendering**: sequence alignments are pre-rendered as matplotlib SVGs (stored in a dict keyed by alignment name) and displayed statically — this replaced an earlier real-time Plotly approach that was too slow.

## Environment

Key dependencies: `shiny`, `biopython`, `plotly`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `requests`, `clustalo` (CLI).

## External services

- **NCBI BLAST** (via EBI REST API) — `scripts/site_selection_util.py:ncbiblast_call()`
- **EBI InterPro API** — `scripts/call_interpro.py`
- **AlphaFold DB** — `scripts/existing_AF_model.py`
- **UniProt taxonomy** — `uniprot_species.flat.txt` (local flat file, ~1.7MB)
