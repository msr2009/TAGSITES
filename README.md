# TAGSITES

A Shiny-based web application for identifying functional sites in protein sequences to guide
protein tagging and CRISPR reagent design.

Given a protein sequence (FASTA or PDB), the app runs a configurable analysis pipeline:
- **BLAST** — finds orthologs across selected organisms to compute conservation
- **pLDDT** — extracts AlphaFold2 confidence scores and solvent-accessible surface area (SASA)
- **Modifications** — identifies PTM sites (phosphorylation, ubiquitination, etc.)
- **Domains** — calls the EBI InterPro API to annotate protein domains
- **Scores** — computes sliding-window amino acid property scores (hydrophobicity, etc.)
- **Reagents** — designs CRISPR knock-in reagents (guide RNAs + HDR homology arms) for every
  potential tag-insertion site, using Genewise to map the protein onto a supplied genomic region

Results are visualized as interactive Plotly plots. Users select a site and the reagents page
provides guide spacers and mutated homology arms ready to order.

---

## Install

### 1. Clone the repository

```bash
git clone https://github.com/msr2009/TAGSITES.git
cd TAGSITES
```

### 2. Create the conda environment

A portable environment spec is provided:

```bash
conda env create -f environment.yml
conda activate tagsites
```

> **Exact-pin reference:** `tagsites.yml` is a full osx-64 build-pinned export of the original
> development environment. Use it on macOS to reproduce the exact environment:
> `conda env create -f tagsites.yml`

### 3. Register your EBI e-mail

All network analyses (BLAST, InterPro, Genewise) are submitted via the EBI REST API and require
a registered e-mail address. Set it in the `global.email` field of your run JSON, or in the
app's Setup page.

---

## Run the app

```bash
conda activate tagsites
python app_modular.py
```

Open the printed URL in a browser. Use the Setup page to upload a sequence, configure analyses,
and submit a run.

---

## Run from the command line

The pipeline can be driven entirely via a JSON config file — useful for batch processing or
testing.

```bash
conda activate tagsites
python scripts/run_tag_sites_from_json.py -i <run.json>
```

A minimal example JSON (using pre-computed Genewise output so no network is needed) is at
[`tests/data/example_run.json`](tests/data/example_run.json). Outputs land in the directory
specified by `global.working_dir`.

---

## Run the tests

The test suite is split into fast offline unit tests and network-dependent integration tests.

```bash
# Fast offline tests only (no network, no EBI account needed)
conda activate tagsites
pytest

# Include network integration tests (requires internet + EBI e-mail)
TAGSITES_EMAIL=your@email.com pytest -m network --run-network
```

---

## EBI webservice clients

The scripts `genewise.py`, `iprscan5.py`, `ncbiblast.py`, `dbfetch.py`, and `clustalo.py` in
`scripts/` are EBI webservice clients vendored from
[ebi-jdispatcher/webservice-clients](https://github.com/ebi-jdispatcher/webservice-clients).
They are used as-is; no separate download is required. If they need refreshing (e.g. after an
EBI API change), replace them with the latest versions from that repository.

---

## Project layout

```
app_modular.py              # main Shiny entry point
config.py                   # species taxonomy, result type config, JSON defaults
default_json.json           # default task parameters and script-to-task mappings

modules/                    # Shiny UI + server components
scripts/                    # core analysis executables + EBI REST clients
  run_tag_sites_from_json.py    # async pipeline orchestrator
  blast_orthologs.py            # NCBI BLAST homolog search
  extract_from_pdb.py           # pLDDT + SASA from AlphaFold PDB
  regex_sites.py                # PTM site identification
  call_interpro.py              # EBI InterPro domain annotation
  calculate_protein_scores.py   # sliding-window property scoring
  crispr_util.py                # guide finding, codon table, PAM mutation logic
  parse_genewise.py             # Genewise GFF parser → insertion-site table
  run_genewise.py               # both-strand Genewise runner + orientation selection
  design_tag_reagents.py        # reagent table: residue × guide + HDR arms
  find_guides.py                # CLI: enumerate guide cut sites
  site_selection_util.py        # shared FASTA/PDB I/O, BLAST API, sequence utils

utils/
  results.py                # load JSON output → DataFrames; Plotly + matplotlib viz

tables/
  modification_sites.txt              # regex patterns for PTM sites
  hydrophobicity_kyle-doolittle.tsv   # amino acid property scores

params/                     # saved parameter presets (JSON)
tests/                      # pytest test suite
  data/                     # committed test fixtures
```
