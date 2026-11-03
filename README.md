# TAGSITES

A Shiny app for finding functional sites in a protein sequence and designing CRISPR knock-in
reagents to tag the protein at those sites.

Given a protein sequence (FASTA, UniProt hit, or AlphaFold PDB), TAGSITES runs a configurable
analysis pipeline:

- **Conservation**: orthologs from BLAST (web) or DIAMOND (local), aligned and scored
- **Structure**: AlphaFold pLDDT and solvent-accessible surface area
- **Modifications**: PTM sites from regex patterns
- **Domains**: InterPro (web) or Pfam via `pyhmmer` (local)
- **UniProt features**: curated lipidation, binding sites, modified residues, and more
- **Scores**: sliding-window amino acid properties
- **Reagents**: guides, homology arms, oligos, genotyping primers and annotated GenBank files
  for every candidate insertion site. The gene model must translate to the input protein, or
  the reagents are not designed (see [CLAUDE.md](CLAUDE.md), "Gene model must translate").

Results are interactive Plotly plots. Select sites, then take the reagent tables to the bench.

TAGSITES runs in four modes. Each has one guide:

| Goal | Guide |
|---|---|
| Analyze one protein using web services only (no downloads) | this page |
| Analyze one protein with every analysis on local hardware | [docs/LOCAL_SETUP.md](docs/LOCAL_SETUP.md) |
| Open precomputed results for a whole proteome in the main app | [docs/LOCAL_DATABASE.md](docs/LOCAL_DATABASE.md) |
| Review and design tags for many proteins at once (batch app) | [docs/BATCH_APP.md](docs/BATCH_APP.md) |
| Build the precomputed database for an organism | [docs/BUILD_LOCAL_DATABASE.md](docs/BUILD_LOCAL_DATABASE.md) |
| Drive the pipeline without the browser | [docs/CLI.md](docs/CLI.md) |

Shared reference: [docs/INSTALL.md](docs/INSTALL.md) (full install) and
[docs/CONFIG.md](docs/CONFIG.md) (the config files and `TAGSITES_BATCH_CONFIG`).

---

## Quick start (web services)

Requires an internet connection. BLAST, Clustal Omega, InterPro and Genewise run through the EBI
REST API, so no local bioinformatics tools are needed.

### 1. Get the code

Clone the repository:

```bash
git clone https://github.com/msr2009/TAGSITES.git
cd TAGSITES
```

Without git, click **⬇ Download App** in the top-right corner of the running app, then
`unzip tagsites-app.zip && cd tagsites-app`.

### 2. Create the environment

```bash
mamba env create -f environment.yml     # conda works the same way
mamba activate tagsites
```

Python 3.10 is required, not incidental: `environment.yml` pins numpy 1.x and old
scikit-learn/lightgbm solely for the `rs3` guide-scoring package. A web-only install can use
`python3.10 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt`
instead; local backends and batch mode need the conda environment. Details and pitfalls are in
[docs/INSTALL.md](docs/INSTALL.md).

### 3. Start the app

```bash
python app_modular.py
```

Open **http://127.0.0.1:8000**. The tabs run left to right:

1. **Setup Analyses**: load a sequence, choose analyses, enter an e-mail address (EBI requires
   one to accept jobs), and save.
2. **Run Analyses**: submit and watch each task.
3. **View Results**: plots, tracks and alignments.
4. **Create Reagents**: pick insertion sites and download the reagent ZIP.

### 4. Check the install (optional)

```bash
pytest                                                   # offline tests
TAGSITES_EMAIL=you@example.edu pytest -m network --run-network   # adds EBI/network tests
```

---

## Going further

- **Faster or offline analysis** on one machine: [docs/LOCAL_SETUP.md](docs/LOCAL_SETUP.md).
- **Precomputed proteome**: [docs/LOCAL_DATABASE.md](docs/LOCAL_DATABASE.md) covers getting a
  database and the main app's *Load precomputed results* card;
  [docs/BATCH_APP.md](docs/BATCH_APP.md) covers the batch review app (`python app_batch.py`,
  port 8001).
- **Another organism**: [docs/BUILD_LOCAL_DATABASE.md](docs/BUILD_LOCAL_DATABASE.md).

The local workflows read a config file chosen with `TAGSITES_BATCH_CONFIG`; start from
`batch.config.example.json` ([docs/CONFIG.md](docs/CONFIG.md)). With no variable set, every
analysis uses web services, which is also how the hosted app behaves.

## Repository map

```
app_modular.py / server.py / ui.py     main app (port 8000)
app_batch.py / server_batch.py / ui_batch.py   batch review app (port 8001, local only)
task_definitions.json                  task registry and default analyses
batch.config.json                      checked-in default config (do not edit for local work)
batch.config.example.json              template for a machine-local config
modules/                               Shiny UI and server modules (setup, progress, results, reagents, batch_*)
scripts/                               analysis scripts; each is also a standalone CLI
utils/                                 result loading and plotting
tables/  params/  tests/               PTM patterns and property tables, saved presets, pytest suite
docs/                                  guides listed above
```

[CLAUDE.md](CLAUDE.md) holds the full per-script map and the design rules.

Developer notes: [docs/DATA_MODEL.md](docs/DATA_MODEL.md) (file types, adding a task),
[docs/GNOMAD_CONSTRAINT_SPEC.md](docs/GNOMAD_CONSTRAINT_SPEC.md),
[docs/ANNOTATE_INTERNAL_TAGS.md](docs/ANNOTATE_INTERNAL_TAGS.md),
[docs/SUGGESTED_SITES_REVISION.md](docs/SUGGESTED_SITES_REVISION.md).

## Vendored clients

`genewise.py`, `iprscan5.py`, `ncbiblast.py`, `dbfetch.py` and `clustalo.py` in `scripts/` are
EBI web service clients from
[ebi-jdispatcher/webservice-clients](https://github.com/ebi-jdispatcher/webservice-clients),
used as-is. Replace them from that repository after an EBI API change.
