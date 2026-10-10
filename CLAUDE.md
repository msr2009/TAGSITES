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
task_definitions.json   # task registry: script-to-task mappings, params, and the
                         # auto-applied default_tasks set (see below)

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
  uniprot_features.py         # curated UniProt feature annotation (lipidation, PTMs, binding sites, modified residues,
                              # mutagenesis positions, ...)
  calculate_protein_scores.py # sliding-window property scoring
  site_selection_util.py      # shared library: FASTA/PDB I/O, BLAST API, sequence utils
  existing_AF_model.py        # search AFDB for existing predictions
  uniprot_api.py               # shared UniProt REST helpers (checksum lookup, entry fetch)
  http_retry.py               # shared bounded-retry + wall-clock-deadline wrapper around
                              # requests; used by ebi_rest.py and ensembl_rest.py
  ensembl_rest.py             # Ensembl REST client for the genomic-sequence auto-fetch
                              # (knobs in ensembl.config.json)
  cds_check.py                # gene-model gate: translates the model's CDS and compares it to
                              # the input protein (length differs = hard fail, substitutions =
                              # warning + {run}.reagents.model_check.json for the UI banner)
  genbank_input.py            # user GenBank (CDS or exon features) used as the gene model in place
                              # of Genewise; writes the same *_genewise.* files, then the same
                              # cds_check gate applies. Sequence-only GenBank falls through to Genewise
  guide_efficiency.py         # RS3 on-target guide scoring (optional; fails soft)
  genbank_export.py           # annotated GenBank (.gb) records for ApE/SnapGene —
                              # two per selected guide, each spanning the WHOLE
                              # genomic region: the WT locus and the same region
                              # carrying the knock-in. Always included in the
                              # reagents download ZIP
  design_guides_across_region.py  # standalone CLI guide finder — NOT used by the app
                                  # (the app uses crispr_util.find_guides via
                                  #  design_tag_reagents.py); its off-target code
                                  #  needs a local blastdb + bedtools and is legacy
                                  #  and unwired — the live screen is the three
                                  #  offtarget_*.py modules below
  offtarget_screen.py         # network-free core of the off-target / primer screen:
                              # classification, post-filter, window scoring,
                              # amplicon prediction. No backend touches it
  offtarget_blat.py           # Screen A via the UCSC REST BLAT endpoint (seconds).
                              # Needs an API key — see offtarget.config.json's blat
                              # block and ucsc.local.json
  offtarget_remote.py         # Screen A via EBI blastn (minutes) plus Screen B, the
                              # concatenated spacer query. Screen B cannot move to
                              # BLAT: a 20 nt query with mismatches is below its
                              # tiling floor
  offtarget_local.py          # Screens A and B plus a genome-wide genotyping-primer
                              # screen on a local BLAST+ database (seconds). Used when
                              # batch.config.json's backends.offtarget_region is "local"
                              # AND blastn + the database exist, else the choice falls
                              # back to BLAT/EBI, so the deployed app (no BLAST+) is
                              # unchanged. Build the database with
                              # `python scripts/reference_data.py --only blastdb`

  proteome_run.py             # batch driver: runs the analysis tasks for every protein of
                              # the proteome (one folder per protein, _status.jsonl)
  build_proteome_db.py        # consolidates a proteome_run into ONE SQLite file (proteins,
                              # features, per-residue tracks, tag-site scores, reagents).
                              # Run --estimate first: it samples, measures and projects
                              # the size and run time; a full build needs --yes
  proteome_db.py              # read-only lookup over that file: find_protein(), get_*(),
                              # query(), regenerate_arms(); CLI `lookup <name>` / `sql`

docs/
  LOCAL_SETUP.md  # step-by-step guide to running the app with all analyses local
                  # (single-protein use; proteome-scale setup is not covered)

utils/
  results.py    # load JSON output → DataFrames; Plotly + matplotlib visualization
  helpers.py    # taxonomy loading, Shiny reactive state helpers

tables/
  modification_sites.txt           # regex patterns for PTM sites
  hydrophobicity_kyte-doolittle.tsv # amino acid property scores

params/
  worm_default.json   # example saved analyses preset (C. elegans) — not auto-loaded
  *.json              # user-saved parameter presets
```

**Off-target backend order**: `design_tag_reagents._region_backend` tries local BLAST+,
then UCSC BLAT, then EBI. Screen B and the primer screen run locally only when Screen A did,
because the intended-locus exclusion is keyed on Screen A's self spans (chromosome names,
not ENA accessions). Local-screen knobs live in `offtarget.config.json`'s `local` block.
Recall of the local spacer screen was measured at 100% (123 planted sites, 1-3 mismatches
in positions 1-15) only with `local.spacer_evalue` >= 1e5: blastn drops a 3-mismatch 20 nt
site at the shared `blast.evalue_spacer` of 1000.

**Proteome database**: `scripts/build_proteome_db.py` builds `proteome.sqlite3` in stages
(`ingest`, `scores`, `reagents`); `scripts/proteome_db.py` reads it (`lookup trxr-1`).
Always `--estimate` first (disk is limited): it measures a sample with SQLite's dbstat and
projects size and time. The file is assembled on LOCAL disk and copied out, since SQLite
locking is unreliable on the network volume. A task's file is ingested only if its last
status is `ok`. `proteins.sequence` is the sequence the batch analysed (the folder's `.fa`),
flagged by `sequence_differs` where it differs from `local_store` (28 proteins). Reagent tables never store the ~1 kb arms; `regenerate_arms()` rebuilds them
from the genome and was verified byte-identical to the reagent TSV. The reagent stage needs
`TAGSITES_BATCH_CONFIG=batch.config.local.json` (local genewise backend) and is resumable.
Knobs live in `batch.config.json`'s `proteome_db` block.

## Dual-use constraint: Shiny app + standalone CLI

**All scripts in `scripts/` must remain runnable from the command line independently of the Shiny app.** When refactoring script internals, never change argparse interfaces, CLI flag names, or `__main__` entry points. Only internal implementation may change (e.g. swapping a subprocess call for an in-process function call). The Shiny app calls script `main()` functions directly via `scripts/task_runners.py`; the CLI calls the same scripts as subprocesses via `scripts/run_tag_sites_from_json.py`. Both paths must produce identical outputs.

## Key design patterns

**Configuration-driven pipeline**: `task_definitions.json` defines which scripts map to which analyses and their default parameters. Its `default_tasks` block lists the analyses auto-added to a new session (e.g. a BROAD blast searching all species) — entries can set `requires_organism: true` to be added only once a species is selected, and `taxid_from_rank` (e.g. `"order"`) to auto-fill their `taxid` from that rank in the selected organism's lineage. `config.py` defines available species (with NCBI taxonomy IDs) and which result types are continuous vs. range-based.
**Shiny reactive state**: analysis parameters are stored in a shared reactive dict (`shared_dict`) passed between modules. `utils/helpers.py:update_shared_dict()` handles updates.
**Async job submission**: `run_tag_sites_from_json.py` spawns analysis scripts as subprocesses and polls for completion, enabling parallel execution of independent analyses.
**Scoring never withholds reagents**: guide and site scores (RS3, isoform/topology
restrictions, conservation) are display-and-ranking signals only. They must never filter a
guide or site out of the reagent design UI — annotate instead, so the user always sees every
option and decides. Guides are ordered by distance from cut to insertion site; RS3 rides
along as a badge and never participates in selection or sorting. A regression check for this
is that the pre-existing columns of `{run}_reagents.tsv` are byte-identical with and without
`rs3` installed.

**Never call a blocking function directly from a Shiny effect**: Shiny's asyncio event
loop serves *every* session, so one slow `requests.get` freezes the whole app — including
its ability to notice a pasted sequence or enable a button. Route remote work through
`modules/setup_server.py:_off_loop()` (`asyncio.to_thread`, plus an optional
`asyncio.wait_for` budget) or `@reactive.extended_task` + `run_in_executor`
(`modules/progress_server.py:293`). Note also that `requests`' `timeout=` is a per-socket
read timeout, not a deadline: a response arriving a byte at a time never trips it at any
setting. `scripts/http_retry.py` enforces a real wall-clock budget. See issue #64.

**Gene model must translate to the input protein**: residue numbers in the reagent table come
from translating the gene model's CDS, so an exon-structure error shifts every downstream site
(DBL-1: Q239 reported as F239). `design_tag_reagents` therefore runs `cds_check.check_model`
against the input protein (`--protein_fasta`, also accepted as `--input_file`): a length
difference (any alignment gap) raises with an indel summary and points to the GenBank upload;
substitutions alone warn (positions listed in the log, the sidecar and the Reagents tab) and
reagents are still designed, unless they exceed 2% of the protein length, which also fails. The cause of the DBL-1 error: EBI Genewise defaults to flat GT/AG
splicing and reads short worm introns through as coding sequence; `genewise_remote.py` now
sends `splice=model, init=global` (sweep of 11 genes x 17 settings; 10/11 exact, vs 3/11 at the
default; EBI exposes no gap penalties). A GenBank with CDS/exon features uploaded as the genomic
region replaces Genewise entirely (`genbank_input.py`) and is held to the same check.
With `backends.genewise` set to `bulk` (a local run) there is no Genewise and no genomic upload:
`genome_regions.find_transcripts_by_sequence` matches the input protein to the GFF3 transcript whose
CDS translates to it exactly, and `genewise_bulk.py` cuts that gene's region from the local genome.

**Alignment rendering**: sequence alignments are pre-rendered as matplotlib SVGs (stored in a dict keyed by alignment name) and displayed statically — this replaced an earlier real-time Plotly approach that was too slow.

## Environment

Key dependencies: `shiny`, `biopython`, `plotly`, `pandas`, `numpy`, `scipy`, `matplotlib`,
`requests`, `primer3-py`, `pyhmmer` (local Pfam scanning), `reportlab` (alignment PDFs), and
`rs3` + `lightgbm` (RS3 guide scoring). `diamond` and `mafft` are used only by the batch-mode
local backends. No local bioinformatics binaries are required for the interactive app —
Clustal Omega, BLAST, Genewise, and InterPro all run via the EBI REST API.

### ⚠️ Version pinning — read before adding any dependency

**The environment is deliberately held at Python 3.10 with numpy 1.x.** `environment.yml` pins:

| Pin | Why |
|---|---|
| `python=3.10` | required by the caps below |
| `numpy=1.26.*` | `rs3` requires `numpy<=1.26.4` |
| `scikit-learn=1.0.2` | `rs3` requires `scikit-learn<=1.0.2` |
| `lightgbm=3.3.5` | `rs3` requires `lightgbm<=3.3.5` |

The **only** reason for all four is the `rs3` package (Rule Set 3 guide scoring,
`scripts/guide_efficiency.py`). rs3 0.0.18 is the latest release and has not been updated
since Feb 2024. Nothing else in TAGSITES needs these versions — note in particular that
**`scikit-learn` is imported nowhere in this codebase**; it is present only to satisfy
rs3/lightgbm. scikit-learn and lightgbm come from conda-forge rather than pip because the
pinned old versions have no osx-arm64 wheels.

**If a new dependency needs `numpy>=2` or a newer Python, drop RS3 rather than fighting the
pins.** RS3 is designed to fail soft: `guide_efficiency.load_rs3()` catches broad `Exception`,
so with `rs3` absent the pipeline logs "RS3 unavailable", writes blank
`rs3_score`/`rs3_percentile` columns, and leaves every other column untouched. The UI renders
"RS3 —". Reverting costs only the score. See commit 929fe8e for the full rationale and the
evidence behind this tradeoff.

## Deploying to shinyapps.io

```bash
./deploy.sh              # deploy (updates the existing app in place)
./deploy.sh --manifest   # dry run: list what would be bundled, upload nothing
```

**Always deploy via `deploy.sh`, never a bare `rsconnect deploy`.** `.rscignore` is an
R-only feature that rsconnect-python does not read, so every exclusion lives as a
`-x` flag in that script. A bare deploy bundles ~10 GB (`data/` is 9.7 GB and
`webservice-clients/` another 188 MB); `deploy.sh` brings it to 181 files / 3.8 MB.
`requirements.txt` is the dependency manifest — `environment.yml` is excluded.

Note rsconnect only ignores `__pycache__/` at the top level, so nested ones need
`-x '**/*.pyc'`; the `**/__pycache__/` form silently fails to match.

## External services

- **NCBI BLAST** (via EBI REST API) — `scripts/site_selection_util.py:ncbiblast_call()`
- **EBI InterPro API** — `scripts/call_interpro.py`
- **AlphaFold DB** — `scripts/existing_AF_model.py`
- **UniProt REST** (entry fetch + checksum lookup) — `scripts/uniprot_api.py`, used by `scripts/uniprot_features.py` and `scripts/existing_AF_model.py`
- **UniProt taxonomy** — `uniprot_species.flat.txt` (local flat file, ~1.7MB)
