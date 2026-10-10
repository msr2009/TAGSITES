# Running TAGSITES locally (single-protein analysis)

This guide sets up the Shiny app so that analyses run on the local machine instead of through
EBI, UCSC and Ensembl web services. It covers one-protein-at-a-time use in the browser. The
default install in `README.md` needs no local tools and sends every analysis to a web service;
follow this guide only to move those calls onto local hardware.

The walkthrough uses *C. elegans*. Notes for other organisms are in [step 9](#9-other-organisms).

## What "local" covers

The app reads a backend mode per analysis from a config file (`scripts/providers.py`). Every
analysis script resolves its backend at run time, so the same Shiny app runs local or remote
depending only on which config file it is started with.

| Analysis | Local mode | Replaces | Reference data needed |
|---|---|---|---|
| Domains | `scan` | EBI InterProScan | Pfam-A HMMs (~2.5 GB) |
| Structure (pLDDT/SASA source) | `bulk` | UniProt checksum lookup, EBI BLAST, AFDB fetch | AFDB proteome PDBs (~3.8 GB) |
| Gene model + genomic region | `bulk` (under `genewise`) | EBI Genewise and the genomic-region upload | Ensembl genome + GFF3 (~60 MB) |
| Conservation (BLAST tasks) | `local` | EBI BLAST + Clustal Omega | Swiss-Prot + Rhabditida TrEMBL DIAMOND databases (~1.8 GB) |
| Off-target / primer screens | `local` | UCSC BLAT / EBI blastn | BLAST+ genome database (~55 MB) |

The `genewise: bulk` mode runs no Genewise. CDS exons are read from the GFF3, cut from the
genome, and translated; the transcript whose translation equals the input protein is the gene
model, and its genomic window (plus flank) is the region. No genomic file is uploaded.

Already local in every mode: pLDDT/SASA extraction from a PDB, PTM regex sites, sliding-window
property scores, and RS3 guide scoring.

Features that still call web services even in local mode are listed in [step 8](#8-calls-that-remain-remote).

## 1. Prerequisites

- macOS or Linux, with [mamba](https://mamba.readthedocs.io/) (or conda).
- About 9 GB of free disk for reference data. Put it on a local disk; SQLite indexes are
  unreliable on network volumes.
- Network access for the one-time downloads in step 3. Nothing else needs the network once the
  remaining-remote features in step 8 are avoided.

## 2. Create the environment

```bash
git clone https://github.com/msr2009/TAGSITES.git
cd TAGSITES
mamba env create -f environment.yml
mamba activate tagsites
```

Use `environment.yml`, not `requirements.txt`. The requirements file is the shinyapps.io
manifest and omits the local tools (`diamond`, `mafft`, `blast`, `pyhmmer`). Python 3.10 and
numpy 1.x are deliberate pins; see "Version pinning" in `CLAUDE.md`.

Check that the tools are on the path:

```bash
which diamond mafft blastn makeblastdb blastdbcmd
python -c "import pyhmmer; print(pyhmmer.__version__)"
```

Every command should print a path or version. A missing `blastn` does not break the app (off-target
screens fall back to BLAT/EBI), but the run is then no longer fully local.

## 3. Download reference data

All downloads go to `data/reference/` (set by `reference_data.out_dir` in `batch.config.json`).
They are resumable, and a step whose output already exists is skipped; add `--force` to redo one.

Run the steps in this order. Each command can be run on its own, so large steps can be done
later if only some backends are wanted.

```bash
# Structure + gene-model lookup tables (UniProt proteome, AlphaFold PDBs)
python scripts/reference_data.py --only uniprot,afdb
python scripts/local_store.py --build

# Gene model (Ensembl genome + GFF3), then its SQLite index (exons + a translation of every CDS)
python scripts/reference_data.py --only genome,gff3
python scripts/genome_regions.py --build

# Domains
python scripts/reference_data.py --only pfam

# Conservation (largest download: Rhabditida TrEMBL is ~2 million sequences)
python scripts/reference_data.py --only swissprot,rhabditida_trembl

# Off-target / primer screens (needs the genome step above)
python scripts/reference_data.py --only blastdb
```

| Step | Size on disk | Enables |
|---|---|---|
| `uniprot` + `local_store.py --build` | ~55 MB | accession/checksum lookup for `structure: bulk` |
| `afdb` | ~3.7 GB (2.6 GB tar + extracted PDBs) | `structure: bulk` |
| `genome`, `gff3` + `genome_regions.py --build` | ~55 MB | `genewise: bulk` (gene model and region; needs nothing else) |
| `pfam` | ~2.5 GB | `domains: scan` |
| `swissprot`, `rhabditida_trembl` | ~1.8 GB | `conservation: local` |
| `blastdb` | ~25 MB | local off-target screens |

`--stats` on `local_store.py` and `genome_regions.py` prints row counts and confirms that an
index was built. For `genome_regions.py`, the "translated transcripts" count should equal the
transcript count. An index built earlier gains the translation table when `--build` is rerun;
no `--force` is needed.

Not needed for single-protein use: the `interpro` step (a 13 GB stream, only for the
`domains: bulk` mode), the WormBase protein FASTA, the topology cache and DeepTMHMM. Those serve
the proteome-scale pipeline.

## 4. Create the backend config

`batch.config.json` is the checked-in default. Do not edit it in place: it is read for every
caller, so changes would silently alter other runs. Create a separate file instead.
`batch.config.local.json` is gitignored, so each machine keeps its own.

```json
{
  "backends": {
    "domains": "scan",
    "structure": "bulk",
    "genewise": "bulk",
    "conservation": "local",
    "offtarget_region": "local"
  },
  "conservation_local": {
    "search_databases": [
      { "name": "swissprot", "evalue": 1e-10, "max_target_seqs": 50 },
      { "name": "rhabditida_trembl", "evalue": 1e-5, "max_target_seqs": 100 }
    ],
    "threads": 8,
    "render_alignment_pdf": true
  }
}
```

Notes on the file:

- Keys missing here take the `batch.config.json` defaults.
- `threads` is the DIAMOND thread count. Set it to the number of cores to spare.
- `render_alignment_pdf: true` produces the alignment view that the Results tab shows. The
  batch-oriented local config turns it off; keep it on for the app.
- `topology` is omitted on purpose. The app has no topology runner, so DeepTMHMM does not run
  from the UI.
- Any single analysis can stay `"remote"`. Mixed setups are valid.

## 5. Optional cache builds

Two caches can speed up repeated work. Neither is required for single-protein runs.

- **Pfam cache** (`python scripts/build_pfam_cache.py`): pre-scans a whole proteome. Without it,
  `domains: scan` scans the submitted sequence on demand.
- **Conservation pre-search** (`python scripts/conservation_presearch.py`): batches DIAMOND
  searches. Without it, each protein runs its own search, taking roughly a minute because the
  databases load each time.

Both build steps need the WormBase protein FASTA (`reference_data.protein_fasta`), which
WormBase's portal blocks for scripted downloads. Download it in a browser and place it in
`data/reference/` if either cache is wanted.

## 6. Start the app

```bash
mamba activate tagsites
TAGSITES_BATCH_CONFIG=batch.config.local.json python app_modular.py
```

Open http://127.0.0.1:8000.

`TAGSITES_BATCH_CONFIG` takes a path, relative to the repository root or absolute. The config is
read once per process, so restart the app after any edit to it.

To confirm that the local modes are active before opening the browser:

```bash
TAGSITES_BATCH_CONFIG=batch.config.local.json python -c "
import sys; sys.path.insert(0, 'scripts')
import providers
for a in ['domains', 'structure', 'genewise', 'conservation', 'offtarget_region']:
    print(a, providers.backend_mode(a))"
```

The printed modes should match step 4. The off-target backend actually used is recorded as
`_meta.backend` (`local` or `remote`) in the run's off-target sidecar JSON next to the reagents
output.

## 7. Run one protein

1. **Setup tab, sequence.** Paste a protein FASTA, or upload a FASTA or an AlphaFold PDB.
   `structure: bulk` finds a model only for canonical *C. elegans* proteins whose sequence
   matches the proteome exactly; any other sequence needs an uploaded PDB.
2. **Setup tab, analyses.** The default tasks are pre-filled. The BROAD and NARROW BLAST tasks
   run the same local search, because `conservation: local` ignores the task's taxid and
   database. Keep one of them.
3. **Genomic region.** Upload nothing. With `genewise: bulk`, the Reagents task is created
   without a genomic file: the gene model and region come from the local GFF3 and genome, with
   flank size set by `batch_run.genomic_flank_bp` (default 2000). The app confirms this with a
   notification on save. An uploaded GenBank file with CDS/exon features still takes precedence
   as the gene model.
4. **E-mail.** The app still requires an address to save an analysis, although no EBI job is
   submitted.
5. **Submit.** The Progress tab shows each task. After completion, the Results tab shows plots
   and alignments, and the Reagents tab lists guides and homology arms.

Expect the first run to be slower: Pfam models and DIAMOND databases load from disk, and the
conservation step is the longest task.

## 8. Calls that remain remote

No local backend exists for the following, so they still use the network in local mode.

| Feature | Service | Workaround |
|---|---|---|
| UNIPROT features task (default task) | rest.uniprot.org | remove the task from the session |
| Gene search and isoform expansion on the Setup tab | UniProt | paste the sequence directly |
| Organism / taxonomy search on the Setup tab | UniProt | not needed for local conservation, which ignores organism |
| AlphaFold PDB download after picking a search hit | EBI dbfetch | rely on `structure: bulk`, or upload a PDB |
| Genomic "Fetch" button | Ensembl | not needed with `genewise: bulk` |
| Off-target screens, when `blastn` or the database is missing | UCSC BLAT, then EBI blastn | install BLAST+ and run the `blastdb` step |

## 9. Other organisms

Organism-agnostic backends work with any sequence:

- `domains: scan` (Pfam covers all organisms).
- `conservation: local`, with Swiss-Prot as the universal target. The second target,
  `rhabditida_trembl`, is the nematode clade; swap it for a relevant TrEMBL clade by changing
  `reference_data.rhabditida_taxid` and the `search_databases` entries.

Worm-specific backends need per-organism data:

- `structure: bulk` uses one AlphaFold proteome tarball (`reference_data.afdb_tarball_url`) and one
  UniProt proteome (`uniprot_proteome_id`).
- `genewise: bulk` uses the Ensembl genome and GFF3 (`ensembl_release`, `ensembl_species`,
  `ensembl_assembly`) and the WormBase cross-references in the UniProt entries, which other
  organisms do not carry in the same form.
- `offtarget_region: local` uses a genome database registered under `local.genomes` in
  `offtarget.config.json`, keyed by NCBI taxid.

For another organism, the safest route is to keep `domains` and `conservation` local, and leave
`structure`, `genewise` and `offtarget_region` at `"remote"`.

## 10. Failure modes and troubleshooting

- **Backend change has no effect.** The config is cached per process. Restart the app.
- **"no model" from the structure task.** The sequence is not a canonical proteome entry.
  Upload a PDB, or set `structure` back to `"remote"`.
- **`LookupError` in the Reagents task.** No annotated transcript translates exactly to the input
  protein (a different strain or species, an unannotated isoform, or a mitochondrial protein).
  Upload a GenBank file with the gene model, or set `genewise` back to `"remote"` and upload a
  genomic region.
- **Off-target backend shows `remote`.** `blastn` or the database is missing, so the screen fell
  back to BLAT/EBI. Re-check step 2 and rerun `reference_data.py --only blastdb`.
- **Slow conservation.** Each query searches two databases from scratch. Raise
  `conservation_local.threads`, or build the pre-search cache (step 5).
- **Pfam settings.** `scripts/pfam_scan.py` reads `domains_scan` from the checked-in
  `batch.config.json` and ignores the override file. Edit that block in `batch.config.json`
  to change scan settings.
- **No local fallback.** `structure: bulk`, `domains: bulk` and `genewise: bulk` never fall back
  to their remote versions automatically. Switch the individual backend in the config instead.
