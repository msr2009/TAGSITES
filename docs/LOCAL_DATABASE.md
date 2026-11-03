# Running TAGSITES with the local precomputed database

The proteome database (`proteome.sqlite3`) holds precomputed analyses and reagent designs for
every protein of one organism (C. elegans, or yeast). Two apps read it:

| App | Command | Port | What the database does |
|---|---|---|---|
| Main app | `python app_modular.py` | 8000 | Setup recognises a loaded sequence and offers **Load precomputed results** |
| Batch review app | `python app_batch.py` | 8001 | Setup resolves a protein list against the database; Review designs tags across many proteins |

Both are local only. The deployed shinyapps.io app has no database file, so the card and the
batch app never appear there.

Installing the code, creating the environment and the single-protein workflow are covered in
the main [README](../README.md). Local analysis backends (not needed just to read the
database) are covered in [LOCAL_SETUP.md](LOCAL_SETUP.md).

---

## Getting the database

The database file is not in the repository. Either:

- **Request it from the maintainer** (Matt Rich). The package is the database plus the batch run
  folder (the AlphaFold structures and conservation alignments of every protein), both required.
  This is the simplest route.
- **Build it locally** for any organism by following [BUILD_LOCAL_DATABASE.md](BUILD_LOCAL_DATABASE.md)
  (config, reference downloads, indexes, per-protein analyses, then the database build).
  [LOCAL_SETUP.md](LOCAL_SETUP.md) covers the environment and local backends it builds on.

Once the file exists, continue below.

---

## 1. What is needed

| Item | Needed for | Config key |
|---|---|---|
| `proteome.sqlite3` | everything below | `proteome_db.path` |
| Genome FASTA, GFF3 and `genome_regions.sqlite3` | regenerating the ~1 kb homology arms (the database stores none) | `reference_data.out_dir` |
| **Required:** the batch run folder (`proteome_v1/`, ~12 GB, one subfolder per protein) | the AlphaFold structure and conservation alignment of every loaded run (Results needs both) | `proteome_db.run_dir` |
| `blastn` plus a local BLAST database | fast off-target screens (optional; falls back to UCSC BLAT, then EBI) | `backends.offtarget_region` |

Get the genome files for the organism with:

```bash
conda activate tagsites
python scripts/reference_data.py --only genome,gff3
python scripts/genome_regions.py --build
```

The run folder must sit beside `proteome.sqlite3` and `run_dir` must point at it. Without it,
structures and alignments are reported as "not in the batch folder" and are absent from Results,
so an install is incomplete until the folder is in place.

The database is a read-only file opened with SQLite's `immutable` mode, so it can sit on a
network volume. Do not edit it while an app is running.

---

## 2. Point the config at the database

Every path comes from one JSON config. The apps read `batch.config.json` unless the
`TAGSITES_BATCH_CONFIG` environment variable names another file (relative paths resolve against
the repo root). `app_batch.py` defaults that variable to `batch.config.worm_local.json`.

Create a config for the organism, or edit `batch.config.json`. The blocks that matter:

```json
{
  "reference_data": {
    "out_dir": "/path/to/reference",
    "ensembl_species": "caenorhabditis_elegans",
    "ensembl_assembly": "WBcel235",
    "ensembl_release": "114",
    "uniprot_proteome_id": "UP000001940",
    "protein_fasta": "/path/to/reference/c_elegans.canonical_bioproject.current.protein.fa.gz"
  },
  "proteome_db": {
    "path": "/path/to/proteome.sqlite3",
    "app_match": true,
    "run_dir": "/path/to/runs/proteome_v1"
  },
  "batch_run": { "taxid": "6239" },
  "batch_app": {
    "batch_dir": "/path/to/batches",
    "workers": 6,
    "cache_entries": 1500,
    "port": 8001
  }
}
```

| Key | Meaning |
|---|---|
| `proteome_db.path` | The database file. Relative paths resolve against the repo root. |
| `proteome_db.app_match` | `true` lets the main app look a loaded sequence up in the database. Inert when the file is missing. |
| `proteome_db.run_dir` | The `proteome_run.py` output folder the database was built from. |
| `reference_data.out_dir` | Folder with the genome, GFF3 and `genome_regions.sqlite3`. Must belong to the same organism and assembly as the database. |
| `reference_data.ensembl_*` | Select the genome and GFF3 file names inside `out_dir`. |
| `batch_run.taxid` | Taxonomy ID written into loaded runs (default 6239, C. elegans). Set it for any other organism. |
| `batch_app.batch_dir` | Local-disk folder where each batch saves its files. |
| `batch_app.cache_entries` | Parsed runs kept in memory for Review (about 3 per protein). |
| `batch_app.port` | Port of the batch app. |

Keep `reference_data` and `proteome_db` in the same config. Pairing a C. elegans database with a
yeast genome gives wrong arms, because both are read from this one file.

For yeast, start from `batch.config.yeast.json` (UP000002311, R64-1-1, `taxid` 559292) and
leave `proteome_db.app_match` as the maintainer set it.

---

## 3. Check that the database is found

```bash
export TAGSITES_BATCH_CONFIG=batch.config.worm_local.json   # or your own config
python scripts/proteome_db.py lookup trxr-1
```

A table of the protein's products means the path resolves. To print the summary the app shows
(size, build date, protein and gene counts):

```bash
python -c "import sys; sys.path.insert(0, 'scripts'); import proteome_db; print(proteome_db.describe())"
```

An error saying `... not found; build it with scripts/build_proteome_db.py` means
`proteome_db.path` is wrong or the config in effect is not the one edited. Another config in
the environment variable wins over `batch.config.json`.

---

## 4. Main app with precomputed results

```bash
conda activate tagsites
export TAGSITES_BATCH_CONFIG=batch.config.worm_local.json   # omit to use batch.config.json
python app_modular.py
```

Open **http://127.0.0.1:8000**. On Setup, load a protein (UniProt hit, upload or paste). When
the sequence matches a database protein exactly, a **Load precomputed results** card appears.

Loading writes the run's `.run.json` at once, then exports in the background (3 to 10 s, mostly
arm regeneration from the genome). Progress shows two tasks, Analyses and Reagents, with Run
disabled. Results and Reagents fill in as each finishes.

Things to know:
- A track the batch never ingested (for example modifications whose task did not finish `ok`)
  is absent.
- Site scores are recomputed by the app with the reagent rows, so they can differ slightly from
  the stored `site_scores`.
- No card means `app_match` is off, the file is missing, or the sequence differs from the
  analysed one (check by pasting the exact database sequence).

---

## 5. Batch review app

```bash
conda activate tagsites
python app_batch.py                                          # default: batch.config.worm_local.json
TAGSITES_BATCH_CONFIG=batch.config.yeast.json python app_batch.py   # another organism
```

Open **http://127.0.0.1:8001**. The app can run beside the main app.

1. **Setup**: panel 0 shows the database in use (path, size, protein count). **Use database**
   there points the process at another `.sqlite3` file (checked before use). Paste or upload a
   protein list. A gene with several products gets a dropdown. The default product is the id
   equal to the search term, else the UniProt entry, then reagents designed, then structure, then
   the longest. Set the site filters (suggested sites per protein, position and annotation
   masks) as batch defaults.
2. **Prepare**: exports each protein to `{batch_dir}/{batch}/{key}/` and fills the in-memory
   cache. About 10 s for 10 proteins. Alignment PDFs render afterwards in a background process.
3. **Review**: the single-protein Results view plus prev/next, a status-coloured list and a
   per-site guide picker. **Accept → Next** accepts and moves on. **Skip** marks a protein
   skipped; it is listed in `{batch}.skipped.txt` and never reaches Output.
4. **Output**: a ZIP of reagents for accepted proteins, plus combined TSVs.

A batch saves and loads as one ZIP. It records the database it used, and reopening it switches
back to that file.

---

## 6. Troubleshooting

| Symptom | Likely cause |
|---|---|
| Panel 0 shows no database | `proteome_db.path` wrong, or the wrong config is active (`echo $TAGSITES_BATCH_CONFIG`) |
| No "Load precomputed results" card | `app_match` false, file missing, or sequence not an exact match |
| Reagents tab errors on arms | `reference_data.out_dir` lacks the genome, GFF3 or `genome_regions.sqlite3`, or is another organism |
| Structure or alignment "not in the batch folder" | `proteome_db.run_dir` wrong or not mounted |
| Slow first protein switch | one-time imports and genome load (about 2 s); later switches take 25 to 500 ms |
| Off-target screen slower than expected | no local `blastn` or BLAST database; the app falls back to BLAT/EBI |
