# Command-line use

Every script in `scripts/` runs without the Shiny app. Use this guide to run the pipeline from a
shell, script many proteins, or produce reagents from a finished run.

Before starting: [INSTALL.md](INSTALL.md). For local backends, also set
`TAGSITES_BATCH_CONFIG` ([CONFIG.md](CONFIG.md)); the CLI reads the same config as the app.

```bash
mamba activate tagsites
```

## 1. One protein from a run JSON

A run JSON lists the sequence, the tasks and their parameters. The app writes one
(`{run}.run.json`) when an analysis is saved; one can also be written by hand. A minimal example
that uses precomputed Genewise output, so it needs no network, is `tests/data/example_run.json`.

```bash
python scripts/run_tag_sites_from_json.py -i path/to/run.json
python scripts/run_tag_sites_from_json.py -i path/to/run.json --force   # redo finished tasks
```

Outputs go to the run JSON's `global.working_dir`. Task status is tracked in `{run}.status.json`,
so an interrupted run resumes without resubmitting EBI jobs. Independent tasks run in parallel.
Set the e-mail in `global.email`.

## 2. Many proteins from a manifest (web services)

`scripts/batch_run_tag_sites.py` builds a run JSON per protein from a manifest and a shared task
template, then runs them with bounded concurrency so EBI is not flooded.

Manifest (CSV or TSV), one row per protein: `run_name`, `input_file`, and optional `pdb`,
`genomic_file`, `working_dir`, `email`. An example is `tests/data/batch_manifest_example.tsv`.

```bash
python scripts/batch_run_tag_sites.py \
    --manifest proteins.tsv \
    --params params/worm_default.json \
    --email you@example.edu \
    --output-root ./batch_runs \
    --max-concurrent 3
```

Other flags: `--taxid`, `--flank-bp` (genomic flank, default 2000), `--force`. Rerunning skips
proteins and tasks that already succeeded.

This is separate from the *batch review app*, which works from a precomputed database
([BATCH_APP.md](BATCH_APP.md)).

## 3. Reagents for chosen sites

`scripts/reagent_bundle.py` produces the same download the Reagents tab does, from a finished
run that has a reagents table.

```bash
python scripts/reagent_bundle.py \
    --run_json data/x/x.run.json \
    --sites 45,120 \
    --insert FLAG \
    --repair_type ssodn \
    --out x_reagents.zip
```

`--insert` takes a name from the tags table or a raw DNA sequence. `--repair_type` is `ha`,
`ssodn`, `amplicon` or `plasmid`; `--arm_length` and `--prefix` adjust arm size and naming. The
ZIP holds the guide, oligo, genotyping and (for SapTrap) plasmid TSVs plus annotated GenBank
records.

## 4. Proteome-scale analysis and the database

These need local backends and reference data. The full procedure, with sizes and checks, is in
[BUILD_LOCAL_DATABASE.md](BUILD_LOCAL_DATABASE.md).

```bash
python scripts/reference_data.py --only uniprot,afdb,genome,gff3   # downloads
python scripts/proteome_run.py --out-dir <run_dir> --limit 20      # analyses, one folder per protein
python scripts/build_proteome_db.py --estimate                      # size and time projection
python scripts/build_proteome_db.py --stage ingest --yes            # build
```

## 5. Query a finished database

```bash
python scripts/proteome_db.py lookup trxr-1                         # products of a gene
python scripts/proteome_db.py sql "SELECT track, COUNT(*) FROM features GROUP BY track"
python scripts/proteome_db.py export Q17745 /path/to/run_folder     # one protein as an ordinary run folder
python scripts/proteome_db.py --db /path/to/other.sqlite3 lookup unc-54
```

The database path comes from `proteome_db.path` in the active config unless `--db` is given.
The exported folder opens in the main app and works with `reagent_bundle.py`.

## Troubleshooting

| Symptom | Cause |
|---|---|
| Local backends ignored | `TAGSITES_BATCH_CONFIG` not exported in this shell. |
| EBI rejects jobs | Missing or invalid `global.email`. |
| `LookupError` from the reagents task | No gene model translates to the protein; supply a GenBank file or a genomic region ([LOCAL_SETUP.md](LOCAL_SETUP.md#troubleshooting)). |
