# Configuration

Every local workflow (local backends, the precomputed database, the batch app, database builds)
is controlled by one JSON file. This page explains it once; the other guides link here.

## Which file is read

| Situation | File read |
|---|---|
| `TAGSITES_BATCH_CONFIG` unset | `batch.config.json` (checked in; all analyses use web services) |
| `TAGSITES_BATCH_CONFIG=<path>` | that file. A relative path resolves against the repository root. |
| `python app_batch.py`, variable unset | `batch.config.worm_local.json` (the app sets this default; the file is gitignored, so create it or export the variable) |
| Hosted app (shinyapps.io) | `batch.config.remote.json`, bundled as `batch.config.json` by `deploy.sh` |

Rules:

- **Never edit `batch.config.json` for local work.** Every caller reads it, so a change silently
  alters other runs. Use a separate file.
- The config is read **once per process**. Restart the app after editing.
- The active file is used on its own: it is *not* merged over `batch.config.json`. A missing
  key falls back to the reader's built-in default (a missing `backends` entry means `remote`),
  so start from the template below rather than a partial file.
- Use **absolute paths**. A leading `~` is not expanded in `out_dir`, `path` or `run_dir`.
- One config describes **one organism**. Use a separate config and a separate `out_dir` for each.
  Pairing a database with another organism's genome gives wrong homology arms.

## Create a local config

```bash
cp batch.config.example.json batch.config.local.json
# replace the @DATA_DIR@ placeholder with an absolute folder on local disk
sed -i.bak "s|@DATA_DIR@|$HOME/tagsites_local/worm|g" batch.config.local.json && rm batch.config.local.json.bak
export TAGSITES_BATCH_CONFIG=batch.config.local.json
```

`batch.config.local.json` and `batch.config.yeast*.json` are gitignored, so machine paths stay
out of the repository. The template holds *C. elegans* values; the organism table in
[BUILD_LOCAL_DATABASE.md](BUILD_LOCAL_DATABASE.md#1-write-the-organisms-config) lists what to
change for another species.

Check which file and which backends are in effect:

```bash
python -c "
import sys; sys.path.insert(0, 'scripts')
import providers
print(providers._config_path())
for a in ['domains', 'structure', 'genewise', 'conservation', 'topology', 'offtarget_region']:
    print(a, providers.backend_mode(a))"
```

## Blocks

| Block | Used by | Key settings |
|---|---|---|
| `backends` | every analysis | Per analysis, `"remote"` (default, EBI/UCSC/Ensembl) or a local mode: `domains: scan` (Pfam via pyhmmer; `bulk` for a prebuilt InterPro table), `structure: bulk`, `genewise: bulk` (gene model from local GFF3 and genome), `conservation: local` (DIAMOND + MAFFT), `topology: deeptmhmm`, `offtarget_region: local` (BLAST+). Any single analysis may stay remote. |
| `reference_data` | `scripts/reference_data.py`, local backends, arm regeneration | `out_dir` (download folder), UniProt proteome id, AlphaFold tarball URL, Ensembl release/species/assembly, `protein_fasta`, TrEMBL clade (`trembl_taxid`, `trembl_name`), UCSC genome and phastCons track, `ucsc_chrom_alias` |
| `conservation_local` | `conservation: local` | `search_databases` (one DIAMOND database per entry, with `evalue` and `max_target_seqs`), `threads`, `render_alignment_pdf` (keep `true` for the apps) |
| `domains_scan` | `domains: scan` | `bit_cutoffs`, `clan_filter`, `cache_profiles`, `cpus` (0 = all cores). Read from the checked-in `batch.config.json`, not the override file (see [LOCAL_SETUP.md](LOCAL_SETUP.md#troubleshooting)). |
| `deeptmhmm` | database builds | `install_dir`, `python` (a separate licensed install), batch sizes |
| `proteome_db` | database build and both apps | `path` (the SQLite file), `app_match` (offer precomputed results in the main app), `run_dir` (per-protein run folders), `build_dir` (local disk for assembly), `scratch_dir`, `reagent_chunk`, optional `curated_sites` |
| `batch_run` | `proteome_run.py`, the apps | `taxid`, `max_workers`, `max_concurrent_network_calls`, `genomic_flank_bp` |
| `batch_app` | `app_batch.py` | `batch_dir` (local disk), `workers`, `ebi_max_concurrent`, `cache_entries`, `port` |

## Other config files

| File | Purpose |
|---|---|
| `offtarget.config.json` | Off-target and primer screens: BLAT assemblies by taxid, post-filter thresholds, `local.genomes` (taxid → BLAST database and FASTA). A new organism needs an entry here. |
| `ucsc.local.json` | Gitignored UCSC BLAT API key. |
| `ensembl.config.json` | Genomic auto-fetch sources and order (NCBI, UCSC, Ensembl). |
| `scores.config.json` | Site-scoring weights and masks. |
| `task_definitions.json` | Task registry and `default_tasks`; see [DATA_MODEL.md](DATA_MODEL.md). |
| `params/*.json` | Saved analysis presets (for example `worm_default.json`). |

## Which blocks each workflow needs

| Workflow | Blocks to set |
|---|---|
| Local backends, one protein | `backends`, `reference_data`, `conservation_local` |
| Use a database in the main app | `proteome_db`, `reference_data` (genome files for arms), `batch_run.taxid` |
| Batch app | the database blocks plus `batch_app` |
| Build a database | all of the above plus `deeptmhmm` (optional) and `proteome_db.build_dir` |
