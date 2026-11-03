# Building a local proteome database

Use this to build `proteome.sqlite3` for one organism: precomputed analyses, tag-site scores and
CRISPR reagent designs for every protein, in one file the apps read (see
[LOCAL_DATABASE.md](LOCAL_DATABASE.md) for using it). The steps are the same for any organism
with a UniProt reference proteome, an Ensembl (or compatible) genome and GFF3, and an AlphaFold
DB proteome. C. elegans (`batch.config.example.json`) and *S. cerevisiae*
(`batch.config.yeast.json`) are the worked examples, and the example values below are C. elegans
unless marked.

Allow disk, time and care: the reference data for one organism is tens of GB (the shared
Swiss-Prot / TrEMBL / Pfam files dominate) and a full run takes hours. Sizes in this guide are
approximate; `--estimate` projects the database size for the organism at hand. Run `--estimate` before every build (step 8).

Before starting:

- The environment from [INSTALL.md](INSTALL.md), with `diamond`, `mafft` and BLAST+ on the path
  (`which diamond mafft blastn makeblastdb`).
- How configs work: [CONFIG.md](CONFIG.md). Each local backend is described in
  [LOCAL_SETUP.md](LOCAL_SETUP.md).
- A machine with many cores, tens of GB of local disk, and hours of run time. A small organism
  such as yeast is a good first build.

---

## Pipeline overview

```
 1 config  →  2 downloads  →  3 indexes  →  4 caches  →  5 analyses  →  6 database
 (organism)   reference_data   local_store    Pfam,         proteome_run    build_proteome_db
                               genome_regions topology      (per protein)   ingest, scores,
                                                                            reagents
```

---

## 1. Write the organism's config

Everything organism-specific lives in one JSON file chosen with `TAGSITES_BATCH_CONFIG`
([CONFIG.md](CONFIG.md)). Copy the template and change the values below.
`batch.config.local.json` and `batch.config.yeast*.json` are gitignored, so machine paths stay
out of the repository. A config for another organism can use any name ending in `.json`, but
only those patterns stay untracked unless `.gitignore` is extended.

```bash
cp batch.config.example.json batch.config.local.json      # C. elegans values
sed -i.bak "s|@DATA_DIR@|$HOME/tagsites_local/worm|g" batch.config.local.json && rm batch.config.local.json.bak
export TAGSITES_BATCH_CONFIG=batch.config.local.json
```

| Block / key | What to set | C. elegans | *S. cerevisiae* |
|---|---|---|---|
| `reference_data.out_dir` | Folder for all downloads (large, ideally local disk) | `data/reference` | `~/tagsites_local/sacCer3` |
| `reference_data.uniprot_proteome_id` | UniProt reference proteome | `UP000001940` | `UP000002311` |
| `reference_data.afdb_tarball_url` | AlphaFold DB proteome tarball for that organism | `…/UP000001940_6239_CAEEL_v6.tar` | `…/UP000002311_559292_YEAST_v6.tar` |
| `reference_data.ensembl_species`, `ensembl_assembly`, `ensembl_release` | Select the genome FASTA and GFF3 file names | `caenorhabditis_elegans`, `WBcel235`, `114` | `saccharomyces_cerevisiae`, `R64-1-1`, `114` |
| `reference_data.protein_fasta` | Protein FASTA of the same annotation (WormBase or Ensembl `pep.all`) | WormBase protein FASTA | Ensembl `pep.all.fa.gz` |
| `reference_data.trembl_taxid`, `trembl_name` | Clade for the TrEMBL ortholog database and its file name | `6236`, `rhabditida_trembl` | `4891`, `saccharomycetales_trembl` |
| `reference_data.ucsc_genome`, `ucsc_phastcons_track` | UCSC assembly and conservation track | `ce11`, `phastCons26way` | `sacCer3`, `phastCons7way` |
| `reference_data.ucsc_chrom_alias` | Map Ensembl/UCSC mitochondrial names that differ | `{"chrM": "MtDNA"}` | `{"chrM": "Mito"}` |
| `conservation_local.search_databases` | One entry per DIAMOND database; the TrEMBL entry's `name` must equal `trembl_name` | `swissprot`, `rhabditida_trembl` | same, with the yeast name |
| `batch_run.taxid` | NCBI taxonomy ID written into runs | `6239` | `559292` |
| `batch_run.max_workers`, `conservation_local.threads`, `domains_scan.cpus`, `deeptmhmm.batch_size` | Parallelism for the machine | `16` | `8` |
| `backends` | Local backends (below) | | |
| `proteome_db.path`, `run_dir`, `build_dir` | Final file, the per-protein run folders, and **local** disk where the database is assembled | `…/proteome_v1/proteome.sqlite3` | same pattern |
| `proteome_db.curated_sites`, `wormtagdb_strains`, `wormtagdb_notes` | Optional curated tag sites and strain records matched in by flanking sequence (see step 7) | | |

Backends for a fully local build (the template sets all but `topology`, which needs a separate
DeepTMHMM install, step 4):

```json
"backends": {
  "domains": "scan",
  "structure": "bulk",
  "genewise": "bulk",
  "conservation": "local",
  "topology": "deeptmhmm",
  "offtarget_region": "local"
}
```

The off-target screen also needs the organism's genome listed in `offtarget.config.json` under
`local.genomes` (taxid → database name and FASTA path); add an entry and re-run the `blastdb`
step. C. elegans (6239) and yeast (559292) are already there.

Check the config is the one in effect:

```bash
python -c "import sys; sys.path.insert(0,'scripts'); import providers; print(providers._config_path())"
```

> **Mixing organisms.** One config describes one organism. Use a separate config and a separate
> `out_dir` per organism, never reuse a run folder, and never read a database with another
> organism's genome (arms are regenerated from `reference_data.out_dir`).

---

## 2. Download the reference data

`scripts/reference_data.py` downloads each piece into `reference_data.out_dir`. Steps are
resumable (`curl -C -`) and skipped if the file exists; `--force` redoes a step.

```bash
# core: proteome, structures, genome and annotation
python scripts/reference_data.py --only uniprot,afdb,genome,gff3

# domains (Pfam HMMs, ~2 GB expanded) and conservation databases
python scripts/reference_data.py --only pfam,swissprot,trembl

# gene models and phastCons from UCSC, and the local off-target database
python scripts/reference_data.py --only ucsc,blastdb
```

| Step | What it fetches | Approx. size (C. elegans) | Organism-specific? |
|---|---|---|---|
| `uniprot` | Proteome FASTA and JSON | ~40 MB | yes |
| `afdb` | AlphaFold DB proteome tarball | ~2.8 GB | yes |
| `genome`, `gff3` | Ensembl genome FASTA and GFF3 | ~40 MB | yes |
| `ucsc` | UniProt-to-genome alignments and phastCons bigWig, plus an accession index | ~170 MB | yes |
| `blastdb` | BLAST+ database of the genome | ~25 MB | yes |
| `pfam` | Pfam-A HMMs and clans | ~2.7 GB | no, shared |
| `swissprot` | All reviewed UniProt, as a DIAMOND database | ~0.4 GB | no, shared |
| `trembl` | Unreviewed UniProt entries of one clade, as a DIAMOND database | ~1.5 GB for Rhabditida | clade only |
| `interpro` | `protein2ipr.dat.gz` streamed and filtered (13 GB download) | ~11 MB kept | only for `domains: bulk`; not needed with `scan` |

The shared files (`pfam`, `swissprot`) can be placed in one `out_dir` and reused: point each
organism's `out_dir` at its own folder and copy or symlink them, or reuse one folder if only the
clade file differs.

Check the downloads before going on:

```bash
ls -lh $(python -c "import sys; sys.path.insert(0,'scripts'); import reference_data as r; print(r._out_dir(r._load_config()))")
```

---

## 3. Build the local indexes

```bash
python scripts/local_store.py --build        # accession / CRC64 index of the proteome
python scripts/genome_regions.py --build     # GFF3 transcript and CDS index
python scripts/local_store.py --stats
python scripts/genome_regions.py --stats
```

`local_store.sqlite3` links each UniProt entry to its gene-model transcript through the entry's
`orfNames` cross-reference, and `genome_regions.sqlite3` holds the exons the gene models use. A
check that the cross-reference works for a new organism: `--stats` should show most proteins with
a transcript, and the "translated transcripts" count from `genome_regions.py` should equal its
transcript count.

---

## 4. Build the Pfam and topology caches (optional but recommended)

Both pre-scan the whole proteome so the per-protein run only reads cache files.

```bash
python scripts/build_pfam_cache.py           # hmmsearch over reference_data.protein_fasta
python scripts/build_topology_cache.py       # DeepTMHMM predictions (slow, GPU strongly preferred)
```

DeepTMHMM is licensed separately and installed outside the repository. Set
`deeptmhmm.install_dir` and `deeptmhmm.python` in the config to its checkout and its dedicated
conda environment. Without the cache, `topology_deeptmhmm.py` predicts per protein on demand, which needs DeepTMHMM
configured. Topology is not among the default `proteome_run.py` tasks.

---

## 5. Run the analyses for every protein

```bash
python scripts/proteome_run.py --out-dir <run_dir>              # canonical proteins
python scripts/proteome_run.py --out-dir <run_dir> --isoforms   # isoforms and annotation-only products
```

`<run_dir>` must equal `proteome_db.run_dir`. Each protein gets one folder; a line per protein is
appended to `_status.jsonl`, so a killed run resumes by skipping proteins already marked
`success`.

Default tasks are `domains,plddt,modifications,uniprot,scores,blast`:

| Task | What it computes | Local source |
|---|---|---|
| `domains` | Pfam domains | Pfam cache, else `pyhmmer` scan |
| `plddt` | pLDDT and SASA | AFDB PDBs (skipped where no AFDB model matches) |
| `modifications` | PTM motif sites | `tables/modification_sites.txt` |
| `uniprot` | Curated UniProt features | UniProt REST (canonical proteins only) |
| `scores` | Hydrophobicity windows | `tables/` |
| `blast` | Conservation | DIAMOND against both databases, then MAFFT + JSD |

Options: `--tasks domains,plddt` limits the tasks, `--limit N` and `--accessions-file` take a
subset, `--workers`, `--force`, `--reagents` (also design reagents; needs the gene-model backend),
and `--no-presearch-conservation`. Try `--limit 20` first and
read the failures in `_status.jsonl`:

```bash
tail -3 <run_dir>/_status.jsonl
```

A protein with no AFDB model logs `plddt: skipped: no AFDB model`; that is normal.

---

## 6. Estimate, then build the database

`build_proteome_db.py` always starts with an estimate. It builds a small random sample, measures
every table with SQLite's `dbstat`, and projects size and run time. A full build refuses to
start without a fresh estimate (under a week old), `--yes`, and enough free space.

```bash
python scripts/build_proteome_db.py --estimate
python scripts/build_proteome_db.py --estimate --reagent-sample 100   # also times reagent design
python scripts/build_proteome_db.py --stage ingest --yes
python scripts/build_proteome_db.py --stage scores --yes
python scripts/build_proteome_db.py --stage reagents --yes
```

`--stage all` and `--stage core` both run `ingest` then `scores`; `reagents` is always run
explicitly because it is the long stage. `--stage wormtagdb` refreshes only the curated strain
table (see step 7).

| Stage | Tables | Notes |
|---|---|---|
| `ingest` | `proteins`, `features`, `residues`, `conservation_summary`, `isoforms`, `run_status`, `meta` | A task's file is ingested only if its last status is `ok`. |
| `scores` | `score_configs`, `site_scores`, `tag_sites`, curated sites | Uses the app's own scorer, so values match the app. |
| `reagents` | `regions`, `insertion_sites`, `guides`, `site_guides`, `edits`, `genotyping_primers` | Runs the gene model, RS3 scoring and off-target screens in chunks. The ~1 kb arms are never stored. |

The reagents stage needs a config with `backends.genewise: "bulk"` and the local gene-model
sources, and is resumable. It writes transient per-protein TSVs (about 6 MB each) to
`scratch_dir`, deleted after each chunk of `reagent_chunk` proteins.

The database is assembled on local disk (`build_dir`) and then copied to `proteome_db.path`,
because SQLite locking is unreliable on network volumes. Useful options: `--stage all`,
`--update` with `--accessions` to refresh selected proteins, `--limit N` for a test build,
`--out` to write somewhere other than the config path.

---

## 7. Curated sites (optional)

`proteome_db.curated_sites` takes JSON lists of known tag sites (`gene`, `wormbase_id`, `site`,
`flank_n`, `flank_c`) that are matched to proteins by flanking sequence and verified against the
protein. `wormtagdb_strains` and `wormtagdb_notes` add a table of existing tagged strains. For an
organism without these, leave them out and the tables stay empty. The yeast build used curated sites and strain records; the build log reports "curated tag sites: N/M verified against sequence".

---

## 8. Verify the finished database

```bash
python scripts/proteome_db.py lookup <gene-name>
python -c "import sys; sys.path.insert(0,'scripts'); import proteome_db; print(proteome_db.describe())"
```

Then point `proteome_db.path` at the file and follow [LOCAL_DATABASE.md](LOCAL_DATABASE.md) to
run the apps against it.

Checks worth making on a new organism:
- Pick a protein with a known gene model and compare its residue numbering with the input
  protein. The build's gene-model gate (`cds_check`) fails clustered substitutions and long gaps
  rather than shifting residue numbers.
- Spot-check `proteins.sequence_differs` (a sequence that differs between the proteome index and
  the analysed run, usually from mixed UniProt releases).
- Confirm the reagent stage produced guides for a handful of proteins, and that arms
  regenerate: `python scripts/proteome_db.py export <name> <dir>`.

---

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| `local UCSC index not built` in the reagent log | Run `python scripts/reference_data.py --only ucsc`; without it, gene models fall back to the NCBI and Genewise paths. |
| No transcript for many proteins | The UniProt `orfNames` do not match the GFF3 transcript ids; compare a few in `local_store.sqlite3` against the GFF3. |
| Mitochondrial or plastid proteins missing coordinates | Add a `ucsc_chrom_alias` entry that maps the assembly's name for that chromosome. |
| `proteome_run` reports `conservation_local` errors | A `search_databases` name has no matching `{name}.dmnd` in `out_dir`. |
| Build refuses to start | No fresh `--estimate`, no `--yes`, or less free space than projected. |
| Off-target screens slow | The genome is not in `offtarget.config.json`'s `local.genomes`, or `blastdb` was not built; the run falls back to BLAT and EBI. |
| Database built from folders with stale files | Only a task's last status `ok` is ingested; re-run the failed tasks and rebuild the affected stage. |

For the batch review app and the single app's precomputed-results card, see
[LOCAL_DATABASE.md](LOCAL_DATABASE.md).
