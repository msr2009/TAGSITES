# Batch review app

Use this to review and choose tag sites for many proteins in one session, then download the
reagents for the ones accepted. The app reads a precomputed proteome database, so there are no
per-protein analysis runs.

It is local only (the hosted app has no database). It runs on **port 8001**, beside the main
app on 8000.

Before starting:

1. [INSTALL.md](INSTALL.md): the `tagsites` environment.
2. A database and its companion files: [LOCAL_DATABASE.md](LOCAL_DATABASE.md) (obtain one) or
   [BUILD_LOCAL_DATABASE.md](BUILD_LOCAL_DATABASE.md) (build one). Required pieces are the
   `proteome.sqlite3` file, the per-protein run folder, and the organism's genome files.
3. A config pointing at them ([CONFIG.md](CONFIG.md)).

## 1. Configure

Start from the template and fill in the paths:

```bash
cp batch.config.example.json batch.config.local.json
# edit: reference_data.out_dir, proteome_db.path, proteome_db.run_dir, batch_app.batch_dir
export TAGSITES_BATCH_CONFIG=batch.config.local.json
```

The `batch_app` block:

| Key | Meaning | Default |
|---|---|---|
| `batch_dir` | Local-disk folder where each batch saves its files (`~` is expanded) | `~/tagsites_local/batches` |
| `workers` | Parallel database exports during Prepare | 6 |
| `ebi_max_concurrent` | Proteins analysed on EBI at once (only for proteins outside the database) | 3 |
| `cache_entries` | Parsed runs kept in memory for Review, about 3 per protein | 500 in code, 1500 in the template |
| `port` | App port | 8001 |

`app_batch.py` falls back to `batch.config.worm_local.json` when `TAGSITES_BATCH_CONFIG` is
unset; the export above overrides that.

## 2. Check, then start

```bash
python scripts/proteome_db.py lookup trxr-1      # a table of products means the database is found
python app_batch.py
```

Open **http://127.0.0.1:8001**. The tabs are **(1) Setup Batch**, **(2) Prepare** and
**(3) Review Sites**.

## 3. Setup Batch

Work down the numbered panels.

0. **Proteome database**: shows the file in use (path, size, protein count). To use another file,
   enter its path and click **Use database**; the file is checked first. A batch records its
   database, and reopening the batch switches back.
1. **Open a saved batch**: resume earlier work (see [Save and reopen](#6-save-reopen-share)).
2. **Proteins**: name the batch, paste or upload a list, then **Look up proteins**. Accepted
   input: gene names, accessions or WormBase ids one per line, an optional insert column, or
   FASTA. A gene with several products gets a dropdown. The default product is the one whose id
   equals the search term, else the UniProt entry, then the one with reagents designed, then with
   a structure, then the longest. Duplicates are flagged.
3. **Reagents**: repair template (Homology arms, ssODN, PCR amplicon, SapTrap plasmid), arm
   length, and the insert (a named tag or **Custom insert DNA**). These apply to every accepted
   site.
4. **Suggested sites**: how many sites to suggest per protein (default 2) and the masks.
   *Position*: no N-terminus, no C-terminus, guide cut within N bp. *Annotation*: constitutive
   sites only, signal / transit peptide and transmembrane (on by default), Pfam domains, all
   UniProt sites. These are defaults for switches in Review's Rescoring panel and only mask
   suggestions; they never remove a site or a guide.
5. **Proteins not in the database**: optionally run the web-services pipeline for them (needs an
   e-mail address). Each such protein takes minutes.
6. **Save batch & continue**: saves the batch and starts Prepare.

## 4. Prepare

Prepare starts automatically and shows one status row per protein.

- Each protein is exported to `{batch_dir}/{batch}/{key}/` as an ordinary run folder, without
  homology arms (those are regenerated only for accepted sites).
- Proteins outside the database run on EBI when panel 5 was enabled.
- The in-memory cache fills so Review switches proteins in 25 to 500 ms.
- Alignment PDFs render afterwards in a background process, so Review can start before they
  finish.

About 10 s for 10 proteins on the reference machine (mostly one-time imports and a 2 s genome
load). Prepare is resumable: finished folders are skipped.

## 5. Review Sites

The single-protein Results view, plus a navigation bar:

- **◀ Prev / Next ▶** and a colour-coded protein list.
- Per-site **guide picker**, defaulting to the guide nearest the insertion site. RS3 scores are
  shown as badges and do not affect selection.
- **Accept → Next** accepts the protein with its chosen sites and moves on.
- **Skip** (next to Accept) marks the protein skipped (red italic). Skipped proteins are listed
  in `{batch}.skipped.txt` beside the batch file and never reach Output. **Unskip** restores the
  earlier status.
- **Save batch** writes the batch file.

Statuses run pending → ready → in_progress (opened) → accepted. Editing an accepted protein
returns it to in_progress. Suggested sites are recomputed whenever a non-accepted protein is
opened; only accepted proteins keep their sites.

## 6. Output, save, reopen, share

- **⬇ Output reagents (.zip)**: reagent files for accepted proteins plus combined TSVs.
  Reagent names are `{gene}_{uniprot}_{aa}{res}_{tag}`, or the product id when there is no gene.
- **⬇ Batch (.zip)**: the whole batch as one file (`BATCH_MANIFEST.json` plus one ordinary run
  bundle per protein). Reopen it in panel 1 of Setup on any machine with the same database. The
  main app can open any single protein's bundle; it recomputes reagents (about 9 s) for database
  proteins because the arm-less TSV is not included.

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| Panel 0 shows no database | `proteome_db.path` is wrong, or a different config is active (`echo $TAGSITES_BATCH_CONFIG`) |
| Prepare is slow | `batch_dir` is on a network volume; use local disk |
| Structure or alignment "not in the batch folder" | `proteome_db.run_dir` is wrong or not mounted |
| Arms fail to generate on Output | `reference_data.out_dir` lacks the genome, GFF3 or `genome_regions.sqlite3`, or is another organism |
| Config edits have no effect | The config is read once per process; restart |
| Port 8001 in use | Another batch app is running, or change `batch_app.port` |
