# Planning the proteome database rebuild

Both `proteome.sqlite3` files were built on 2026-11-03 at commit `61d0824` with stages
`ingest,scores,reagents`. Several changes have landed since that the stored tables do not
reflect. This note records what is stale, what each stage would fix, and the order to do it in.
Nothing here has been run yet.

## What is stale

| Stored data | Why it is out of date | Stage that fixes it |
|---|---|---|
| `tag_sites` | the picker now offers one site per score-equivalent region, not a fixed three, and `f52a8cc` bridges dips | dropped, not rebuilt — see below |
| `sequence_aliases` | table does not exist in either file | `ingest` |
| `proteins.structure_source` | only `afdb` / NULL; never `afdb:ACC` | `ingest` |
| `features` (`secondary_structure`) | helix/strand track added after the build | `ingest`, **after a pre-pass** |
| `features` (`modification`) | PTS1, H/KDEL, KKXX, FFAT added and CAAX widened | re-run the modifications task |

Note the modification re-run is **corrective as well as additive**. A stored motif call is not
re-derived at ingest: `features` rows are copied from whatever `{id}_mods.txt` holds, so a hit
produced by a regex that no longer exists survives every rebuild until the task itself is re-run.
DBL-1 is the worked example — `data/DBL-1/DBL-1.MODS_modifications.txt` carries
`CAAX prenylation 362-366` on the C-terminal `...VESCGCR`, which matches neither the current
pattern nor its predecessor. It was produced when CAAX was the unrestricted `C...$`, tightened in
`035a2cd`. Because `prenylation` is one of `in_highconf_modification`'s
`include_descriptions`, that stale line hard-masks the residue, and the terminus ladder then
reports DBL-1's C-terminus as `masked` — on a secreted TGF-beta ligand whose C362/C364 are the
curated cystine-knot disulfide partners.

The app recomputes suggested sites live from the exported track files, so a stale `tag_sites`
affects stored queries only, not the UI. The other four are genuinely missing data.

## Current state (audited 2026-11-05)

```
afdb_pdbs tarball: 19,701 models

worm: 28,626 proteins
  has_structure = 1 in the database : 25,569
  .AF.pdb present in the run folder : 25,569      <- exact match, no drift
  model in tarball under own acc    : 21,266      <- isoform rows share a canonical accession
  model only under an alias acc     :    113
  NEW structures the aliases recover:    103

yeast: 6,742 proteins
  has_structure = 1 in the database :  6,055
  .AF.pdb present in the run folder :  6,055      <- exact match
  model in tarball under own acc    :      0
```

Every protein the database claims has a structure has its file on disk, in both organisms. The
4,303 worm models not covered by the tarball were fetched at run time into the run folders.

Two findings that shape the plan:

- **The alias work recovers 103 worm proteins** that currently have no structure. Modest (0.4%),
  and it requires the structure task to actually run — `ingest` can only record `afdb:ACC` for a
  protein whose run folder already holds the `{id}.uniprot_accession.txt` sidecar. It cannot
  create a structure that was never fetched.
- **Yeast has no local AFDB coverage at all.** `data/reference/afdb_pdbs` and `local_store` are
  the worm reference, so the yeast alias stage will insert nothing, and the 6,055 `.AF.pdb` files
  in the yeast run folders are the only local copy of those structures. They are intact but not
  re-derivable. A yeast alias index would need its own reference data.

## The reagents stage does not need re-running

This is the expensive stage (genome, local backends, hours). It can be preserved:
`insert_proteins` reuses each protein's existing `pid`, keyed on `id`, and the reagent tables are
in neither `STAGE_TABLES["ingest"]` nor `["scores"]`. Running `ingest,scores` with `--update`
against a **copy** of the existing file keeps `guides`, `site_guides`, `regions`,
`insertion_sites`, `edits` and `genotyping_primers` intact.

## Stop storing `tag_sites`

**Decided 2026-11-08: the rebuild should drop the precomputed suggested sites.** They are a
derived value with a short shelf life and no reader in the app.

- **Nothing in the app or the batch app queries the table.** The only writer is
  `build_proteome_db._scores_stage`, and the only readers are offline analysis scripts
  (`scripts/export_tag_site_table.py`). Every path a user takes — Results, the batch Review tab,
  `proteome_export` — recomputes suggestions from the scores frame on load.
- **Recomputing is trivial.** `find_regions` + `pick_in_region` is one pass over a few hundred
  residues, microseconds per protein against the tens of milliseconds a protein switch already
  costs. There is nothing to amortise.
- **Storing them manufactures staleness.** The rows are frozen against the masks, the weights and
  the region rules in force at build time. Every later change to `scores.config.json` silently
  desynchronises them from what the app shows, with no mechanism that notices — commit `f52a8cc`
  (dip bridging via `score_tolerance` / `region_gap`) already invalidated every stored row, and
  the suggestion machinery is the part of the config that changes most often.
- **`site_scores` is the useful artefact and stays.** Per-residue score, mask flag and criterion
  bits are what the offline analyses actually need, and the region rules can be re-applied to them
  after the fact — that is exactly how the dip-bridging measurements were run.

What this means for the rebuild: keep `site_scores` in the `scores` stage, drop the `tag_sites`
writes (both the `source = 'suggested'` rows at `build_proteome_db.py:708` and the curated rows at
`:763`, or keep the curated ones alone if an analysis depends on them), and point
`export_tag_site_table.py` at `site_scores` instead. If the table is kept regardless, it needs a
`score_config_hash` column so a stale row is at least detectable.

## Decide first: what earns a hard mask

**Settle this before the rebuild if option 3 is in play.** With `tag_sites` dropped (above), the
only stored thing a mask policy can change is `site_scores`. Option 2 leaves the per-residue score
untouched and so is no longer blocking — it can be adopted any time by editing the `terminus`
block. Option 3 changes `site_scores` for ~1,039 proteins and therefore means re-running the
`scores` stage, which is why it is cheap to decide now and expensive to revisit.

Three kinds of evidence currently collapse into the single tier `masked`, and they are not of
equal strength:

| path | evidence | worm termini | defensible as tier 1? |
|---|---|---|---|
| `topology` | `in_membrane_or_cleaved`: DeepTMHMM / Phobius / UniProt signal, TM, propeptide, transit | 819 alone, 5,230 with curated | yes — independent predictors, usually several agreeing |
| `curated` | a UniProt / UniProt_site feature in the terminal window matching `curated_blocking` | 1,631 alone | yes — a curator asserted the function |
| `mod-mask` | `in_highconf_modification`: a **regex** from `tables/modification_sites.txt` whose description contains myristoylation, palmitoylation or prenylation | **1,039 alone** | **no** — this is the open question |

The third path is the problem. A regex hit with no curated annotation behind it is exactly what
the terminus ladder defines as `putative`, yet it arrives at `masked` through the scores frame,
bypassing the tier split entirely. The measured precision of those motifs against UniProt's own
terminal annotations says the distinction is real: H/KDEL 100%, CAAX 38.6%, PTS1 20.0%, RAB
prenylation 9.1%. Only CAAX, RAB and the palmitoylation motifs feed the mask, and two of those
three are the least precise in the table.

### The worked example

Three Rab GTPases, all with a genuinely prenylated C-terminus, land in three different tiers:

| gene | acc | C-terminus verdict | path |
|---|---|---|---|
| rab-11.1 | O01803 | `masked` | curated — `Lipidation: S-geranylgeranyl cysteine` |
| rab-7 | Q23146 | `masked` | **regex alone** — `RAB prenylation` |
| rab-5 | P91857 | `putative` | regex, `di-cysteine palmitoylation site` (not in `include_descriptions`) |

Same biology, three answers, and which one a protein gets depends on how well UniProt happens to
have curated it and on which motif name the regex carries. Two more that make the cost concrete:
**sel-12** (P52166), a nine-pass intramembrane protease hard-masked by a CAAX match, and
**unc-25** (Q8I4E0), cytosolic GAD hard-masked by `RAB prenylation` at 9.1% precision.

### Options

1. **Status quo.** Any mask is tier 1. Simple; keeps 1,039 regex-only worm termini behind a
   confirm click they probably should not need, two thirds of them wrongly if CAAX's 38.6% is
   representative.
2. **`curated_mask_keys` in the `terminus` block.** A list naming which mask criteria count as
   curated evidence — `["in_membrane_or_cleaved"]` — with every other mask demoting the terminus
   to `putative`. One config key, no new matching logic, and it leaves the per-residue score
   untouched: the residue stays masked for scoring and only its *terminus tier* changes. This is
   the recommended option.
3. **Split the criterion.** Separate `in_highconf_modification` into a curated-backed and a
   regex-only criterion. More faithful, but it changes the per-residue score for ~1,039 proteins
   and so changes `site_scores` as well as `tag_sites`.

Option 2 changes what the `scores` stage stores for the regex-only population and nothing else,
which is why it wants deciding before the stage runs rather than after.

### Proteins to check in the app before deciding

| gene | acc | N | C | case |
|---|---|---|---|---|
| gly-8 | O45293 | masked (curated) | masked (curated) | both ends curated |
| dpy-7 | P34688 | masked (curated) | ok | N only |
| daf-22 | G5EDP2 | ok | masked (curated) | C only, PTS1 motif annotation |
| lgg-1 | Q09490 | ok | masked (curated) | the only PubMed-cited reason in the set — tests `curated_requires_experimental` |
| sel-12 | P52166 | ok | masked (**regex only**) | tier-1 leak |
| unc-25 | Q8I4E0 | ok | masked (**regex only**) | tier-1 leak |
| unc-54 | P02566 | putative | ok | putative at N |
| rab-5 | P91857 | ok | putative | putative at C |
| let-23 | A0A0K3AX35 | masked (curated) | putative | one tier at each end |
| dbl-1 | G5EEL5 | masked (curated) | ok | control: the database says C is fine; the old `data/DBL-1/` folder says masked |

Database-wide for worm: 14,810 flagged termini — N 7,123 masked / 3,968 putative, C 1,672 masked
/ 2,047 putative.

## Order of work

1. **Structure re-run**, worm only, for the ~3,057 proteins with `has_structure = 0`. Mostly
   local file checks plus some AFDB fetches; writes the models and their accession sidecars, so
   the 103 alias-recoverable proteins gain a structure and `ingest` can label them `afdb:ACC`.
2. **Secondary-structure pre-pass.** No `{id}.ss.txt` exists anywhere, but 25,569 worm and 6,055
   yeast run folders still hold their `.AF.pdb`, so `dssp_numpy.write_ss_file` regenerates all of
   them locally — no downloads, no re-running pLDDT or SASA. `build_proteome_db.py:522` already
   ingests the file once it is there.
3. **Modifications re-run.** Pure regex over the sequence, no network, cheap.

   ```bash
   python scripts/proteome_run.py --out-dir <run_dir> --tasks modifications --force
   ```

   `--force` is required: every protein is already marked complete in `_status.jsonl`, so the
   resume logic would otherwise skip all of them. A task-subset re-run is safe for the other
   analyses — `read_last_status` takes each task's status from the last line that ran *that*
   task, so a modifications-only pass does not make the untouched tasks look unrun. The task
   overwrites `{id}_mods.txt` in place, and `features` is in `CHILD_TABLES["ingest"]`, so the
   later `--update` deletes each protein's feature rows before re-inserting; stale calls
   disappear rather than accumulating.

   Scale, measured over the first 4,000 worm run folders by re-running the current table against
   each stored sequence:

   | motif | stored | with the current table | delta |
   |---|---|---|---|
   | CAAX prenylation | 7 | 42 | +35 |
   | KKXX ER retention | 0 | 40 | +40 |
   | PTS1 peroxisome targeting | 0 | 20 | +20 |
   | FFAT / phospho-FFAT | 0 | 2 | +2 |
   | everything else | unchanged | unchanged | 0 |

   **Nothing is lost in the proteome run.** It was built with the tightened
   `C[AILMV][AILMV].$`, so the stored CAAX calls are a subset of what the current pattern finds;
   the promiscuous `C...$` artifact is confined to older single-protein run folders under `data/`
   (see below). Extrapolated to 28,626 worm proteins the re-run adds roughly 250 CAAX, 286 KKXX
   and 143 PTS1 calls.

   One consequence to decide **before** the `scores` stage: of those, only CAAX feeds
   `in_highconf_modification`, so ~250 worm proteins gain a hard-masked C-terminal residue from a
   motif measured at 38.6% precision. KKXX, PTS1, H/KDEL and FFAT are not in that criterion's
   `include_descriptions`, so they cost a point via `not_in_modification` and place the terminus
   in the `putative` tier, which is the intended weight for a regex-only call. Whether a
   regex-only prenylation hit should likewise be demoted from tier 1 to tier 2 — a
   `curated_mask_keys` knob in the `terminus` block — is open; see **Decide first: what earns a
   hard mask** above, which must be settled before step 4.
4. **Settle the hard-mask question** (previous section). Blocking only for option 3, which
   changes `site_scores`; option 2 is a config change that can land before or after.
5. **`--estimate`, then `ingest,scores` with `--update`** on a local copy of each file, then copy
   back to the NAS. SQLite locking is unreliable on the network volume, so the file is always
   assembled on local disk.

Disk at the time of writing: 261 GB free; worm is 5.9 GB and yeast 1.2 GB. The secondary-structure
features table is the only material addition.

## Old single-protein run folders under `data/`

These are not part of the proteome run and no rebuild stage touches them, but the app loads them
directly and they carry the oldest annotations — DBL-1's `C...$` prenylation call is one. They
are test and demo folders rather than data of record.

**Decided 2026-11-08: delete them all, then reconstitute a few.** Re-running the modifications
task over them would fix only the motif rows while leaving every other track as old as the folder
is, and a folder that is partly current is worse than one that is plainly fresh: it looks
trustworthy and is not. The cost of regenerating one is now small, because the database can
produce a complete run folder directly.

1. **Take the list first.** Record the folder names and what each was kept for before anything is
   removed; some are the only trace of a bug that motivated a fix (DBL-1 and the `C...$` call are
   the example this document already leans on). The names are the record, not the files. The
   inventory as of 2026-11-08, 53 entries, of which four are **not** single-protein folders and
   must survive: `reference/` (genome, AFDB tarball, UCSC and phastCons copies), `cache/`,
   `runs/`, `SCORE_TESTING/`, plus the batch folders `20260727_batch/`, `20260728_batch/` and
   `internal_batch/`:

   ```
   A0A7I4YC85_20260804  abl-1  apm2  CHA-1  CHA1_20261011  CHA1_20261012  cha1
   clr-1_20260803  col-103  DBL-1  eevs_20260804  frh-1  frh1  his3  hlh-30  hras
   ifb-2_20260803  kin-20_20260803  knl-1_20260804  let-23  leu2  lgg-1  LHDH_worm
   madd2  mut16  mut16_new  Ppa-mbd-2  ret1  snap29  snt-1_20260803  src-1
   src-1_20261011  src1  TMEM240_20260804  ubiquitin  UNC-17  unc-33  unc-33_20260803
   unc33  unc-40  unc-47  unc17  UNC25  UNC44_ens  vab-1_20261011  white_20260804
   ```

   Note the duplicates with different spellings (`src-1` / `src1` / `src-1_20261011`,
   `CHA-1` / `cha1` / `CHA1_20261011` / `CHA1_20261012`, `unc-33` / `unc33` / `unc-33_20260803`,
   `frh-1` / `frh1`, `mut16` / `mut16_new`, `UNC-17` / `unc17`), which is itself the argument for
   clearing them: no one can now say which copy is the good one.
2. **Delete every single-protein folder under `data/`.** Not the reference data —
   `data/reference/` holds the genome, the AFDB tarball and the UCSC/phastCons copies and must
   stay. Nothing in the app or the test path depends on a specific folder existing; the app
   creates its own run folders, and the batch app works out of `batch_app.batch_dir`.
3. **Reconstitute a handful from the rebuilt database**, via `proteome_db.py export <name> <dir>`,
   once the rebuild is verified. Each comes back as an ordinary run folder with every track
   current and regenerated arms, so the app loads it unchanged.

Which proteins to keep is the only open part. The set should cover the cases that get tested by
hand rather than being a sample of the proteome: a clean protein with both ends OK, a
curated-masked terminus, a regex-only masked terminus (the tier-1 leak case), a protein whose
gene model only warns, a multi-TM protein, and one with no structure. The hard-mask check table
above already names candidates for most of these.

Proteins loaded through the precomputed-results card are exported from the database and are
unaffected by any of this.

## Verification after the rebuild

- `tag_sites` is gone (or holds only curated rows), and `site_scores` for a sample of proteins
  reproduces `utils.scoring.score_tag_sites` on the same exported run folder.
- `sequence_aliases` is non-empty for worm and `structure_source` contains `afdb:ACC` values.
- `features` contains a `secondary_structure` source, and the new modification motifs appear.
- Re-running the current regex table against each stored sequence reproduces every
  `modification` row exactly — that is, no stored call the table can no longer generate. This is
  the check that would have caught DBL-1.
- `guides` / `site_guides` row counts are unchanged from before the rebuild.
- `has_structure` rises by roughly 103 for worm.
- The ten proteins in the hard-mask check table give the tiers that section predicts, under
  whichever option was chosen — in particular sel-12 and unc-25, which move from `masked` to
  `putative` under option 2 and stay `masked` under option 1.
