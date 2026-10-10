# Precomputed-results wrapper for the app (browse the proteome database)

## Context

`proteome.sqlite3` now holds every per-residue track, feature, score, isoform map and tag
site for all 28,626 proteins, but the app cannot show any of it: its Results and Reagents
tabs read files named in a run JSON. Today every protein needs a fresh analysis run (minutes,
remote services) even when it is already precomputed.

Outcome: after a user searches for a protein or pastes a sequence in the Setup tab, the app
lists the precomputed proteins that match it and loads one straight into the Results tab,
with no analysis run, no network and no change to the Results/Reagents code. Where the
database is absent (shinyapps.io) the new panel does not render and the app is unchanged.

Decisions made with you: matches are **identical + other isoforms of the same gene +
contained** (no homology search); the panel lives **inside Setup, under the sequence
input**; alignments and AlphaFold PDBs are **folded into the database** so one file is
enough to host or send.

## Integration seam (found in the code)

- All four tabs talk through one reactive value, `shared_values` (the active run JSON path)
  in `server.py`. `results_server.auto_load_results` reloads whenever it changes, and
  `reagents_server` reads the same run JSON. So the wrapper only has to **write a run folder
  and set that path**; the Results and Reagents tabs then work unchanged.
- Tab switching already exists: `ui.update_navs("main_tabs", selected="results",
  session=session.root_scope())` (`progress_server.py:580`).
- `load_data_from_json` positions continuous tracks by **row order** (row i -> pos i+1) and
  reads file names from the task registry (`.jsd`, `.txt`, `.tsv` with companions `.aln`,
  `.isoforms.json`, `.sasa.txt`, `.hydro.txt`, `.patches.txt`), so exported files must be
  contiguous per residue with `-1000` for gaps.
- The alignment pane embeds `{aln}.pdf`, which the batch run never rendered, so the export
  renders it (`build_heatmap_reportlab.plot_alignment_reportlab`). 3D view needs
  `{run}.AF.pdb`.
- Setup already has the three sequence sources (UniProt hit > uploaded file > pasted
  sequence), resolved with that precedence at save time (`setup_server.py` ~1304-1350).
- Measured: alignments total ~2.6 GB raw (90 KB each); AlphaFold PDBs 1.1 GB gzip.

## Design

### 1. Database additions (`scripts/build_proteome_db.py`)
- Index `proteins(crc64)` for exact-match lookup (applied by any `--update`; the immutable
  read-only file cannot create indexes at runtime).
- New `assets(pid, kind, data)` table (WITHOUT ROWID): `aln` = zlib-compressed alignment,
  `pdb` = the AlphaFold `.pdb.gz` bytes as shipped. New stage `--stage assets` reads each
  protein's `{id}_conservation.aln` (only when blast status is ok) and, for
  `has_structure=1`, `data/reference/afdb_pdbs/AF-{accession}-F1-model_v6.pdb.gz`.
- Size test first, as before: `--estimate` gains an assets sample (compress ~300 alignments,
  stat the PDB gz files) and projects the added size; the full stage needs `--yes` and
  passes the free-space gate. Expected: ~0.5 GB (alignments) + 1.1 GB (PDBs), so
  0.84 GB -> ~2.5 GB.

### 2. Matching (`scripts/proteome_db.py`)
- `normalize_sequence(text)`: strip FASTA header/whitespace/`*`, uppercase.
- `sequence_matches(seq, conn)` returns one row per proteome protein with a `relation`:
  - `identical`: CRC64 index lookup (includes duplicate-sequence twins);
  - `isoform of matched gene`: other proteins sharing the matched protein's `wb_gene` or
    base accession, shown with their length difference;
  - `contains query`: `instr(sequence, :q) > 0` for queries of >= 20 residues;
  - `contained in query`: `instr(:q, sequence) > 0` for proteins of >= 30 residues (a
    tagged or fusion construct).
  Each row also carries `gene_name`, `kind`, `seq_length`, `has_structure`, whether
  reagents exist, and where the query sits in the protein. All of it is SQL over 28.6k rows
  (milliseconds).
- `search_proteins(term, limit)`: case-insensitive prefix match on gene name, accession,
  id and WormBase names, for the panel's own search box (`find_protein` stays exact).

### 3. Run-folder export (`scripts/proteome_export.py`, new shallow file)
`export_run(pid, out_dir, conn=None) -> path to {run}.run.json`, writing the files the
loaders expect, with task names the app already colours by substring
(`PRECOMPUTED_blast`, `PRECOMPUTED_plddt`, `PRECOMPUTED_scores`, `DOMAINS_domains`,
`MODS_modifications`, `UNIPROT_uniprot`, `TOPOLOGY_topology`):
- `{run}.fa` from `proteins.sequence`.
- `.jsd`: `#` header lines, then one row per residue `pos score aa`, unscaled, `-1000`
  where conservation is NULL. `.scores.tsv`: contiguous rows 1..N of the K-D track.
  `.txt`/`.sasa.txt`/`.hydro.txt`/`.patches.txt` only when `has_structure`.
- Range files from `features`: domains, uniprot, topology, modifications. **Modification
  stops are written as stop+1**, the app's own convention, so the live Results/scores match
  both the app and `site_scores` (the DB stores the corrected stop; see follow-ups).
- `.isoforms.json` rebuilt from `isoforms` + `isoform_spans`; `.aln` and `{run}.AF.pdb` from
  `assets`; the alignment PDF rendered lazily and cached in the export folder.
- `global` block: run name, working dir, input file, pdb, taxid 6239,
  `selected_sites` = curated sites for that protein if any.
- Reagent files (phase 2 below) only when the reagent tables hold rows for the protein.

### 4. UI (`modules/precomputed_ui.py` + `modules/precomputed_server.py`, nested in Setup)
- A new panel "Precomputed results" under the protein-sequence input (`setup_ui.py`
  panel 2), rendered only if `proteome_db.available()`.
- `setup_server` exposes a small reactive `current_protein()` (name, sequence) using the
  existing UniProt > file > paste precedence, refactored out of the save-time code, and
  passes it with `shared_json` and the session temp-dir getter (`server.py:_get_session_dir`)
  to `precomputed_server`.
- The panel shows a status line ("3 precomputed proteins match"), a table of matches with
  relation badges (identical / isoform / contains / contained), and a "View results" button
  per row (the `use_hit_{id}` button pattern already used for UniProt hits). It has its own
  name/accession search box for direct browsing.
- Clicking runs `export_run` through `setup_server._off_loop` (blocking file writes and PDF
  render must not run on the shared event loop, per CLAUDE.md), sets `shared_json` to the
  new run JSON (a fresh folder per load, so the path always changes), then switches to the
  Results tab. No match: "No precomputed results (the database covers C. elegans)".
- `server.py` passes the session-dir getter into `setup_server`; temp folders are already
  removed by `_cleanup_session`.

### 5. Reagents (phase 2, once the reagent stage has run)
The reagent tables are empty today, so Results works first and the Reagents tab shows its
existing "not found" state. When filled, the export also writes `{run}.reagents.tsv` with
arms regenerated by one `regenerate_all_arms(pid)` pass (region and frame lookup computed
once, not per row), `.genotyping.tsv`, and the genewise/genomic files from `regions` so the
GenBank export works. The detailed off-target sidecar is not stored, which the Reagents tab
already tolerates.

### 6. Configuration and hosting
- `batch.config.json` `proteome_db.path` stays the single setting; an optional
  `TAGSITES_PROTEOME_DB` environment variable overrides it for hosts.
- `proteome_db.available()` is false when the file is missing, so shinyapps.io and any
  install without the database behave exactly as today. Reads use the existing
  `mode=ro&immutable=1` connection, one per call, so concurrent sessions are safe.

## Critical files
- New: `scripts/proteome_export.py`, `modules/precomputed_ui.py`,
  `modules/precomputed_server.py`
- Modified: `scripts/proteome_db.py` (matching, search, `available`), `scripts/build_proteome_db.py`
  (crc64 index, `assets` table and stage, estimator), `modules/setup_server.py` and
  `modules/setup_ui.py` (panel, `current_protein`), `server.py` (pass session dir),
  `batch.config.json`, `CLAUDE.md`
- Reused: `utils/results.load_data_from_json` and `scripts/task_registry.companion_path`
  (file naming), `modules/setup_server._off_loop`, the `use_hit_*` button pattern,
  `proteome_db.regenerate_arms`, `build_heatmap_reportlab.plot_alignment_reportlab`

## Verification
1. **Estimate before building assets:** `--estimate` reports the projected assets size;
   build only after you confirm.
2. **Export parity:** for ~30 random proteins (with and without structure, an isoform, a
   WormBase-only id, a selenoprotein, one of the 28 re-run proteins), run
   `load_data_from_json` on the exported run and compare `aa_df`/`range_df` with the same
   call on the original batch files (values equal to the stored precision, features
   identical), and `score_tag_sites` with the database's `site_scores` (expect ~100%
   identical; any difference only where a value sits within 1e-5 of a threshold).
3. **Isoform round trip:** rebuilt `.isoforms.json` equals the original for the sample.
4. **Matching:** a known protein pasted as FASTA returns `identical`; a UniProt hit for a
   selenoprotein matches; a 60-aa fragment returns `contains query`; a query with a
   Strep-tag appended returns `contained in query`; an isoform returns its siblings; a human
   protein returns "no match" without error.
5. **End to end in the running app** (`python app_modular.py`): search `trxr-1`, pick a
   match, land on Results with tracks, alignment pane (PDF renders) and 3D structure; the
   Reagents tab shows the not-precomputed state. Time the load (target a few seconds) and
   the alignment PDF render.
6. **Responsiveness:** with two browser sessions, a load in one must not freeze the other
   (event loop not blocked).
7. **Absent database:** with `proteome_db.path` pointing nowhere, the Setup tab renders
   without the panel and the existing tests pass (625 + 6 skipped).
8. **Tests** (only if approved, per your preference): export-vs-loader parity on synthetic
   folders, the four match relations, and `available()`.

## Follow-ups, not in this plan
- `regex_sites.py` writes modification stops one residue too long; fixing it needs a rebuild
  of the affected outputs and would remove the stop+1 special case in the export.
- Option B (running new analyses offline) and non-C. elegans species are out of scope.
- Homology (similar-sequence) matches were deliberately left out; a DIAMOND search against
  the 28.6k proteome could be added as a second tier later.
