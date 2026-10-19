# Plan: gnomAD constraint range track

Build a single `GnomAD_RMC` range track: gnomAD v2 regional missense constraint (RMC) from
every Alliance human ortholog, merged into one consensus set of worm-coordinate blocks.
Background and prototype results are in `docs/GNOMAD_CONSTRAINT_SPEC.md`. The
per-residue MTR track and the evaluation code stay in the prototype and are out of scope here.

## Scope

In:
- New `gnomad` task (`result_type: range`), runnable from the app and the CLI.
- A consensus merge of the per-ortholog RMC regions into one track.
- Plot support: one row, three class colours, details in the hover text.
- Shipping a pre-filtered worm-human orthology table.

Out (later):
- Per-residue MTR as a continuous track.
- A structure colour scheme from the track.
- Species other than C. elegans.
- A proteome-scale batch path (the API is too slow; needs gnomAD bulk downloads).

## Decisions already made (from the prototype)

- **Orthologs**: all Alliance human orthologs, not reciprocal-best only. The reciprocal-best
  rule dropped KRAS and MYH7, which carry most of the ClinVar signal.
- **Consensus**: median o/e across the orthologs covering a segment.
- **Classes**: constrained o/e <= 0.4, intermediate, tolerant >= 0.8 (absolute; RMC is
  already modelled against a mutation-rate expectation, unlike raw MTR).
- **Never filters**: the track annotates sites and guides; it never removes one (same rule as
  RS3 and conservation).

## Steps

### 1. Orthology table — `tables/alliance_worm_human_orthologs.tsv`

- Filter `ORTHOLOGY-ALLIANCE_TSV_COMBINED.tsv.gz` to worm -> human rows (15,460 rows,
  7,432 worm genes; about 1 MB). Keep the columns the script reads and a header line recording the
  Alliance version (9.1.0) and download date.
- Add a `--build_orthology <combined.tsv.gz>` mode to the script that writes this table, so
  it can be regenerated when Alliance releases.
- `load_orthologs()` reads the small table by default and still accepts the combined file.

### 2. Worm protein -> WBGene from the input FASTA

The task gets a sequence, not a gene symbol. Resolution order:
1. The `accession` param, if set (same convention as the `uniprot` task).
2. The UniProt accession in the FASTA header.
3. `uniprot_api.checksum_lookup(seq, taxid=6239)`.
Then `uniprot_api.fetch_entry()` and the WormBase `GeneId` xref (`_wbgene_from_entry`).
In batch mode, use `local_store`'s `wormbase_gene` column when the store exists.

If no WBGene or no ortholog is found, write an empty output file and log the reason. Missing
data is not an error (par-2 and trp-4 behave this way now).

### 3. Script — rename the prototype to `scripts/gnomad_rmc.py`

Keep `scripts/gnomad_constraint.py` as the research prototype. The new script reuses its
fetch, numbering-check and alignment functions by importing them, not copying them.
- `main(fasta_in, output, accession="", orthology=..., cache_dir="", report=None)`, with a
  CLI that mirrors `uniprot_features.py` (`--input_file/--fasta`, `--output`, `--accession`).
- Per ortholog: gnomAD GRCh37 RMC query, then the GRCh37 protein, mapped to MANE when they
  differ (KRAS), then the MANE protein, then projected through the worm alignment. Only the
  GRCh37 query and the GRCh38 `mane_select_transcript` field are needed; the multi-MB
  variant list is not, which cuts payload and runtime sharply.
- Numbering checks stay hard gates: an ortholog with any RMC boundary mismatch is skipped
  and listed in the sidecar.
- Network calls go through `http_retry.request_with_retries` with a wall-clock deadline
  (`deadline_from`). The task runs off the event loop like the other tasks; nothing new
  is called from a Shiny effect.
- Cache raw responses under the run folder (`gnomad_cache/`) so a rerun costs nothing.

### 4. Consensus merge — `merge_rmc_regions(per_ortholog_regions, worm_len, ...)`

1. Each ortholog contributes worm-coordinate regions `(start, end, oe)` from
   `project_regions()`.
2. Breakpoints = union of every region's start and end + 1; this splits 1..worm_len into
   elementary segments.
3. For each segment: the orthologs covering it, the median o/e, and each ortholog's class.
   Segments covered by none are dropped (the plot shows nothing there, read as "no data").
4. Class = class of the median. If the covering orthologs disagree on class, label the
   segment `mixed`.
5. Merge neighbouring segments with the same label and the same set of covering orthologs.
   The merged o/e is the length-weighted mean of the segment medians.
6. Pure function, no I/O, unit-testable on hand-built inputs (ask before adding tests).

### 5. Output format

Main output `{run}.GNOMAD_rmc.txt`, the standard 4-column range format:

```
GnomAD_RMC	743	784	constrained|oe=0.06|n=2/7|MYH3:0.06,MYH7:0.05
```

The description is pipe-delimited, following the `hydrophobic_patch` precedent
(`_parse_patch_description`):
- field 1 = class label, the only part used for colour and legend;
- `oe` = merged median;
- `n` = orthologs covering / orthologs with RMC;
- then per-ortholog values for the hover text.

Sidecar `{run}.GNOMAD_rmc.json`:
- Alliance rows and orthologs used, with each one's status (`ok`, `no_rmc`,
  `no_mane_transcript`, `numbering_mismatch`, `error: ...`);
- the per-ortholog projected regions before the merge;
- the alignment identity and the aligned fraction of the worm protein.

### 6. Registration

- `task_definitions.json`: new `gnomad` task (`script: gnomad_rmc.py`,
  `result_type: range`, `output_suffix: .txt`, `selectable: true`, `accession` param with a
  tooltip). Not in `default_tasks` until it has been used on real runs.
- `scripts/task_runners.py`: `run_gnomad(args, report=None, ...)`, then add it to the runner map.
- `run_tag_sites_from_json.py` needs no change if it dispatches by the registry. Confirm the
  subprocess and in-process paths give byte-identical output (dual-use rule).

### 7. Plot — `utils/results.py`

- Add `"GnomAD_RMC"` to `_FEATURE_SOURCES` and give it the next free row in `y_positions`.
- In `_annotation_color_map` / `plot_results`, key the colour on the class label only
  (`desc.split("|")[0]`) and use a fixed class palette rather than the rotating one:
  constrained dark, intermediate mid, tolerant light, mixed hatched or grey.
- Hover text from the parsed fields: `Constrained (o/e 0.06, 2 of 7 orthologs): MYH3 0.06,
  MYH7 0.05`.
- Check that other consumers of `range_df` (`residue_colors_for_annotations`,
  `residues_for_sidechains`, reagents ranking) ignore or handle the new source.

### 8. Docs

- `CLAUDE.md`:
  - add `gnomad_rmc.py` to the layout;
  - add gnomAD GraphQL and Ensembl GRCh37 REST under External services;
  - note that the orthology table must be regenerated per Alliance release.
- Script-level usage notes go in the script docstring, not CLAUDE.md.

## Verification

1. Run the six test genes from the prototype (sma-4, let-60, trxr-1, dbl-1, unc-54, par-2)
   plus trp-4 via the CLI. Expected:
   - unc-54 shows a constrained block near worm 743-784;
   - sma-4 shows MH1 and MH2 constrained with an intermediate linker;
   - par-2 and trp-4 produce empty files and say why in the sidecar.
2. The per-ortholog regions in the sidecar match the prototype's `GNOMAD_regions.txt`.
3. The same genes through the app: one `GnomAD_RMC` row, three or four legend entries,
   correct hover text, and no event-loop stall while the task runs.
4. In-process vs. subprocess outputs are byte-identical.
5. `./deploy.sh --manifest` includes the new orthology table and nothing large.

## Open questions

- Whether `mixed` segments need their own colour or just a hatch over the median's class.
- Whether to show segments covered by only one of several orthologs, or require a minimum
  fraction of orthologs.
- Cache location for the deployed app (run folder, or a shared cache across runs).
