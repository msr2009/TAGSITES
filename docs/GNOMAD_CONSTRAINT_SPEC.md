# Spec: human-population constraint track (gnomAD via Alliance orthologs)

Status: feasibility test, not yet a pipeline task. Goal: decide whether a worm-residue track
derived from human gnomAD missense depletion is worth adding, and at what resolution
(per-residue vs. ranges).

## 1. Question the test must answer

For a worm protein, which residues sit in regions where the human ortholog tolerates few
missense variants? Three outputs are compared:

| Tier | Source | Resolution | Expected value |
|---|---|---|---|
| A. Gene | `gnomad_constraint.oe_mis`, `mis_z` | whole protein | context badge only |
| B. Regions | gnomAD v2 regional missense constraint (RMC) | 2-10 ranges / gene | "loose ranges", robust |
| C. Windows | computed from the gnomAD v4 variant list | per-residue, smoothed | "tight", noisy |

The test passes if B and/or C, projected onto worm coordinates, enrich for known functional
positions (section 6) beyond what the existing BLAST conservation track already gives.

## 2. Inputs (confirmed 2026-10-15)

**Orthology**: `~/Downloads/ORTHOLOGY-ALLIANCE_TSV_COMBINED.tsv.gz` (Alliance 9.1.0,
"Stringent" filter). Worm-human rows: 15,460 pairs per direction, 7,432 worm genes, of which
3,083 have more than one human ortholog. Columns used: `Gene1ID` (WB:WBGene...),
`Gene2ID` (HGNC:...), `Gene2Symbol`, `AlgorithmsMatch`, `IsBestScore`, `IsBestRevScore`.
The older `ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz` (Aug) is superseded.

**gnomAD**: GraphQL at `https://gnomad.broadinstitute.org/api`. One `gene` query returns
everything needed (verified on SMAD4, MYH7, KRAS):

```graphql
gene(gene_symbol: "SMAD4", reference_genome: GRCh38) {
  hgnc_id canonical_transcript_id mane_select_transcript { ensembl_id }
  gnomad_constraint { oe_mis mis_z pli }
  variants(dataset: gnomad_r4) { variant_id consequence hgvsp transcript_id
                                 exome { ac an } genome { ac an } }
  clinvar_variants { variant_id clinical_significance major_consequence hgvsp }
}
# RMC is only on GRCh37 (returns null on GRCh38):
gene(gene_symbol: "SMAD4", reference_genome: GRCh37) {
  gnomad_v2_regional_missense_constraint { passed_qc has_no_rmc_evidence
    regions { aa_start aa_stop obs_mis exp_mis obs_exp chisq_diff_null } }
}
```

SMAD4 example: 424 missense, 261 synonymous variants (v4); RMC splits it into
Met1-Ser150 (o/e 0.20), Ser150-Asn316 (0.68), Asn316-Asp552 (0.20), i.e. MH1 / linker / MH2.

**Human protein sequence**: Ensembl REST `/sequence/id/{ENST}?type=protein` for the MANE
transcript (v4 variants) and for the GRCh37 canonical transcript (RMC), via the existing
`scripts/ensembl_rest.py`.

**Worm protein -> WBGene**: UniProt entry cross-reference (`scripts/uniprot_api.py`
already fetches entries), or WormBase gene ID when the input came from the local store.

## 3. Pipeline under test (prototype script: `scripts/gnomad_constraint.py`)

1. **Resolve orthologs.** Look up the WBGene in the Alliance file. Keep pairs with
   `IsBestScore == IsBestRevScore == Yes`; if none, keep the max `AlgorithmsMatch` set.
   Keep all ties: one-to-many is common (unc-54 -> 10 MYH genes; dbl-1 -> BMP2/4/6/7).
2. **Fetch gnomAD** for each human ortholog (two queries, above). Cache raw JSON per gene
   under the run folder; gnomAD rate-limits, so requests go through `http_retry.py` with a
   1-2 s spacing and, in the app, off the event loop (`_off_loop`).
3. **Check human residue numbering.** Parse `hgvsp` (`p.Arg361His`) for missense variants on
   the MANE transcript; the reference residue must match the fetched human sequence at that
   position for 100% of variants, else the gene fails (same lesson as `cds_check`). Do the
   same for RMC `aa_start`/`aa_stop` against the GRCh37 transcript protein. Note RMC regions
   come in genomic order, so minus-strand genes list them reversed (MYH7: `Glu1935 -> Ala850`,
   ...); normalise to `min/max`.
4. **Align worm to human** pairwise with Biopython `PairwiseAligner` (global, BLOSUM62,
   open -10 / extend -0.5); no network needed. Build a residue map `worm_pos -> human_pos`
   (None in gaps). Record per-position local identity (e.g. 11-aa window).
5. **Per-residue human score (tier C).** For human residue i, in a window of w residues
   (test w = 11, 21, 31):
   `mis_rate = distinct missense alleles / w`, normalised by the same ratio for synonymous
   alleles in the window (MTR-style; synonymous counts act as a crude mutability control).
   Report `depletion = window_ratio / gene_ratio` (<1 = constrained). Also try allele-count
   weighting vs. distinct-allele counting; singleton-only is a third variant.
6. **Project to worm.**
   - *Tight*: worm residue gets the human score of its aligned partner only if aligned and
     local identity >= threshold (test 20%, 30%, 40%); else NaN.
   - *Loose*: RMC regions projected through the map: worm range = min/max worm positions whose
     partners fall in the human region; regions covering < 50% aligned worm residues are
     flagged low-confidence.
7. **Combine paralogs.** With several human orthologs, report each separately and a
   consensus: median depletion across orthologs per worm residue, plus a count of how many
   call it constrained. Pick the rule after looking at unc-54 and dbl-1.
8. **Classify.** Worm residue -> `constrained` (o/e or depletion <= 0.4),
   `intermediate`, `tolerant` (>= 0.8), or `unmapped`. Thresholds are a test variable.

Output files match existing track formats so `utils/results.py` can load them unchanged:
- `{run}.GNOMAD_depletion.txt` — continuous, `pos<TAB>value` (NaN where unmapped)
- `{run}.GNOMAD_regions.txt` — range, `GnomAD_RMC<TAB>start<TAB>end<TAB>o/e=0.20 (SMAD4)`
- `{run}.GNOMAD_map.json` — ortholog choice, sequence-check results, alignment identity

## 4. Test genes

| Worm | Human | Case covered |
|---|---|---|
| sma-4 | SMAD4 | clean 1:1, strong RMC (MH1/MH2), ClinVar hotspots (R361) |
| let-60 | HRAS, KRAS, NRAS | 1:3, very high identity, small protein (few variants: noise floor) |
| trxr-1 | TXNRD1, TXNRD3 | selenocysteine (`U`) near C-term: alignment and hgvsp edge case |
| dbl-1 | BMP2/4/6/7, NODAL | 1:many, low identity outside mature domain: tests the loose fallback |
| unc-54 | MYH1-8, 13, 15, 7B | 1,966 aa, many-to-many, minus-strand RMC, large payload |
| (negative) | — | a nematode-specific gene absent from the file (e.g. an srh/srx receptor) |

The `worm_human_orthologs_noGOterm.tsv` list (667 genes) is a pool for a later scale test.

## 5. Mapping checks (must pass)

- Human `hgvsp` reference-residue match: 100% per gene.
- RMC boundary residues match the GRCh37 protein; report whether GRCh37 canonical and MANE
  proteins are identical (if not, RMC needs its own human->human map).
- Fraction of worm residues mapped, per gene and per ortholog.
- Alliance ortholog vs. top human hit from the app's own BLAST task: agreement rate. This
  settles whether Alliance is needed or the existing BLAST hits suffice.

## 6. Signal checks (decide go / no-go)

- **ClinVar positive control**: human residues with Pathogenic/Likely-pathogenic missense
  variants, projected to worm. Measure enrichment (odds ratio, AUC) of those positions in
  `constrained` vs. other classes, tier B vs. tier C, per window size.
- **Worm-side control**: known worm loss-of-function missense alleles and UniProt
  active/binding sites already in `{run}.UNIPROT_uniprot.txt`.
- **Redundancy**: Spearman correlation with the existing BLAST conservation track. A high
  correlation means little added value; the useful case is constraint that conservation
  misses (e.g. conserved-but-tolerant surface loops, or human-specific constrained motifs).
- **Stability**: how much the classification changes across window size and identity
  threshold; prefer the setting with the fewest flips at equal enrichment.
- **Tag-site relevance**: for each test gene, check that N/C-termini and loops chosen by
  the existing tag-site scoring fall in `tolerant` or `unmapped`, not `constrained`.

## 7. Operational checks

- Time and bytes per gene (MYH7 is the worst case); response time when gnomAD is uncached.
- Behaviour on gnomAD errors, rate limiting, and genes with `passed_qc: false` or
  `has_no_rmc_evidence: true` (fail soft: track absent, run continues).
- Proteome scale (7,432 worm genes, ~12k human genes): the API is not suitable; the batch path
  would need the gnomAD v4 exome sites table (or a pre-extracted missense/synonymous per-codon
  table) and the RMC download. Estimate size before committing, as with `build_proteome_db`.

## 8. Integration if it passes

- New task `gnomad` in `task_definitions.json` (`result_type: continuous`, plus regions as a
  range track), script runnable standalone per the dual-use rule.
- Display-and-ranking only: constraint annotates sites and never filters a site or guide
  (same rule as RS3/conservation).
- Alliance file: ship a pre-filtered worm-human TSV (~1 MB) rather than the 21 MB combined
  file; record the Alliance version in the map JSON.

## Open decisions

1. Consensus rule for one-to-many orthologs (median vs. best-identity ortholog).
2. Whether tier C beats RMC enough to justify per-residue display, or ranges alone ship.
3. Whether to extend beyond human (Alliance covers mouse/fly/zebrafish, but only human has
   population data at this depth).

## Prototype results (2026-10-15, `scripts/gnomad_constraint.py`)

Settings: MTR window 31, local identity >= 0.3 for per-residue MTR, MTR classes by
within-gene quartile, RMC classes absolute (<= 0.4 / >= 0.8). AUC = how well each signal ranks
ClinVar P/LP missense partner residues above all other scored worm residues (0.5 = chance).
"identity" is the worm/human local identity, a stand-in for the BLAST conservation track.

| Worm | Orthologs (rule) | P/LP res. | AUC MTR | AUC RMC | AUC identity |
|---|---|---|---|---|---|
| sma-4 | SMAD4 | 19 | 0.57 | 0.38* | 0.69 |
| let-60 | HRAS (recip.) | 12 | 0.61 | 0.64 | 0.63 |
| let-60 | HRAS, KRAS, NRAS (all) | 39 | 0.65 | 0.64 | 0.67 |
| dbl-1 | BMP2/4/7 (recip.) | 9 | 0.61 | 0.64 | 0.70 |
| dbl-1 | + BMP6, NODAL (all) | 11 | 0.63 | 0.64 | 0.75 |
| unc-54 | MYH1, MYH3 (recip.) | 39 | 0.64 | 0.73 | 0.71 |
| unc-54 | 9 MYH genes (all) | 266 | 0.63 | 0.68 | 0.56 |
| trxr-1 | TXNRD1, TXNRD3 | 0 | — | — | — |
| par-2, trp-4 | none in Alliance | — | — | — | — |

\* SMAD4 has 3 RMC regions with near-equal o/e (0.196 vs 0.202) for MH1 and MH2, so tiny
differences decide the ranking; 95% of P/LP residues do sit in the constrained class.

Findings:
- Numbering checks passed everywhere: 0 hgvsp reference mismatches over 36,876 variants in 20 human genes, all
  RMC boundaries matched, and every CDS translated (after allowing Sec and the CTG start of TXNRD3).
- Only unc-54 + MYH7 (266 P/LP) is large enough to be decisive; there gnomAD beats sequence
  identity (RMC 0.68, MTR 0.63 vs 0.56). With few P/LP residues, identity is as good or better.
  ClinVar P/LP labels are biased toward conserved residues, so this test favours identity.
- RMC is at least as good as per-residue MTR in 5 of 6 evaluable runs, and gives clean
  ranges (unc-54: MYH3 740-776, o/e 0.06, a 40-residue block; SMAD4 MH1/linker/MH2).
- Dropping singletons (`--min_ac 2`) raised MTR AUC slightly for let-60, dbl-1, unc-54
  (+0.02-0.03) and lowered it for sma-4 (0.57 -> 0.44). Not adopted.
- The reciprocal-best rule drops the most informative paralogs (KRAS, MYH7). Use all Alliance
  orthologs with a median consensus.
- RMC coverage is patchy: absent for TXNRD1/3, MYH1, MYH13, MYH15, NRAS, NODAL. MYH7B has no
  MANE transcript. For KRAS the GRCh37 and MANE proteins differ, so RMC is re-mapped.
- Worm termini are mostly unmapped or tolerant, so the track rarely argues against a terminal tag.
- Cost: 20-100 s per gene uncached (gnomAD rate-limit spacing); unc-54 with 10 orthologs 340 s.
  Cached reruns < 1 s. MYH payloads are several MB each.
- par-2 has no orthologs in any Alliance species. trp-4 (fly NOMPC only) was added as a
  second negative; both exit cleanly with `no_human_ortholog`.
