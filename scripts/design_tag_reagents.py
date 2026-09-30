"""
design_tag_reagents.py

For every potential tag-insertion site in a protein, identify the N nearest
CRISPR guide RNAs and construct HDR homology arms with PAM-disrupting mutations.

Pipeline:
  1. Parse Genewise output → CDS exon map
  2. Enumerate all insertion sites (one per protein residue)
  3. Find all guide cut sites in the genomic sequence (both strands)
  4. For each (residue × guide) pair within arm_length of the cut:
       - Build left and right homology arms centered on the insertion site
       - Block re-cutting with protein-preserving edits only. Ladder:
         (a) insertion within the PAM seed already blocks re-cutting;
         (b) 1 synonymous codon change disrupting the PAM;
         (c) 1 silent intronic base change within the PAM;
         (d) synonymous codon change(s) placing >=2 mismatches in the 15 nt
             protospacer seed proximal to the PAM.
         Guides that can't be blocked without a non-synonymous change (or whose
         PAM lies outside both arms) are discarded, and the next-closest guide
         is tried in their place.

Output TSV (one row per residue × guide):
  residue_index     1-based protein residue number
  amino_acid        residue identity
  insert_pos        0-based genomic insertion coordinate
  exon_index        exon containing (or adjacent to) this residue
  is_split_codon    True if codon straddles an exon–intron boundary
  dist_to_5p_splice distance to 5′ exon boundary (bp)
  dist_to_3p_splice distance to 3′ exon boundary (bp; -1 for last exon)
  guide_strand      '+' or '-'
  spacer            guide spacer sequence
  pam_seq           PAM sequence on guide strand
  pam_fwd_start     0-based start of PAM in fwd genomic coords
  cut_pos           0-based DSB position in fwd coords
  distance          |cut_pos - insert_pos| (bp)
  pam_in_arm        'left' / 'right' / 'both' / 'none'
  recut_block_method  'insertion' / 'syn_1' / 'mut_1' / 'syn_seed'
  mutation_desc     human-readable description of the PAM mutation
  left_arm          left homology arm: exonic bases uppercase, intronic lowercase
  right_arm         right homology arm: exonic bases uppercase, intronic lowercase
  left_arm_wt       left arm before PAM disruption (WT); same case encoding
  right_arm_wt      right arm before PAM disruption (WT); same case encoding
  rs3_score         Rule Set 3 on-target cutting score (blank if unavailable).
                    Display-only: never affects guide choice or ordering.
  rs3_percentile    rank percentile of this guide's RS3 score among all
                    candidate guides in this region (blank if unavailable)
  offtarget_count   PAM-bearing near-matches found elsewhere in the genome.
                    BLANK (not 0) when no screen ran, so an unscreened guide is
                    never mistaken for a clean one. Display-only.
  offtarget_identical  of those, how many are perfect matches with a PAM — the
                    dangerous case, since a guide cuts both copies equally
  offtarget_detail  human-readable breakdown, e.g. '2 sites (1 identical, 1 near)'
  offtarget_status  'screened' / 'not_checked' / 'failed' / 'pending'

The companion .genotyping.tsv additionally carries offtarget_amplicons (count of
predicted spurious products where BOTH primers bind one duplicated locus) and
offtarget_detail (their accessions and sizes; '*' marks a perfect match).

A companion <output>.genotyping.tsv is always written (genomic sequence is
always available): one row per residue x amplicon_type. 'external' alone
when --insert_sequence is empty; plus '5p_junction'/'3p_junction' when a tag
sequence is given. Columns: residue_index, amino_acid, insert_pos,
amplicon_type, fwd_seq, fwd_tm, rev_seq, rev_tm, product_size.
See reagent_sequences.design_genotyping_primers.

Matt Rich, 2025
"""

import sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from crispr_util import find_guides, build_frame_lookup, disrupt_pam
from parse_genewise import parse_genewise, enumerate_insertion_sites, \
    parse_genewise_score, parse_genewise_gff_score, cds_coverage
from guide_efficiency import guide_key, rs3_percentiles, score_guides
from offtarget_screen import (
    load_config as load_offtarget_config,
    predict_amplicons,
    screen_guide,
    summarise as offtarget_summarise,
)
from progress import report as _report, resolve_reporter
from reagent_sequences import annotate_region_spans, design_genotyping_primers


def _case_arm(arm_seq, arm_start, frame_lookup):
    """Return arm_seq with exonic (coding) positions uppercase, intronic lowercase."""
    return ''.join(
        ch.upper() if (arm_start + i) in frame_lookup else ch.lower()
        for i, ch in enumerate(arm_seq)
    )


def _design_genotyping_rows(sites, dna, L, insert_sequence, internal_threshold,
                            primer_opt_tm, product_opt_size, flank_min, flank_max,
                            offtarget_hsps=None, offtarget_cfg=None):
    """One row per (residue x amplicon_type) genotyping primer pair.

    Pulls flanking sequence directly from the full genomic record around each
    site's insert_pos — independent of arm_length, since genotyping primers may
    need to sit outside the (possibly much shorter) homology arms.
    """
    margin = 60   # extra room beyond flank_max so primer3 has a real window to search
    w = flank_max + margin
    rows = []
    for _, site in sites.iterrows():
        insert_pos = int(site['insert_pos'])
        left_flank  = dna[max(0, insert_pos - w):insert_pos]
        right_flank = dna[insert_pos:min(L, insert_pos + w)]

        primers = design_genotyping_primers(
            left_flank, right_flank, insert_sequence,
            internal_threshold=internal_threshold,
            primer_opt_tm=primer_opt_tm,
            product_opt_size=product_opt_size,
            flank_min=flank_min, flank_max=flank_max,
        )
        # place each primer in genomic coordinates so it can be screened for
        # co-amplification; a primer inside the insert maps to None and is skipped
        annotate_region_spans(primers, len(left_flank), len(insert_sequence),
                              insert_pos - len(left_flank))
        for amplicon_type, p in primers.items():
            amps = []
            if offtarget_hsps and p.get('fwd_region_span') and p.get('rev_region_span'):
                amps = predict_amplicons(offtarget_hsps, p['fwd_region_span'],
                                         p['rev_region_span'], offtarget_cfg)
            rows.append({
                'residue_index': int(site['residue_index']),
                'amino_acid':    site['amino_acid'],
                'insert_pos':    insert_pos,
                'amplicon_type': amplicon_type,
                'fwd_seq':       p['fwd_seq'],
                'fwd_tm':        p['fwd_tm'],
                'rev_seq':       p['rev_seq'],
                'rev_tm':        p['rev_tm'],
                'product_size':  p['product_size'],
                'offtarget_amplicons': len(amps),
                'offtarget_detail':    _amplicon_detail(amps),
            })
    return pd.DataFrame(rows)


def _exon_intervals(cds_df):
    """Annotated exon intervals in region coords, for transcript-vs-genomic hit calling."""
    if cds_df is None or len(cds_df) == 0:
        return []
    return [(int(r['start']), int(r['stop']) + 1) for _, r in cds_df.iterrows()]


def _offtarget_preconditions(taxid, email, reporter):
    """True when an off-target screen can run at all (needs a species and an email)."""
    if not taxid or str(taxid).strip() in ('1', '1.0', 'None'):
        _report(reporter, 'Off-target screen skipped: no species taxid for this run',
                stage='offtarget', level='warning')
        return False
    if not email:
        _report(reporter, 'Off-target screen skipped: EBI submissions require an email',
                stage='offtarget', level='warning')
        return False
    return True


def _run_region_screen(dna, taxid, email, cds_df, cfg, reporter, job_id_cb, resume_job_ids):
    """Screen A. Returns (status, region_result, pending_sentinel)."""
    if not _offtarget_preconditions(taxid, email, reporter):
        return 'not_checked', {}, None

    # imported here rather than at module top so the CLI still runs when no network
    # stack is available, matching how the other remote backends are reached
    import offtarget_remote

    try:
        result = offtarget_remote.run_region_screen(
            dna, email, taxid, exons=_exon_intervals(cds_df), cfg=cfg, report=reporter,
            job_id_cb=job_id_cb, resume_job_ids=resume_job_ids)
    except Exception as e:
        # A failed screen must never fail the whole reagents run, and must never
        # look like "no off-targets found" — the status stays unscreened
        _report(reporter, 'Off-target region screen failed ({}: {}); guides will show as '
                          'not checked'.format(type(e).__name__, e),
                stage='offtarget', level='warning')
        return 'failed', {}, None
    if isinstance(result, dict) and 'ebi_status' in result:
        return 'pending', {}, result
    return 'screened', result, None


def _run_spacer_screen(spacers, taxid, email, pam, cfg, reporter, job_id_cb,
                       resume_job_ids, self_spans=None, own_accessions=None):
    """Screen B. Returns (spacer_result, pending_sentinel)."""
    import offtarget_remote
    try:
        result = offtarget_remote.run_spacer_screen(
            spacers, email, taxid, pam=pam, cfg=cfg, report=reporter,
            job_id_cb=job_id_cb, resume_job_ids=resume_job_ids, self_spans=self_spans,
            own_accessions=own_accessions)
    except Exception as e:
        _report(reporter, 'Off-target spacer screen failed ({}: {}); only the region '
                          'screen contributes'.format(type(e).__name__, e),
                stage='offtarget', level='warning')
        return {}, None
    if isinstance(result, dict) and 'ebi_status' in result:
        return {}, result
    return result, None


def _amplicon_detail(amps):
    """Compact description of predicted spurious products for the genotyping TSV."""
    if not amps:
        return ''
    parts = ['{}{} ~{}bp'.format(a['acc'], '*' if a['perfect'] else '', a['product_size'])
             for a in amps[:3]]
    if len(amps) > 3:
        parts.append('+{} more'.format(len(amps) - 3))
    return '; '.join(parts)


# ── Core logic ────────────────────────────────────────────────────────────────

def design_reagents(
    genewise_out,
    genomic_fasta,
    protein_length=None,
    n_guides=5,
    arm_length=1000,
    pam='NGG',
    guide_length=20,
    cut_offset=3,
    insert_sequence='',
    internal_threshold=500,
    primer_opt_tm=60.0,
    product_opt_size=200,
    flank_min=50,
    flank_max=150,
    rs3=True,
    rs3_tracr='Hsu2013',
    offtarget=True,
    taxid='',
    email='',
    offtarget_sidecar='',
    report=None,
    job_id_cb=None,
    resume_job_ids=None,
):
    """
    Full pipeline: Genewise output + genomic FASTA → reagent table (DataFrame).

    Parameters
    ----------
    genewise_out   : str  path to Genewise .out.txt
    genomic_fasta  : str  path to the genomic region FASTA
    protein_length : int or None  amino acid count of the query protein; used to
                     check CDS coverage.  If None, coverage check is skipped.
    n_guides       : int  max guides to report per insertion site
    arm_length     : int  homology arm length (bp) on each side
    pam            : str  IUPAC PAM (default 'NGG')
    guide_length   : int  spacer length (default 20)
    cut_offset     : int  cut distance from PAM (SpCas9=3)
    insert_sequence : str  default insert/tag DNA sequence. Genotyping primers are
                     always designed (external pair only when empty; plus 5'/3'
                     junction pairs when a tag sequence is given). See
                     reagent_sequences.design_genotyping_primers.
    internal_threshold, primer_opt_tm, product_opt_size, flank_min, flank_max :
                     passed through to design_genotyping_primers.
    rs3            : bool  compute Rule Set 3 on-target scores (default True).
                     Display-only: scores never affect which guides are kept or
                     their order. Silently blank if the optional rs3 package is
                     unavailable, or if pam/guide_length aren't SpCas9 NGG/20.
    rs3_tracr      : str  tracrRNA scaffold assumed by RS3 ('Hsu2013'/'Chen2013').

    Returns
    -------
    (df, genotyping_df) : genotyping_df is a DataFrame of genotyping primer
                     pairs, one row per residue x amplicon_type.

    Raises
    ------
    ValueError if the Genewise alignment looks like a wrong or truncated
    genomic sequence (score < 50 bits OR CDS coverage < 90%), or if any
    homology arm in the output is shorter than half of arm_length.
    """
    reporter = resolve_reporter(report)

    # 1. Load genomic sequence
    records = list(SeqIO.parse(genomic_fasta, 'fasta'))
    if not records:
        raise ValueError('No sequences in {}'.format(genomic_fasta))
    dna = str(records[0].seq)
    L = len(dna)

    # 2. Parse Genewise → CDS exon table
    cds_df = parse_genewise(genewise_out)
    _report(reporter, '{} CDS exons parsed'.format(len(cds_df)), stage='parse_genewise')

    # Bad-alignment guard: raise an error rather than silently producing garbage
    # reagents.  The EBI REST API does not write a "Score NNN bits" header line,
    # so we read the score from the GFF match row instead.
    score = parse_genewise_gff_score(genewise_out)
    if score is None:
        score = parse_genewise_score(genewise_out)
    coverage = None
    if protein_length and protein_length > 0:
        coverage = cds_coverage(cds_df, protein_length)

    bad_score    = score    is not None and score    < 50.0
    # 90% threshold: catches a partial isoform (e.g. 3-exon short DNA aligned against
    # a 5-exon long protein = 81.6% coverage) while passing correct full-gene alignments
    bad_coverage = coverage is not None and coverage < 0.90

    if bad_score or bad_coverage:
        parts = []
        if score is not None:
            parts.append('score {:.1f} bits'.format(score))
        if coverage is not None:
            parts.append('CDS coverage {:.0%}'.format(coverage))
        detail = ', '.join(parts)
        raise ValueError(
            'Mismatch between DNA and protein sequences ({}).\n'
            'Confirm genomic sequence covers (1) correct gene isoform and '
            '(2) contains sufficient flanking sequence to make full-length '
            'homology arms (>{} bp).'.format(detail, arm_length)
        )

    # 3. Enumerate insertion sites
    sites = enumerate_insertion_sites(cds_df, dna)
    _report(reporter, '{} residues / insertion sites'.format(len(sites)), stage='sites')

    # 4. Build per-position frame lookup for synonymous mutation
    frame_lookup = build_frame_lookup(cds_df, dna)

    # 5. Find all guide cut sites
    guides = find_guides(dna, pam=pam, guide_length=guide_length,
                         cut_offset=cut_offset)
    _report(reporter, '{} guide sites found (PAM={})'.format(len(guides), pam), stage='guides')

    # 5b. RS3 on-target scores, computed once for the whole region before the per-site
    # loop: find_guides() has already enumerated every candidate, and many insertion
    # sites reuse the same guide, so cost is independent of the number of sites.
    # Display-only — nothing below selects or orders guides by these values.
    rs3_scores = {}
    rs3_pct = {}
    if rs3:
        rs3_scores, rs3_note = score_guides(dna, guides, pam=pam,
                                            guide_length=guide_length, tracr=rs3_tracr)
        rs3_pct = rs3_percentiles(rs3_scores)
        if rs3_scores:
            _report(reporter, 'RS3 scored {} guides (tracr={})'.format(
                len(rs3_scores), rs3_tracr), stage='guides')
        if rs3_note:
            _report(reporter, rs3_note, stage='guides')

    # 5c. Off-target / primer-specificity screen. Two EBI blastn jobs: the whole
    # region (duplicated segments -> primer co-amplification and guides in repeats)
    # and all spacers concatenated into one query (scattered guide near-matches).
    # Annotation only: nothing below selects or orders guides by these results.
    ot_cfg = load_offtarget_config()
    ot_hsps = []
    ot_region = {}
    ot_status = 'not_checked'
    if offtarget:
        ot_status, ot_region, pending = _run_region_screen(
            dna, taxid, email, cds_df, ot_cfg, reporter, job_id_cb, resume_job_ids)
        if pending:
            return pending
        ot_hsps = ot_region.get('duplicates', [])

    pam_len = len(pam)
    rows = []

    for _, site in sites.iterrows():
        insert_pos = int(site['insert_pos'])

        # Sort guides by distance of their cut to the insertion site
        scored = sorted(guides, key=lambda g: abs(g['cut_pos'] - insert_pos))

        # Keep up to n_guides guides (nearest first) whose re-cutting can be
        # blocked without changing the protein.  Guides that can't be blocked
        # protein-preservingly are discarded, so we keep scanning past them
        # until n_guides survive (or we run out of guides within arm_length).
        _SEED = 15
        n_kept = 0
        for g in scored:
            if abs(g['cut_pos'] - insert_pos) > arm_length:
                break   # sorted ascending → no closer guides remain
            if n_kept >= n_guides:
                break

            cut_pos       = g['cut_pos']
            pam_fwd_start = g['pam_fwd_start']
            distance      = abs(cut_pos - insert_pos)

            # Build raw arms
            left_start  = max(0, insert_pos - arm_length)
            right_end   = min(L, insert_pos + arm_length)
            left_arm_raw  = dna[left_start:insert_pos]
            right_arm_raw = dna[insert_pos:right_end]

            # Determine which arm(s) contain the PAM
            pam_end = pam_fwd_start + pam_len
            pam_in_left  = pam_end > left_start and pam_fwd_start < insert_pos
            pam_in_right = pam_fwd_start < right_end and pam_end > insert_pos

            if pam_in_left and pam_in_right:
                pam_arm = 'both'
            elif pam_in_left:
                pam_arm = 'left'
            elif pam_in_right:
                pam_arm = 'right'
            else:
                pam_arm = 'none'

            # Check whether the tag insertion itself disrupts re-cutting.
            # Any insertion within the 15 bp protospacer seed region immediately
            # 5' of the PAM (on the guide strand) prevents Cas9 from re-binding.
            # PAM mutations are only needed when the insert falls outside this window.
            if g['strand'] == '+':
                # seed region = 15 bp 5' of PAM + PAM itself (fwd coords)
                insertion_blocks = (pam_fwd_start - _SEED <= insert_pos < pam_end)
            else:
                # seed region = PAM + 15 bp 3' of PAM (= 15 bp 5' on − strand)
                insertion_blocks = (pam_fwd_start <= insert_pos < pam_end + _SEED)

            recut_method  = 'none'
            mutation_desc = ''
            left_arm  = left_arm_raw
            right_arm = right_arm_raw

            if insertion_blocks:
                recut_method  = 'insertion'
                mutation_desc = 'insert within {} bp seed region of PAM; re-cutting impossible'.format(_SEED)

            elif pam_arm == 'none':
                # PAM lies outside both homology arms → no protein-preserving edit
                # can protect this allele from re-cutting.  Discard the guide.
                continue

            else:
                # disrupt_pam works on the full genomic sequence and returns a
                # modified copy; we then re-extract the arms from it.  A None
                # result means re-cutting can't be blocked without a
                # non-synonymous change → discard the guide.
                result = disrupt_pam(dna, pam, pam_fwd_start,
                                     g['strand'], frame_lookup, seed_len=_SEED)
                if result is None:
                    continue
                mutated_seq, desc, method = result
                left_arm  = mutated_seq[left_start:insert_pos]
                right_arm = mutated_seq[insert_pos:right_end]
                recut_method  = method
                mutation_desc = desc

            # encode exonic vs intronic in arm case (uppercase = coding, lowercase = intronic)
            left_arm  = _case_arm(left_arm, left_start, frame_lookup)
            right_arm = _case_arm(right_arm, insert_pos, frame_lookup)

            # also emit the unmutated (WT) arms so the reagents display can show a
            # true wild-type row and highlight every mutated base by diffing — not
            # just the PAM positions (a synonymous PAM-disrupting codon change often
            # alters a base outside the 3 bp PAM).
            left_arm_wt  = _case_arm(left_arm_raw, left_start, frame_lookup)
            right_arm_wt = _case_arm(right_arm_raw, insert_pos, frame_lookup)

            gkey = guide_key(g)
            # Combine both screens for this guide: duplicated-segment sites from the
            # region query plus scattered near-matches from the spacer query
            ot_n = ot_ident = ot_unver = 0
            if ot_status == 'screened':
                region_hit = screen_guide(ot_hsps, {'pam_fwd_start': pam_fwd_start,
                                                    'guide_strand': g['strand']},
                                          guide_length, pam, ot_cfg)
                ot_n = region_hit['n_total']
                ot_ident = region_hit['n_identical']
                ot_unver = region_hit['n_pam_unverified']
            rows.append({
                # carried so the spacer screen, which runs once the kept guides are
                # known, can add its counts without re-deriving the guide identity
                '_gkey':               gkey,
                '_spacer_pam':         g['spacer'] + g['pam_seq'],
                '_ot_unver':           ot_unver,
                'residue_index':       int(site['residue_index']),
                'amino_acid':          site['amino_acid'],
                'insert_pos':          insert_pos,
                'exon_index':          int(site['exon_index']),
                'is_split_codon':      site['is_split_codon'],
                'dist_to_5p_splice':   int(site['dist_to_5p_splice']),
                'dist_to_3p_splice':   int(site['dist_to_3p_splice']),
                'guide_strand':        g['strand'],
                'spacer':              g['spacer'],
                'pam_seq':             g['pam_seq'],
                'pam_fwd_start':       pam_fwd_start,
                'cut_pos':             cut_pos,
                'distance':            distance,
                'pam_in_arm':          pam_arm,
                'recut_block_method':  recut_method,
                'mutation_desc':       mutation_desc,
                'left_arm':            left_arm,
                'right_arm':           right_arm,
                'left_arm_wt':         left_arm_wt,
                'right_arm_wt':        right_arm_wt,
                'rs3_score':           rs3_scores.get(gkey, ''),
                'rs3_percentile':      rs3_pct.get(gkey, ''),
                # blank, not 0, when nothing was screened — a 0 would read as "clean"
                'offtarget_count':     ot_n if ot_status == 'screened' else '',
                'offtarget_identical': ot_ident if ot_status == 'screened' else '',
                'offtarget_detail':    offtarget_summarise(ot_n, ot_ident, ot_unver)
                                       if ot_status == 'screened' else '',
                'offtarget_status':    ot_status,
            })
            n_kept += 1

    df = pd.DataFrame(rows)

    # 5d. Screen B runs HERE, not before the loop: only now are the guides that
    # actually reach the output known. Screening every candidate in the region would
    # be a ~4x longer query (measured: 570 candidates vs 131 kept for snt-1) for hits
    # on guides the user never sees.
    ot_spacer = {}
    if ot_status == 'screened' and not df.empty:
        wanted = {}
        for gkey, spacer_pam in zip(df['_gkey'], df['_spacer_pam']):
            wanted.setdefault(gkey, spacer_pam)
        spacer_result, pending = _run_spacer_screen(
            list(wanted.values()), taxid, email, pam, ot_cfg, reporter,
            job_id_cb, resume_job_ids, self_spans=ot_region.get('self_spans'),
            own_accessions=ot_region.get('own_accessions'))
        if pending:
            return pending
        ot_spacer = spacer_result
        hits = spacer_result.get('spacer_hits') or {}
        order = list(wanted)
        by_gkey = {order[i]: v for i, v in hits.items() if i < len(order)}
        # fold the spacer counts into the per-row totals from the region screen
        df['offtarget_count'] = [
            c + len(by_gkey.get(k, [])) if c != '' else c
            for c, k in zip(df['offtarget_count'], df['_gkey'])]
        df['offtarget_identical'] = [
            c + sum(1 for s in by_gkey.get(k, []) if s['perfect'] and s['pam_ok'])
            if c != '' else c
            for c, k in zip(df['offtarget_identical'], df['_gkey'])]
        df['offtarget_detail'] = [
            offtarget_summarise(
                n, i, u + sum(1 for s in by_gkey.get(k, []) if s['pam_unverified']))
            if n != '' else ''
            for n, i, u, k in zip(df['offtarget_count'], df['offtarget_identical'],
                                  df['_ot_unver'], df['_gkey'])]

    if ot_status == 'screened' and offtarget_sidecar:
        import offtarget_remote
        offtarget_remote.write_sidecar(offtarget_sidecar, ot_region, ot_spacer,
                                       taxid, len(dna), pam, ot_cfg)
    df = df.drop(columns=[c for c in ('_gkey', '_spacer_pam', '_ot_unver') if c in df],
                 errors='ignore')

    # Short-arm check: only runs when protein_length is provided (same gate as
    # the coverage check above).  If any reagent has an arm shorter than half
    # the configured arm_length, the genomic sequence lacks enough flanking
    # sequence for HDR — typically because the user supplied a short isoform
    # transcript instead of a proper genomic region.
    if protein_length and protein_length > 0 and not df.empty:
        min_arm   = arm_length // 2
        min_left  = df['left_arm'].str.len().min()
        min_right = df['right_arm'].str.len().min()
        if min_left < min_arm or min_right < min_arm:
            raise ValueError(
                'Mismatch between DNA and protein sequences '
                '(shortest arm: left={} nt, right={} nt).\n'
                'Confirm genomic sequence covers (1) correct gene isoform and '
                '(2) contains sufficient flanking sequence to make full-length '
                'homology arms (>{} bp).'.format(
                    min_left, min_right, arm_length)
            )

    # Genomic sequence is always available here, so genotyping (screening) primers
    # are always designed: an external pair alone when insert_sequence is empty,
    # or external + 5'/3' junction pairs when a tag sequence is given.
    genotyping_df = _design_genotyping_rows(
        sites, dna, L, insert_sequence.upper(), internal_threshold,
        primer_opt_tm, product_opt_size, flank_min, flank_max,
        offtarget_hsps=ot_hsps if ot_status == 'screened' else None,
        offtarget_cfg=ot_cfg,
    )
    _report(reporter, '{} genotyping primer pairs designed'.format(len(genotyping_df)),
           stage='genotyping_primers')
    if ot_status == 'screened' and 'offtarget_amplicons' in genotyping_df:
        n_flagged = int((genotyping_df['offtarget_amplicons'] > 0).sum())
        if n_flagged:
            _report(reporter, '{} genotyping pair(s) may also amplify a duplicated '
                              'locus'.format(n_flagged),
                    stage='genotyping_primers', level='warning')

    return df, genotyping_df


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(genewise, genomic_fasta, output, protein_length=None, n_guides=5,
         arm_length=1000, pam='NGG', guide_length=20, cut_offset=3,
         insert_sequence='', internal_threshold=500, primer_opt_tm=60.0,
         product_opt_size=200, flank_min=50, flank_max=150, rs3=True,
         rs3_tracr='Hsu2013', offtarget=True, taxid='', email='',
         report=None, job_id_cb=None, resume_job_ids=None):
    """Entry point for in-process calls from task_runners."""
    result = design_reagents(
        genewise_out       = genewise,
        genomic_fasta      = genomic_fasta,
        protein_length     = protein_length,
        n_guides           = n_guides,
        arm_length         = arm_length,
        pam                = pam,
        guide_length       = guide_length,
        cut_offset         = cut_offset,
        insert_sequence    = insert_sequence,
        internal_threshold = internal_threshold,
        primer_opt_tm      = primer_opt_tm,
        product_opt_size   = product_opt_size,
        flank_min          = flank_min,
        flank_max          = flank_max,
        rs3                = rs3,
        rs3_tracr          = rs3_tracr,
        offtarget          = offtarget,
        taxid              = taxid,
        email              = email,
        # sidecar sits beside the reagents TSV so the UI can screen its own
        # on-demand primers without resubmitting a BLAST
        offtarget_sidecar  = str(Path(output).with_suffix('')) + '.offtarget.json',
        report             = report,
        job_id_cb          = job_id_cb,
        resume_job_ids     = resume_job_ids,
    )
    # a queued EBI job returns the sentinel instead of DataFrames; propagate it so
    # progress_server.py marks the task resumable rather than failed
    if isinstance(result, dict) and 'ebi_status' in result:
        return result
    df, genotyping_df = result
    df.to_csv(output, sep='\t', index=False)
    if genotyping_df is not None:
        genotyping_out = str(Path(output).with_suffix('')) + '.genotyping.tsv'
        genotyping_df.to_csv(genotyping_out, sep='\t', index=False)


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description=(
            'Design CRISPR knock-in reagents for every potential tag-insertion '
            'site in a protein, given Genewise alignment output.'
        )
    )
    # Required
    parser.add_argument('--genewise', required=True,
                        help='Genewise .out.txt result file')
    parser.add_argument('--genomic_fasta', required=True,
                        help='Genomic region FASTA (same sequence/orientation submitted to Genewise)')
    parser.add_argument('--output', required=True,
                        help='Output TSV file path')
    parser.add_argument('--protein_fasta',
                        help='Protein FASTA; used to compute protein length for CDS coverage check')
    # Optional guide / arm parameters
    parser.add_argument('--n_guides', type=int, default=5,
                        help='Max guide RNAs to report per insertion site (default: 5)')
    parser.add_argument('--arm_length', type=int, default=500,
                        help='Homology arm length in bp on each side (default: 500)')
    parser.add_argument('--PAM', default='NGG',
                        help='IUPAC PAM sequence (default: NGG)')
    parser.add_argument('--guide_length', type=int, default=20,
                        help='Guide spacer length in nt (default: 20)')
    parser.add_argument('--cut_offset', type=int, default=3,
                        help='Bases upstream of PAM where DSB occurs (SpCas9=3)')
    # Optional genotyping-primer parameters — genotyping primers are always
    # designed and written to a companion <output>.genotyping.tsv; setting
    # --insert_sequence additionally designs 5'/3' junction primer pairs.
    parser.add_argument('--insert_sequence', default='',
                        help='Insert/tag DNA sequence — when set, also designs 5\'/3\' '
                             'junction genotyping primers spanning the tag boundaries')
    parser.add_argument('--internal_threshold', type=int, default=500,
                        help='Unused; kept for CLI compatibility')
    parser.add_argument('--primer_opt_tm', type=float, default=60.0,
                        help='primer3 PRIMER_OPT_TM for genotyping primers (default: 60)')
    parser.add_argument('--product_opt_size', type=int, default=200,
                        help='primer3 PRIMER_PRODUCT_OPT_SIZE, the optimal amplicon length '
                             '(default: 200)')
    parser.add_argument('--flank_min', type=int, default=50,
                        help='Minimum distance (bp) of a genotyping primer from the insert site '
                             '(default: 50)')
    parser.add_argument('--flank_max', type=int, default=150,
                        help='Maximum distance (bp) of a genotyping primer from the insert site '
                             '(default: 150)')
    parser.add_argument('--no_rs3', action='store_true',
                        help='Skip Rule Set 3 on-target scoring. Scores are display-only and '
                             'never change which guides are chosen, so skipping only blanks '
                             'the rs3_score/rs3_percentile columns')
    parser.add_argument('--rs3_tracr', type=str, default='Hsu2013',
                        choices=['Hsu2013', 'Chen2013'],
                        help='tracrRNA scaffold assumed by the RS3 model (default: Hsu2013)')
    parser.add_argument('--no_offtarget', action='store_true',
                        help='Skip the off-target / primer-specificity BLAST screens. Results '
                             'are display-only and never change which guides are chosen, so '
                             'skipping only blanks the offtarget_* columns')
    parser.add_argument('--taxid', type=str, default='',
                        help='Species taxid scoping the off-target blastn search. Without it '
                             'the screens are skipped (a genome-wide search needs a species)')
    parser.add_argument('--email', type=str, default='',
                        help='E-mail address for EBI job submission, required by their REST '
                             'API; needed only when the off-target screens run')
    args, unknowns = parser.parse_known_args()

    protein_length = None
    if args.protein_fasta:
        recs = list(SeqIO.parse(args.protein_fasta, 'fasta'))
        if recs:
            protein_length = len(recs[0].seq)

    print('Running design_tag_reagents.py', file=sys.stderr)
    print('  genewise      : {}'.format(args.genewise), file=sys.stderr)
    print('  genomic_fasta : {}'.format(args.genomic_fasta), file=sys.stderr)
    print('  protein_length: {}'.format(protein_length), file=sys.stderr)
    print('  PAM           : {}'.format(args.PAM), file=sys.stderr)
    print('  arm_length    : {}'.format(args.arm_length), file=sys.stderr)
    print('  n_guides      : {}'.format(args.n_guides), file=sys.stderr)

    cli_result = design_reagents(
        genewise_out        = args.genewise,
        genomic_fasta       = args.genomic_fasta,
        protein_length      = protein_length,
        n_guides            = args.n_guides,
        arm_length          = args.arm_length,
        pam                 = args.PAM,
        guide_length        = args.guide_length,
        cut_offset          = args.cut_offset,
        insert_sequence     = args.insert_sequence,
        internal_threshold  = args.internal_threshold,
        primer_opt_tm       = args.primer_opt_tm,
        product_opt_size    = args.product_opt_size,
        flank_min           = args.flank_min,
        flank_max           = args.flank_max,
        rs3                 = not args.no_rs3,
        rs3_tracr           = args.rs3_tracr,
        offtarget           = not args.no_offtarget,
        taxid               = args.taxid,
        email               = args.email,
        offtarget_sidecar   = str(Path(args.output).with_suffix('')) + '.offtarget.json',
    )

    # standalone CLI has no resume machinery, so a queued EBI job is just an exit
    if isinstance(cli_result, dict) and 'ebi_status' in cli_result:
        print('EBI job {}: {}'.format(cli_result['ebi_status'], cli_result.get('detail', '')),
              file=sys.stderr)
        sys.exit(1)
    df, genotyping_df = cli_result

    df.to_csv(args.output, sep='\t', index=False)
    print('Wrote {} rows to {}'.format(len(df), args.output), file=sys.stderr)
    if genotyping_df is not None:
        genotyping_out = str(Path(args.output).with_suffix('')) + '.genotyping.tsv'
        genotyping_df.to_csv(genotyping_out, sep='\t', index=False)
        print('Wrote {} genotyping primer rows to {}'.format(len(genotyping_df), genotyping_out),
              file=sys.stderr)
