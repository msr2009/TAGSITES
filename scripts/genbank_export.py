"""
genbank_export.py

Build annotated GenBank (.gb) records for ApE / SnapGene from a reagents row.

Both records span the **entire genomic region** the pipeline searched, not just
the homology arms, so the files open in ApE as the locus in its full context:

  gDNA     the WT region, annotated with exon structure, the crRNA spacer and
           PAM, the Cas9 cut site, the insertion point, the homology arms and
           any genotyping primers.

  knockin  the same region with the repair product spliced in at the insertion
           point — the mutated arms replacing their WT counterparts and the tag
           between them — annotated with the arms, the inserted tag, the
           PAM-disrupting mutation(s) and any genotyping primers.

Coordinates in a reagents row (pam_fwd_start, cut_pos, insert_pos) are 0-based
in that region, so the gDNA record needs no coordinate translation at all.

Splicing the arms in preserves coordinates: an arm is anchored at the junction
and replaces exactly the WT bases it corresponds to, so the only shift in the
knock-in record is the insert itself — everything at or past insert_pos moves by
len(insert) and everything before it is unmoved.

Exon structure comes from the Genewise CDS intervals when they are supplied, and
falls back to the arms' own case encoding (uppercase exonic, lowercase intronic)
when only a reagents TSV is available.

Network-free and Shiny-free: everything takes plain strings / dict-likes.

Matt Rich, 2026
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from crispr_util import reverse_complement

# ApE renders misc_features in the color given by ApEinfo_fwdcolor/revcolor
COLORS = {
    "exon":      "#2e7d32",
    "spacer":    "#d32f2f",
    "pam":       "#f9a825",
    "cds":       "#1b5e20",
    "junction":  "#00838f",
    "insert":    "#e65100",
    "left_arm":  "#90caf9",
    "right_arm": "#a5d6a7",
    "mutation":  "#c2185b",
    "primer_f":  "#1565c0",
    "primer_r":  "#6a1b9a",
}

# genotyping amplicon types, in the order they should appear in a record
PRIMER_TYPES = ("external", "5p_junction", "3p_junction")


def _feat(start, end, label, ftype="misc_feature", strand=1, color=None, note=None):
    """One SeqFeature with an ApE-visible label and color; clamps to >=1 bp."""
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    if end <= start:
        end = start + 1
    quals = {"label": [label]}
    if color:
        quals["ApEinfo_fwdcolor"] = [color]
        quals["ApEinfo_revcolor"] = [color]
    if note:
        quals["note"] = [note]
    return SeqFeature(FeatureLocation(int(start), int(end), strand=strand),
                      type=ftype, qualifiers=quals)


def _clip(feats, length):
    """Drop features falling entirely outside [0, length); trim the rest."""
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    kept = []
    for f in feats:
        # a joined location (the spliced CDS) is built in range already, and
        # rebuilding it from .start/.end would collapse the join across introns
        if len(f.location.parts) > 1:
            kept.append(f)
            continue
        s, e = int(f.location.start), int(f.location.end)
        # a feature wholly off either end carries no information in this window
        if e <= 0 or s >= length:
            continue
        s, e = max(s, 0), min(e, length)
        if e <= s:
            continue
        kept.append(SeqFeature(FeatureLocation(s, e, strand=f.location.strand),
                               type=f.type, qualifiers=dict(f.qualifiers)))
    return kept


def exon_features(seq, offset=0, label="exon"):
    """Features for each run of uppercase (exonic) bases in a case-encoded arm."""
    feats = []
    start = None
    for i, ch in enumerate(seq):
        if ch.isupper() and start is None:
            start = i
        elif not ch.isupper() and start is not None:
            feats.append(_feat(offset + start, offset + i, label,
                               ftype="exon", color=COLORS["exon"]))
            start = None
    # a run reaching the end of the string closes here rather than in the loop
    if start is not None:
        feats.append(_feat(offset + start, offset + len(seq), label,
                           ftype="exon", color=COLORS["exon"]))
    return feats


def guide_features(row, region_start, guide_length=None, cut_offset=3):
    """Spacer, PAM and cut-site features for a guide, in record-local coords.

    Mirrors crispr_util.find_guides: on '+' the spacer sits immediately 5' of
    the PAM, on '-' it sits immediately 3' of it in forward coordinates.
    """
    spacer  = str(row["spacer"]).upper()
    pam     = str(row["pam_seq"]).upper()
    strand  = str(row["guide_strand"])
    glen    = int(guide_length) if guide_length else len(spacer)
    pam_s   = int(row["pam_fwd_start"]) - region_start
    sign    = 1 if strand == "+" else -1

    if strand == "+":
        sp_s, sp_e = pam_s - glen, pam_s
    else:
        sp_s, sp_e = pam_s + len(pam), pam_s + len(pam) + glen

    return [
        _feat(sp_s, sp_e, "crRNA {} ({})".format(spacer, strand),
              ftype="primer_bind", strand=sign, color=COLORS["spacer"],
              note="guide spacer, {} strand; cut {} bp from PAM".format(
                  strand, cut_offset)),
        _feat(pam_s, pam_s + len(pam), "PAM {}".format(pam),
              strand=sign, color=COLORS["pam"]),
    ]


def find_oligo(seq, oligo):
    """Locate an oligo in seq on either strand; returns (start, end, strand) or None.

    Used as the fallback when a primer has no recorded region span, and as the
    only route for the knock-in record (whose coordinates are not genomic).
    Returns None when absent or ambiguous — a primer matching twice in a 2 kb
    window is itself a problem and should not be silently annotated once.
    """
    if not oligo:
        return None
    hay = seq.upper()
    for needle, strand in ((oligo.upper(), 1), (reverse_complement(oligo).upper(), -1)):
        first = hay.find(needle)
        if first < 0:
            continue
        if hay.find(needle, first + 1) >= 0:
            return None
        return (first, first + len(needle), strand)
    return None


def primer_features(primers, seq, types=PRIMER_TYPES, region_start=None,
                    shift_at=None, shift_by=0):
    """Features for the requested genotyping amplicon types.

    primers      {amplicon_type: {fwd_seq, rev_seq, fwd_region_span, ...}}
    region_start record-local offset of region coordinate 0, or None to always
                 locate primers by sequence search instead of by span. Spans are
                 preferred when given because a primer overlapping a
                 PAM-disrupting mutation will not string-match the locus.
    shift_at     region coordinate at or after which spans move by shift_by (the
                 insertion point in a knock-in record). A primer with no span at
                 all — a junction pair's internal primer, which lies inside the
                 tag and so has no genomic coordinates — falls back to search.
    """
    feats = []
    for atype in types:
        pair = (primers or {}).get(atype)
        if not pair:
            continue
        for side, key, color in (("F", "fwd", COLORS["primer_f"]),
                                 ("R", "rev", COLORS["primer_r"])):
            oligo = str(pair.get("{}_seq".format(key), "") or "")
            if not oligo:
                continue
            span = pair.get("{}_region_span".format(key))
            if span and region_start is not None:
                lo, hi = span
                if shift_at is not None and lo >= shift_at:
                    lo, hi = lo + shift_by, hi + shift_by
                start, end = lo - region_start, hi - region_start
                strand = 1 if side == "F" else -1
            else:
                hit = find_oligo(seq, oligo)
                if hit is None:
                    continue
                start, end, strand = hit
            tm = pair.get("{}_tm".format(key))
            note = "Tm {:.1f}C".format(float(tm)) if tm is not None else None
            feats.append(_feat(start, end, "{}_{}".format(atype, side),
                               ftype="primer_bind", strand=strand, color=color,
                               note=note))
    return feats


def mutation_features(arm, arm_wt, offset, side):
    """Per-base features wherever a mutated arm differs from its WT counterpart.

    Compares from the junction outward, since the arms are anchored there and
    may differ in length when a mutation pushed one arm past arm_length.
    """
    feats = []
    n = min(len(arm), len(arm_wt))
    if n == 0:
        return feats
    # left arms end at the junction, right arms start there — align on that end
    a, w = (arm[-n:], arm_wt[-n:]) if side == "left" else (arm[:n], arm_wt[:n])
    base = offset + (len(arm) - n if side == "left" else 0)
    for i, (x, y) in enumerate(zip(a, w)):
        if x.upper() != y.upper():
            feats.append(_feat(base + i, base + i + 1,
                               "{}>{}".format(y.upper(), x.upper()),
                               ftype="variation", color=COLORS["mutation"],
                               note="recut-blocking edit"))
    return feats


def site_label(row):
    """Compact site/guide label, e.g. K45+3 — matches reagents_server._guide_label."""
    return "{}{}{}{}".format(
        str(row["amino_acid"]), int(row["residue_index"]),
        str(row["guide_strand"]), int(row["distance"]),
    )


def _record(seq, rec_id, description, features):
    """Assemble a SeqRecord with a GenBank-legal LOCUS name and clipped features."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(
        Seq(str(seq).upper()),
        id=rec_id,
        name=rec_id.replace(" ", "_")[:16],
        description=description,
        annotations={"molecule_type": "DNA"},
    )
    rec.features = sorted(_clip(features, len(seq)),
                          key=lambda f: (int(f.location.start), int(f.location.end)))
    return rec


def cds_exon_features(exons, shift_at=None, shift_by=0, label="exon"):
    """One feature per Genewise CDS exon, shifted past an insert if given.

    An exon straddling shift_at is split, so the tag is drawn between the two
    halves rather than inside one of them.
    """
    feats = []
    for i, exon in enumerate(exons or []):
        # shift one exon at a time so a split keeps the original exon's number
        parts = shift_exons([exon], shift_at, shift_by, absorb=False)
        for k, (s, e) in enumerate(parts):
            name = "{} {}".format(label, i + 1)
            if len(parts) > 1:
                name += " ({}' part)".format("5" if k == 0 else "3")
            feats.append(_feat(s, e, name, ftype="exon", color=COLORS["exon"]))
    return feats


def shift_exons(exons, shift_at=None, shift_by=0, absorb=False):
    """Genewise CDS intervals as half-open spans, shifted past an insert.

    exons    iterable of (start, stop) 0-based INCLUSIVE, as parse_genewise gives
    shift_at region coordinate at or after which coordinates move, or None
    shift_by how far they move (len(insert))
    absorb   True to swallow the inserted bases into the exon that straddles
             shift_at, making the tag part of the coding sequence; False to
             split that exon around them
    """
    out = []
    for s, e in (exons or []):
        s, e = int(s), int(e) + 1
        if shift_at is None:
            out.append((s, e))
        elif absorb:
            # the tag belongs to whichever exon it lands in OR abuts, so that an
            # insert at an exon boundary (the common N-/C-terminal tag) is still
            # part of the coding sequence rather than falling into an intron
            if e < shift_at:
                out.append((s, e))
            elif s > shift_at:
                out.append((s + shift_by, e + shift_by))
            else:
                out.append((s, e + shift_by))
        elif e <= shift_at:
            out.append((s, e))
        elif s >= shift_at:
            out.append((s + shift_by, e + shift_by))
        else:
            out.append((s, shift_at))
            out.append((shift_at + shift_by, e + shift_by))
    return out


def cds_feature(seq, spans, strand=1, label="CDS", note=None):
    """One joined CDS feature over spans, carrying the conceptual translation.

    Spans are half-open and in ascending coordinate order. ApE and SnapGene both
    render /translation, which is what makes the reading frame across the tag
    junction checkable by eye. Returns None when there is nothing to translate.
    """
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

    spans = [(int(s), int(e)) for s, e in spans if int(e) > int(s)]
    if not spans:
        return None
    parts = [FeatureLocation(s, e, strand=strand) for s, e in spans]
    if strand == -1:
        parts = parts[::-1]
    location = parts[0] if len(parts) == 1 else CompoundLocation(parts)

    coding = "".join(str(seq)[s:e] for s, e in spans).upper()
    if strand == -1:
        coding = reverse_complement(coding)
    # a Genewise CDS need not be a whole number of codons, so trim rather than
    # let Bio.Seq.translate warn and pad
    aa = str(Seq(coding[:len(coding) - len(coding) % 3]).translate())

    quals = {
        "label":            [label],
        "translation":      [aa],
        "codon_start":      ["1"],
        "ApEinfo_fwdcolor": [COLORS["cds"]],
        "ApEinfo_revcolor": [COLORS["cds"]],
    }
    if note:
        quals["note"] = [note]
    return SeqFeature(location, type="CDS", qualifiers=quals)


def _arm_features(left_start, left_len, right_start, right_len):
    """Spanning features for the two homology arms."""
    return [
        _feat(left_start, left_start + left_len, "5' homology arm",
              color=COLORS["left_arm"], note="{} bp".format(left_len)),
        _feat(right_start, right_start + right_len, "3' homology arm",
              color=COLORS["right_arm"], note="{} bp".format(right_len)),
    ]


def build_gdna_record(row, region_seq, left_wt=None, right_wt=None, exons=None,
                      primers=None, types=PRIMER_TYPES, run_name="",
                      guide_length=None, cut_offset=3, cds_strand=1):
    """The whole WT genomic region, annotated for genotyping and cutting.

    region_seq is the genomic FASTA the pipeline searched — row coordinates are
    already in its frame, so nothing is translated. left_wt/right_wt only set
    the extent of the homology-arm features; exons, when given, come from
    parse_genewise and take precedence over the arms' case encoding.
    """
    seq  = str(region_seq)
    ip   = int(row["insert_pos"])
    site = site_label(row)

    if exons:
        feats = cds_exon_features(exons)
        cds = cds_feature(seq, shift_exons(exons), cds_strand,
                          label="CDS (WT)")
        if cds:
            feats.append(cds)
    else:
        # no Genewise output on hand: recover what the arms' case encoding shows,
        # which covers only the arm window rather than the whole region
        feats  = exon_features(str(left_wt or ""), ip - len(left_wt or ""))
        feats += exon_features(str(right_wt or ""), ip)
    if left_wt and right_wt:
        feats += _arm_features(ip - len(left_wt), len(left_wt), ip, len(right_wt))
    feats += guide_features(row, 0, guide_length, cut_offset)
    feats.append(_feat(ip - 1, ip + 1, "insertion site",
                       color=COLORS["junction"],
                       note="tag inserted between these two bases"))
    feats += primer_features(primers, seq, types, region_start=0)

    rid = "{}_{}_gDNA".format(run_name, site) if run_name else "{}_gDNA".format(site)
    return _record(seq, rid,
                   "WT genomic region ({} bp); tag site {} at region position {}".format(
                       len(seq), site, ip),
                   feats)


def splice_knockin(region_seq, insert_pos, left, right, insert):
    """Region with the mutated arms and the insert spliced in at insert_pos.

    The arms replace exactly the WT bases they correspond to, so every
    coordinate below insert_pos is unmoved and every one at or above it shifts
    by len(insert).
    """
    region_seq, left, right = str(region_seq), str(left), str(right)
    insert = str(insert or "")
    lo = max(insert_pos - len(left), 0)
    hi = min(insert_pos + len(right), len(region_seq))
    return region_seq[:lo] + left + insert + right + region_seq[hi:]


def build_knockin_record(row, region_seq, left, right, insert, left_wt=None,
                         right_wt=None, exons=None, primers=None,
                         types=PRIMER_TYPES, insert_name="tag", run_name="",
                         guide_length=None, cut_offset=3, cds_strand=1):
    """The whole genomic region carrying the repair product, annotated.

    left/right are the mutated (PAM-disrupted) arms as they will be ordered.
    left_wt/right_wt, when given, locate the recut-blocking edits.
    """
    left, right, insert = str(left), str(right), str(insert or "")
    ip   = int(row["insert_pos"])
    seq  = splice_knockin(region_seq, ip, left, right, insert)
    L    = len(insert)
    ins1 = ip + L
    site = site_label(row)

    if exons:
        feats = cds_exon_features(exons, shift_at=ip, shift_by=L)
        # absorb=True puts the tag inside the CDS, so the translation shows it
        # in frame with the protein — the whole point of a translated knock-in
        cds = cds_feature(seq, shift_exons(exons, ip, L, absorb=True), cds_strand,
                          label="CDS ({} knock-in)".format(insert_name),
                          note="includes the inserted {} bp {}".format(L, insert_name))
        if cds:
            feats.append(cds)
    else:
        feats = exon_features(left, ip - len(left)) + exon_features(right, ins1)
    feats += _arm_features(ip - len(left), len(left), ins1, len(right))
    if insert:
        feats.append(_feat(ip, ins1, insert_name, color=COLORS["insert"],
                           note="inserted sequence, {} bp".format(L)))
    if left_wt:
        feats += mutation_features(left, str(left_wt), ip - len(left), "left")
    if right_wt:
        feats += mutation_features(right, str(right_wt), ins1, "right")
    # the guide is annotated on the edited sequence only when its PAM survived
    # intact; a disrupted PAM no longer matches and is left off rather than
    # drawn somewhere misleading
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    for f in guide_features(row, 0, guide_length, cut_offset):
        s, e, strand = int(f.location.start), int(f.location.end), f.location.strand
        if s >= ip:
            s, e = s + L, e + L
        elif e > ip:
            continue          # the insert lands inside this feature
        if f.qualifiers["label"][0].startswith("PAM"):
            sub = seq[s:e].upper()
            if (sub if strand == 1 else reverse_complement(sub)) != str(row["pam_seq"]).upper():
                continue
        feats.append(SeqFeature(FeatureLocation(s, e, strand=strand),
                                type=f.type, qualifiers=dict(f.qualifiers)))
    feats += primer_features(primers, seq, types, region_start=0,
                             shift_at=ip, shift_by=L)

    desc = "Knock-in at {}: {} bp region with {} bp {} inserted ({} bp arms)".format(
        site, len(seq), L, insert_name, len(left))
    rid = "{}_{}_knockin".format(run_name, site) if run_name else "{}_knockin".format(site)
    return _record(seq, rid, desc, feats)


def write_record(record, path):
    """Write one SeqRecord as GenBank."""
    from Bio import SeqIO

    with open(path, "w") as fout:
        SeqIO.write(record, fout, "genbank")


# ── CLI ───────────────────────────────────────────────────────────────────────

def load_exons(genewise_out):
    """CDS exon intervals (0-based inclusive) from a Genewise .out.txt, or None.

    run_genewise writes the genomic FASTA already flipped to the winning
    orientation, so these are effectively always on '+'; load_cds_strand reads
    the strand rather than assuming it.
    """
    if not genewise_out or not Path(genewise_out).exists():
        return None
    from parse_genewise import parse_genewise

    try:
        cds = parse_genewise(genewise_out)
    except Exception:
        return None
    return [(int(r["start"]), int(r["stop"])) for _, r in cds.iterrows()]


def load_cds_strand(genewise_out):
    """+1 / -1 for the CDS in a Genewise .out.txt; +1 when unknown or mixed."""
    if not genewise_out or not Path(genewise_out).exists():
        return 1
    from parse_genewise import parse_genewise

    try:
        strands = set(parse_genewise(genewise_out)["strand"])
    except Exception:
        return 1
    return -1 if strands == {"-"} else 1


def load_region(genomic_fasta):
    """First record's sequence from a genomic FASTA, or None."""
    if not genomic_fasta or not Path(genomic_fasta).exists():
        return None
    from Bio import SeqIO

    recs = list(SeqIO.parse(genomic_fasta, "fasta"))
    return str(recs[0].seq) if recs else None


def main(reagents_tsv, outdir, residue=None, guide_id=None, insert_sequence="",
         insert_name="tag", types=PRIMER_TYPES, genotyping_tsv=None,
         run_name="", guide_length=None, cut_offset=3, genomic_fasta=None,
         genewise_out=None, report=None):
    """Write gDNA + knock-in GenBank files for selected rows of a reagents TSV.

    residue/guide_id filter which rows are exported; with neither, every row is
    exported, which on a whole region is thousands of files — so a residue is
    effectively required in practice.
    """
    import pandas as pd

    df = pd.read_csv(reagents_tsv, sep="\t")
    if residue is not None:
        df = df[df["residue_index"] == int(residue)]
    if guide_id is not None and "guide_id" in df.columns:
        df = df[df["guide_id"].astype(str) == str(guide_id)]
    if df.empty:
        raise SystemExit("No reagents rows match residue={} guide_id={}".format(
            residue, guide_id))

    # genotyping primers live in the companion TSV, keyed by residue + type
    primers_by_rid = {}
    if genotyping_tsv and Path(genotyping_tsv).exists():
        gdf = pd.read_csv(genotyping_tsv, sep="\t")
        for rid, grp in gdf.groupby("residue_index"):
            primers_by_rid[int(rid)] = {
                str(r["amplicon_type"]): {
                    "fwd_seq": str(r["fwd_seq"]), "rev_seq": str(r["rev_seq"]),
                    "fwd_tm": r.get("fwd_tm"), "rev_tm": r.get("rev_tm"),
                }
                for _, r in grp.iterrows()
            }

    region = load_region(genomic_fasta)
    exons  = load_exons(genewise_out)
    strand = load_cds_strand(genewise_out)
    if region is None:
        raise SystemExit(
            "--genomic_fasta is required: the records span the whole genomic region.")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for _, row in df.iterrows():
        rid     = int(row["residue_index"])
        primers = primers_by_rid.get(rid)
        left_wt  = str(row.get("left_arm_wt") or row["left_arm"])
        right_wt = str(row.get("right_arm_wt") or row["right_arm"])

        gdna = build_gdna_record(row, region, left_wt, right_wt, exons, primers,
                                 types, run_name=run_name,
                                 guide_length=guide_length, cut_offset=cut_offset,
                                 cds_strand=strand)
        ki = build_knockin_record(row, region, str(row["left_arm"]),
                                  str(row["right_arm"]), insert_sequence,
                                  left_wt, right_wt, exons, primers, types,
                                  insert_name, run_name=run_name,
                                  guide_length=guide_length, cut_offset=cut_offset,
                                  cds_strand=strand)
        for rec in (gdna, ki):
            path = outdir / "{}.gb".format(rec.id)
            write_record(rec, path)
            written.append(str(path))
            if report:
                report("wrote {}".format(path))
    return written


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Write annotated GenBank (ApE) files for TAGSITES reagents.")
    p.add_argument("--reagents", required=True, help="a *_reagents.tsv from design_tag_reagents.py")
    p.add_argument("--outdir", default=".", help="directory for the .gb files")
    p.add_argument("--residue", default=None, help="restrict to this residue_index")
    p.add_argument("--guide_id", default=None, help="restrict to this guide_id")
    p.add_argument("--insert_sequence", default="", help="tag DNA sequence for the knock-in record")
    p.add_argument("--insert_name", default="tag", help="label for the inserted sequence")
    p.add_argument("--genotyping", default=None,
                   help="companion *.genotyping.tsv, to annotate primers")
    p.add_argument("--primer_types", default=",".join(PRIMER_TYPES),
                   help="comma-separated amplicon types to annotate")
    p.add_argument("--run_name", default="", help="prefix for record ids and filenames")
    p.add_argument("--guide_length", type=int, default=None)
    p.add_argument("--cut_offset", type=int, default=3)
    p.add_argument("--genomic_fasta", required=True,
                   help="the genomic region FASTA the reagents were designed against "
                        "(*_genewise.genewise_genomic.fa)")
    p.add_argument("--genewise", dest="genewise_out", default=None,
                   help="the Genewise *.genewise.out.txt, for exon annotation")
    a = p.parse_args()

    out = main(a.reagents, a.outdir, a.residue, a.guide_id, a.insert_sequence,
               a.insert_name, tuple(t for t in a.primer_types.split(",") if t),
               a.genotyping, a.run_name, a.guide_length, a.cut_offset,
               a.genomic_fasta, a.genewise_out)
    print("\n".join(out))
