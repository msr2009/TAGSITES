"""phastcons_track.py — per-residue phastCons conservation from the run's own gene model.

No second gene model is built: the codons are those of the model the reagents use
(gene_model.py's winner, read back from {outprefix}.genewise.out.txt + .genewise_genomic.fa),
renumbered to the input protein through cds_check exactly as design_tag_reagents does.
UCSC is used only to (a) place that region on the assembly and (b) read the scores:

  placement   region.json (UCSC or bulk GFF3 models carry chromosome coordinates), else an
              exact match of the region (or its reverse complement) inside the gene's locus,
              found from the FASTA header's coordinates, the UniProt alignment, or the gene name
  agreement   the codon bases read back from the assembly at the mapped positions must
              translate to the same protein as the model itself; otherwise no track is
              written (the scores would belong to other codons)
  scores      UCSC's phastCons bigWig read over HTTP per CDS window (pyBigWig; the local
              copy in bulk mode, so both give identical values); without pyBigWig, the
              getData/track REST values, which UCSC byte-quantizes (~1/127 steps). A
              residue's score is the mean of its three bases, -1000 when none is scored

Output ({run}.reagents.phastcons.tsv, written beside model_check.json): '#' header lines
(genome, track, placement, gene mean), then 'pos<TAB>score' for every input-protein residue
in order, the layout utils/results.py reads for continuous tracks. Assemblies without a
phastCons track get no file. Fails soft: main() logs and returns None on any problem.

Matt Rich, 2026
"""

import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import cds_check
from crispr_util import reverse_complement
from parse_genewise import codon_positions, parse_genewise
from progress import report as _report, resolve_reporter

# narrow-clade alignment per assembly (closer orthologs, more informative); --track overrides
PHASTCONS_TRACKS = {
    "ce11": "phastCons26way", "hg38": "phastCons100way", "hg19": "phastCons100way",
    "mm39": "phastCons35way", "mm10": "phastCons60way", "dm6": "phastCons27way",
    "sacCer3": "phastCons7way", "galGal6": "phastCons77way",
}
LOCUS_PAD = 5000  # bases searched beyond the region's own length around a locus
WINDOW_GAP = 50000  # codons further apart than this are fetched in separate windows
BIGWIG_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/{track}/{genome}.{track}.bw"


def _complement(base):
    """Complement of one base."""
    return {"A": "T", "C": "G", "G": "C", "T": "A"}.get(base, "N")


def _windows(positions, gap=WINDOW_GAP):
    """Sorted chromosome positions -> [(start0, end)] windows split at gaps > gap."""
    pos = sorted(positions)
    out = [[pos[0], pos[0] + 1]]
    for p in pos[1:]:
        if p - out[-1][1] > gap:
            out.append([p, p + 1])
        else:
            out[-1][1] = p + 1
    return [tuple(w) for w in out]


class _Remote:
    """UCSC REST access for one assembly (sequence, scores, loci), paced like ucsc_rest."""

    def __init__(self, genome):
        import ensembl_rest
        import http_retry

        self.genome = genome
        self.score_source = ""
        self.cfg = ensembl_rest.load_config()
        self.deadline = http_retry.deadline_from(ensembl_rest.timeouts()["total_timeout_s"])

    def sequence(self, chrom, start0, end):
        """Assembly sequence (+ strand, upper case) for [start0, end)."""
        import offtarget_blat
        import ucsc_rest

        ucsc_rest._pause(self.deadline)
        return offtarget_blat.fetch_sequence(self.genome, chrom, start0, end,
                                             deadline=self.deadline).upper()

    def scores(self, track, chrom, start0, end):
        """{0-based position: phastCons} over [start0, end).

        Read from UCSC's bigWig over HTTP when pyBigWig is installed (full precision, the
        same file bulk mode downloads); otherwise from getData/track, whose wiggle values
        are byte-quantized (steps of ~1/127).
        """
        import math

        import ucsc_rest

        try:
            import pyBigWig

            bw = pyBigWig.open(BIGWIG_URL.format(genome=self.genome, track=track))
            try:
                vals = bw.values(chrom, start0, end)
            finally:
                bw.close()
            self.score_source = "bigWig"
            return {start0 + i: v for i, v in enumerate(vals) if not math.isnan(v)}
        except Exception:
            self.score_source = "REST (quantized)"
        ucsc_rest._pause(self.deadline)
        vals = {}
        for it in ucsc_rest.track_items(self.genome, track, chrom, start0, end, self.cfg,
                                        deadline=self.deadline):
            for p in range(int(it["start"]), int(it["end"])):
                vals[p] = it["value"]
        return vals

    def aliases(self):
        """{alias: UCSC chrom}."""
        import ucsc_rest

        ucsc_rest._pause(self.deadline)
        return ucsc_rest.chrom_aliases(self.genome, self.cfg, deadline=self.deadline)

    def loci(self, accession, gene_name):
        """Candidate (chrom, start1, end1) loci for the gene: UniProt alignment, then symbol."""
        import ucsc_rest

        spans = []
        if accession:
            ucsc_rest._pause(self.deadline)
            spans = [(c, s, e) for _, c, s, e, _ in
                     ucsc_rest.search_accession(self.genome, accession, self.cfg,
                                                deadline=self.deadline)]
        if not spans and gene_name:
            ucsc_rest._pause(self.deadline)
            spans = ucsc_rest.search_gene(self.genome, gene_name, self.cfg,
                                          deadline=self.deadline)[1]
        return spans


class _Local:
    """Bulk-mode access: local genome (Ensembl names) and the local phastCons bigWig."""

    def __init__(self, genome):
        self.genome = genome
        self.score_source = "local bigWig"

    def sequence(self, chrom, start0, end):
        """Local genome sequence for a UCSC-named range."""
        import genome_regions
        import ucsc_local

        return genome_regions.extract_sequence(ucsc_local.local_chrom(chrom), start0 + 1, end,
                                               "+").upper()

    def scores(self, track, chrom, start0, end):
        """{0-based position: phastCons} from the local bigWig."""
        import ucsc_local

        return ucsc_local.phastcons_values(chrom, start0, end)

    def aliases(self):
        """No alias table offline; local names map by ucsc_local's rule."""
        return {}

    def loci(self, accession, gene_name):
        """Loci from the local UniProt index."""
        import ucsc_local

        return [(it["chrom"], int(it["chromStart"]) + 1, int(it["chromEnd"]))
                for it in ucsc_local.lookup(accession)] if accession else []


def _ucsc_chrom(name, aliases):
    """A chromosome name (UCSC, Ensembl or RefSeq) as the UCSC name."""
    if name in aliases:
        return aliases[name]
    if name.startswith("chr"):
        return name
    import ucsc_local

    # inverse of the bulk alias rule (I -> chrI, MtDNA -> chrM)
    inverse = {v: k for k, v in ucsc_local.settings()[2].items()}
    return inverse.get(name, "chr" + name)


def _from_region_json(path, n, aliases):
    """(chrom, map index -> chromosome position, strand) from a region.json, or None."""
    try:
        with open(path) as fh:
            r = json.load(fh)
    except (OSError, ValueError):
        return None
    plus = r["strand"] == "+"
    lo1 = int(r["start"]) - int(r["flank_5p"] if plus else r["flank_3p"])
    hi1 = int(r["stop"]) + int(r["flank_3p"] if plus else r["flank_5p"])
    if hi1 - lo1 + 1 != n:
        return None  # the region on disk is not the one region.json describes
    chrom = _ucsc_chrom(r["chrom"], aliases)
    if plus:
        return chrom, (lambda i: lo1 - 1 + i), "+"
    return chrom, (lambda i: hi1 - 1 - i), "-"


def _header_locus(genomic_fa, aliases):
    """(chrom, start1, end1) from an auto-fetch header >src:assembly:chrom:lo:hi:strand."""
    try:
        with open(genomic_fa) as fh:
            parts = fh.readline().lstrip(">").strip().split(":")
        if len(parts) >= 5:
            return _ucsc_chrom(parts[2], aliases), int(parts[3]), int(parts[4])
    except (OSError, ValueError):
        pass
    return None


def _match_region(src, dna, loci):
    """Place dna inside one of the loci by exact match; (chrom, mapper, strand) or None."""
    n = len(dna)
    rc = reverse_complement(dna)
    for chrom, start1, end1 in loci:
        s0 = max(0, start1 - 1 - n - LOCUS_PAD)
        seq = src.sequence(chrom, s0, end1 + n + LOCUS_PAD)
        fwd, rev = seq.count(dna), seq.count(rc)
        # exactly one occurrence, else the placement is ambiguous
        if fwd + rev != 1:
            continue
        if fwd:
            o = s0 + seq.index(dna)
            return chrom, (lambda i, o=o: o + i), "+"
        o = s0 + seq.index(rc)
        return chrom, (lambda i, o=o: o + n - 1 - i), "-"
    return None


def _gene_name(accession):
    """UniProt primary gene name for an accession, '' if unavailable."""
    if not accession:
        return ""
    try:
        from uniprot_api import fetch_entry

        genes = fetch_entry(accession.split("-")[0]).get("genes") or []
        return (genes[0].get("geneName") or {}).get("value", "") if genes else ""
    except Exception:
        return ""


def build_track(genewise_out, genomic_fa, protein_fa, taxid="", region_json=None, track="",
                backend="remote", accession="", report=None):
    """Per-residue phastCons for the input protein; returns a dict, raises on any failure.

    {"genome", "track", "placement", "scores": [score or None per query residue],
     "gene_mean"}.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from uniprot_api import resolve_accession

    reporter = resolve_reporter(report)
    dna = str(next(SeqIO.parse(genomic_fa, "fasta")).seq).upper()
    cds_df = parse_genewise(genewise_out)
    query = cds_check.read_protein(protein_fa).rstrip("*")
    model_prot = cds_check.translate_model(cds_df, dna)
    check = cds_check.compare_to_query(model_prot, query)
    if check["status"] == "fail":
        raise ValueError("gene model fails cds_check; no reagents, so no phastCons track")

    if backend == "local":
        import ucsc_local

        genome = ucsc_local.settings()[0]
        track = track or ucsc_local.settings()[1]
        src = _Local(genome)
    else:
        import offtarget_blat
        import offtarget_screen

        genome = offtarget_blat.assembly_for_taxid(taxid, offtarget_screen.load_config()) \
            if taxid else ""
        if not genome:
            raise ValueError(f"no UCSC assembly for taxid {taxid or '(none)'}")
        track = track or PHASTCONS_TRACKS.get(genome, "")
        src = _Remote(genome)
    if not track:
        raise ValueError(f"no phastCons track known for {genome}")

    # placement: region.json first, else an exact match inside the gene's locus
    aliases = {}
    placed, method = None, ""
    if region_json and os.path.exists(region_json):
        if backend != "local":
            aliases = src.aliases()
        placed, method = _from_region_json(region_json, len(dna), aliases), "region.json"
    if placed is None:
        if backend != "local" and not aliases:
            aliases = src.aliases()
        loci = []
        hdr = _header_locus(genomic_fa, aliases)
        if hdr:
            loci.append(hdr)
        acc = resolve_accession(protein_fa, accession, taxid)[0]
        loci += src.loci(acc, "" if loci else _gene_name(acc))
        if not loci:
            loci = src.loci("", _gene_name(acc))
        placed, method = _match_region(src, dna, loci), "exact region match"
    if placed is None:
        raise ValueError(f"could not place the genomic region on {genome}")
    chrom, to_chrom, strand = placed

    # the model's codons, keyed by input-protein residue (renumbered as design_tag_reagents does)
    rmap = {int(k): v for k, v in (check.get("residue_map") or {}).items()}
    codons = {}
    for i, c in enumerate(codon_positions(cds_df), start=1):
        q = rmap.get(i, i) if rmap else i
        if q is not None:
            codons[q] = [to_chrom(p) for p in c]
    allp = [p for c in codons.values() for p in c]
    windows = _windows(allp)

    # agreement: the assembly's bases at those positions must translate like the model
    genome_seq = {}
    for w0, w1 in windows:
        seq = src.sequence(chrom, w0, w1)
        genome_seq.update({w0 + i: b for i, b in enumerate(seq)})
    cds = "".join(genome_seq.get(to_chrom(p), "N") if strand == "+"
                  else _complement(genome_seq.get(to_chrom(p), "N"))
                  for c in codon_positions(cds_df) for p in c)
    readback = str(Seq(cds).translate()).rstrip("*").replace("*", "U")
    if readback != model_prot:
        raise ValueError(f"{genome} codons at the mapped positions do not translate to the gene "
                         "model's protein; refusing to mix the two")

    vals = {}
    for w0, w1 in windows:
        vals.update(src.scores(track, chrom, w0, w1))
    scores = []
    for q in range(1, len(query) + 1):
        got = [vals[p] for p in codons.get(q, []) if p in vals]
        scores.append(sum(got) / len(got) if got else None)
    present = [s for s in scores if s is not None]
    gene_mean = sum(present) / len(present) if present else None
    _report(reporter, f"phastCons {genome} {track}: {len(present)}/{len(scores)} residues "
                      f"scored, gene mean {gene_mean:.3f}" if present else
                      f"phastCons {genome} {track}: no residue scored", stage="phastcons")
    return {"genome": genome, "track": track, "placement": method, "chrom": chrom,
            "strand": strand, "scores": scores, "gene_mean": gene_mean,
            "score_source": src.score_source}


def write_track(result, output):
    """Write the per-residue TSV (one row per input-protein residue, -1000 = unscored)."""
    with open(output, "w") as fh:
        fh.write(f"# genome\t{result['genome']}\n# track\t{result['track']}\n"
                 f"# placement\t{result['placement']}\n"
                 f"# scores\t{result.get('score_source', '')}\n")
        mean = result["gene_mean"]
        fh.write(f"# gene_mean\t{'' if mean is None else round(mean, 4)}\n")
        for pos, s in enumerate(result["scores"], start=1):
            fh.write(f"{pos}\t{-1000 if s is None else round(s, 4)}\n")


def main(genewise_out, genomic_fa, protein_fa, output, taxid="", region_json=None, track="",
         backend="", accession="", report=None):
    """Build and write the track; returns the output path, or None (logged) on any failure."""
    reporter = resolve_reporter(report)
    if not backend:
        import providers

        backend = "local" if providers.region_from_annotation() else "remote"
    if region_json is None:
        region_json = genewise_out.replace(".genewise.out.txt", ".region.json")
    if not protein_fa or not os.path.exists(protein_fa):
        _report(reporter, "no phastCons track: no protein FASTA given", stage="phastcons",
                level="warning")
        return None
    try:
        result = build_track(genewise_out, genomic_fa, protein_fa, taxid, region_json, track,
                             backend, accession, report=report)
    except Exception as e:
        _report(reporter, f"no phastCons track: {e}", stage="phastcons", level="warning")
        try:
            os.remove(output)  # never leave a previous run's track behind
        except FileNotFoundError:
            pass
        return None
    write_track(result, output)
    return output


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Per-residue phastCons for a run's gene model.")
    parser.add_argument("--genewise", required=True, help="{run}_genewise.genewise.out.txt")
    parser.add_argument("--genomic_fasta", required=True, help="{run}_genewise.genewise_genomic.fa")
    parser.add_argument("--protein_fasta", "--input_file", dest="protein_fasta", required=True)
    parser.add_argument("--output", required=True, help="output TSV")
    parser.add_argument("--taxid", default="")
    parser.add_argument("--accession", "--uniprot_accession", dest="accession", default="")
    parser.add_argument("--region_json", default=None,
                        help="default: beside --genewise ({outprefix}.region.json)")
    parser.add_argument("--track", default="", help="phastCons track (default: narrow clade)")
    parser.add_argument("--backend", choices=("remote", "local"), default="")
    args, _ = parser.parse_known_args()
    out = main(args.genewise, args.genomic_fasta, args.protein_fasta, args.output, args.taxid,
               args.region_json, args.track, args.backend, args.accession)
    sys.exit(0 if out else 1)
