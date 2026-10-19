"""
gnomad_constraint.py

PROTOTYPE (see docs/GNOMAD_CONSTRAINT_SPEC.md): project human population missense
constraint from gnomAD onto a C. elegans protein through its Alliance human orthologs.

Per worm residue, two human-derived signals are reported:
  - MTR  (missense tolerance ratio, Traynelis 2017): within a window of codons on the human
         MANE transcript, observed missense / (missense + synonymous) alleles in gnomAD v4,
         divided by the same ratio over all possible SNVs. 1 = neutral, < 1 = depleted.
  - RMC  (gnomAD v2 regional missense constraint): o/e of the human region the residue's
         aligned partner falls in. Coarse ranges, but modelled with a proper mutation model.

Pipeline: worm gene -> WBGene -> Alliance orthologs (HGNC) -> gnomAD gene query (GRCh38
variants + ClinVar; GRCh37 RMC) -> Ensembl protein + CDS -> residue-numbering checks ->
pairwise worm/human alignment -> projection -> per-ortholog + consensus tracks.

Outputs (in workingdir):
  {name}.GNOMAD_mtr.txt         continuous track, pos<TAB>consensus MTR (NaN where unmapped)
  {name}.GNOMAD_regions.txt     range track, RMC regions projected to worm coordinates
  {name}.GNOMAD_residues.tsv    per-residue, per-ortholog detail
  {name}.GNOMAD_map.json        ortholog choice, checks, evaluation summary
  gnomad_cache/                 raw API responses (reused on rerun)

Usage
-----
    python scripts/gnomad_constraint.py --gene sma-4 --workingdir out/ --name sma-4
    python scripts/gnomad_constraint.py --input_file x.fa --wbgene WBGene00004858 ...
"""

import argparse
import gzip
import json
import math
import re
import statistics
import sys
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).parent))
from http_retry import request_with_retries

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_ORTHOLOGY = Path.home() / "Downloads" / "ORTHOLOGY-ALLIANCE_TSV_COMBINED.tsv.gz"
GNOMAD_SPACING_S = 4.0  # gnomAD rate-limits bursts; spacing between uncached queries

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
       "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
       "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Sec": "U",
       "Ter": "*"}
HGVSP_RE = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=)$")
RESIDUE_RE = re.compile(r"^([A-Z][a-z]{2})(\d+)$")

_last_gnomad_call = [0.0]


# ---------------------------------------------------------------- HTTP + cache

def _cached_json(cache_path, fetch):
    """Return JSON from cache_path, calling fetch() and saving the result on a miss."""
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    data = fetch()
    cache_path.write_text(json.dumps(data))
    return data


def _cached_text(cache_path, fetch):
    """Return text from cache_path, calling fetch() and saving the result on a miss."""
    if cache_path.exists():
        return cache_path.read_text()
    text = fetch()
    cache_path.write_text(text)
    return text


def gnomad_query(query, cache_path):
    """POST a GraphQL query to gnomAD (cached, rate-spaced); raise on GraphQL errors."""
    def fetch():
        # space uncached calls so a multi-ortholog gene does not trip the rate limit
        wait = GNOMAD_SPACING_S - (time.monotonic() - _last_gnomad_call[0])
        if wait > 0:
            time.sleep(wait)
        resp = request_with_retries("post", GNOMAD_API, retries=3, retry_wait=20,
                                    json={"query": query}, timeout=(10, 180))
        _last_gnomad_call[0] = time.monotonic()
        data = resp.json()
        if data.get("errors"):
            raise RuntimeError(f"gnomAD error: {data['errors'][:2]}")
        return data["data"]
    return _cached_json(cache_path, fetch)


def ensembl_sequence(transcript_id, kind, build, cache_dir):
    """Protein or CDS sequence of an Ensembl transcript on the given build (cached)."""
    tid = transcript_id.split(".")[0]
    url = f"{ENSEMBL[build]}/sequence/id/{tid}"

    def fetch():
        resp = request_with_retries("get", url, params={"type": kind},
                                    headers={"Content-Type": "text/plain"}, timeout=(10, 60))
        return resp.text.strip()
    return _cached_text(cache_dir / f"{tid}.{build}.{kind}.txt", fetch)


def fetch_worm_protein(symbol, cache_dir):
    """Best UniProt entry for a worm gene symbol -> (accession, sequence, WBGene id).

    Many worm genes (par-2, trp-4) have only unreviewed entries, so preference is
    reviewed, then carrying a WormBase gene xref, then longest.
    """
    def fetch():
        q = f"gene_exact:{symbol} AND organism_id:6239"
        resp = request_with_retries("get", UNIPROT_SEARCH, timeout=(10, 60),
                                    params={"query": q, "format": "json", "size": 25})
        return resp.json()
    data = _cached_json(cache_dir / f"uniprot_any_{symbol}.json", fetch)
    if not data.get("results"):
        raise ValueError(f"no UniProt entry for C. elegans gene {symbol}")
    entry = max(data["results"], key=lambda e: (
        "reviewed" in e.get("entryType", "").lower() and "unreviewed" not in
        e.get("entryType", "").lower(), _wbgene_from_entry(e) is not None, e["sequence"]["length"]))
    return entry["primaryAccession"], entry["sequence"]["value"], _wbgene_from_entry(entry)


def _wbgene_from_entry(entry):
    """WormBase gene id (WBGene...) from a UniProt entry's cross-references, or None."""
    for x in entry.get("uniProtKBCrossReferences", []):
        if x.get("database") != "WormBase":
            continue
        for p in x.get("properties", []):
            if p.get("key") == "GeneId":
                return p["value"]
    return None


# ---------------------------------------------------------------- orthology

def load_orthologs(orthology_file, wbgene):
    """All Alliance worm->human rows for a WBGene id, as dicts."""
    key = f"WB:{wbgene}"
    rows = []
    opener = gzip.open if str(orthology_file).endswith(".gz") else open
    with opener(orthology_file, "rt") as fh:
        for line in fh:
            # cheap substring test first; the file is ~1.5M lines
            if key not in line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != key or f[6] != "NCBITaxon:9606":
                continue
            rows.append({"hgnc": f[4], "symbol": f[5], "algorithms": f[8],
                         "n_match": int(f[9]), "n_total": int(f[10]),
                         "best": f[11] == "Yes", "best_rev": f[12] == "Yes"})
    return rows


def select_orthologs(rows, all_orthologs=False):
    """Reciprocal-best rows (IsBestScore and IsBestRevScore); else the max-AlgorithmsMatch set."""
    if all_orthologs or not rows:
        return rows
    recip = [r for r in rows if r["best"] and r["best_rev"]]
    if recip:
        return recip
    top = max(r["n_match"] for r in rows)
    return [r for r in rows if r["n_match"] == top]


# ---------------------------------------------------------------- gnomAD parsing

GRCH38_QUERY = """{ gene(gene_symbol: "%s", reference_genome: GRCh38) {
  gene_id hgnc_id symbol mane_select_transcript { ensembl_id }
  gnomad_constraint { oe_mis mis_z pli }
  variants(dataset: gnomad_r4) { variant_id consequence hgvsp transcript_id
                                 exome { ac } genome { ac } }
  clinvar_variants { variant_id clinical_significance major_consequence hgvsp transcript_id }
} }"""

GRCH37_QUERY = """{ gene(gene_symbol: "%s", reference_genome: GRCh37) {
  canonical_transcript_id
  gnomad_constraint { oe_mis }
  gnomad_v2_regional_missense_constraint { passed_qc has_no_rmc_evidence
    regions { aa_start aa_stop obs_mis exp_mis obs_exp chisq_diff_null } }
} }"""


def parse_hgvsp(hgvsp):
    """'p.Arg361His' -> ('R', 361, 'H'); synonymous 'p.Arg361=' -> ('R', 361, 'R'); else None."""
    if not hgvsp:
        return None
    m = HGVSP_RE.match(hgvsp)
    if not m or m.group(1) not in AA3:
        return None
    ref = AA3[m.group(1)]
    alt = ref if m.group(3) == "=" else AA3.get(m.group(3))
    return (ref, int(m.group(2)), alt) if alt else None


def count_observed(variants, transcript_id, protein, min_ac=1):
    """Distinct missense / synonymous alleles per residue (1-based) on the MANE transcript.

    Returns (mis, syn, check) where check counts hgvsp reference residues that do or do
    not match the human protein, i.e. whether gnomAD's numbering is the sequence in hand.
    """
    tid = transcript_id.split(".")[0]
    n = len(protein)
    mis, syn = [0] * (n + 1), [0] * (n + 1)
    check = {"match": 0, "mismatch": 0, "examples": []}
    for v in variants:
        if (v.get("transcript_id") or "").split(".")[0] != tid:
            continue
        cons = v.get("consequence") or ""
        # only the two classes MTR uses; splice/stop/frameshift are not missense signal
        if cons not in ("missense_variant", "synonymous_variant"):
            continue
        ac = ((v.get("exome") or {}).get("ac") or 0) + ((v.get("genome") or {}).get("ac") or 0)
        p = parse_hgvsp(v.get("hgvsp"))
        # min_ac > 1 drops singletons, which carry little selection signal at v4 sample size
        if ac < max(min_ac, 1) or p is None:
            continue
        ref, pos, _ = p
        if pos > n or protein[pos - 1] != ref:
            check["mismatch"] += 1
            if len(check["examples"]) < 5:
                check["examples"].append(v.get("hgvsp"))
            continue
        check["match"] += 1
        (mis if cons == "missense_variant" else syn)[pos] += 1
    return mis, syn, check


def clinvar_positions(clinvar, transcript_id, protein):
    """Human residues with Pathogenic/Likely pathogenic and with Benign/Likely benign missense."""
    tid = transcript_id.split(".")[0]
    plp, blb = set(), set()
    for v in clinvar:
        if v.get("major_consequence") != "missense_variant":
            continue
        if v.get("transcript_id") and v["transcript_id"].split(".")[0] != tid:
            continue
        p = parse_hgvsp(v.get("hgvsp"))
        if p is None or p[1] > len(protein) or protein[p[1] - 1] != p[0]:
            continue
        sig = (v.get("clinical_significance") or "").lower()
        # conflicting/uncertain calls are neither positive nor negative labels
        if "conflicting" in sig or "uncertain" in sig:
            continue
        if "pathogenic" in sig:
            plp.add(p[1])
        elif "benign" in sig:
            blb.add(p[1])
    return plp, blb


def parse_rmc(rmc, protein37):
    """RMC regions as [(start, stop, o/e, chisq)] on the GRCh37 protein, plus a numbering check.

    Regions arrive in genomic order, so minus-strand genes list them reversed with
    aa_start > aa_stop; both are normalised to ascending residue order.
    """
    regions, check = [], {"match": 0, "mismatch": 0}
    for r in rmc.get("regions") or []:
        ends = []
        for key in ("aa_start", "aa_stop"):
            m = RESIDUE_RE.match(r[key] or "")
            if not m:
                break
            aa, pos = AA3.get(m.group(1)), int(m.group(2))
            ok = pos <= len(protein37) and protein37[pos - 1] == aa
            check["match" if ok else "mismatch"] += 1
            ends.append(pos)
        if len(ends) == 2:
            regions.append((min(ends), max(ends), r["obs_exp"], r["chisq_diff_null"]))
    regions.sort()
    return regions, check


# ---------------------------------------------------------------- MTR

def possible_counts(cds, protein):
    """Possible missense / synonymous SNVs per codon of a CDS (index 1-based by residue)."""
    n = len(protein)
    pm, ps = [0] * (n + 1), [0] * (n + 1)
    for i in range(n):
        codon = cds[3 * i:3 * i + 3].upper()
        # selenocysteine (TGA read through) and ambiguous codons have no clean classes
        if protein[i] == "U" or len(codon) < 3 or set(codon) - set("ACGT"):
            continue
        ref = str(Seq(codon).translate())
        for j in range(3):
            for b in "ACGT":
                if b == codon[j]:
                    continue
                alt = str(Seq(codon[:j] + b + codon[j + 1:]).translate())
                if alt == ref:
                    ps[i + 1] += 1
                elif alt != "*":
                    pm[i + 1] += 1
    return pm, ps


def window_mtr(mis, syn, pm, ps, window):
    """MTR per residue over a centred window, truncated at the protein ends; NaN if no data."""
    n = len(mis) - 1
    half = window // 2
    out = [math.nan] * (n + 1)
    for i in range(1, n + 1):
        lo, hi = max(1, i - half), min(n, i + half)
        om, os_ = sum(mis[lo:hi + 1]), sum(syn[lo:hi + 1])
        qm, qs = sum(pm[lo:hi + 1]), sum(ps[lo:hi + 1])
        # no observed alleles means the ratio is undefined rather than zero
        if om + os_ == 0 or qm + qs == 0 or qm == 0:
            continue
        out[i] = (om / (om + os_)) / (qm / (qm + qs))
    return out


# ---------------------------------------------------------------- alignment

def make_aligner():
    """Global BLOSUM62 aligner with free end gaps (worm and human termini often differ)."""
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    al.end_gap_score = 0
    al.mode = "global"
    return al


def residue_map(seq_a, seq_b, aligner=None):
    """Map 1-based positions of seq_a to seq_b (None in gaps) from one optimal alignment."""
    # selenocysteine is not in BLOSUM62; align it as its sulfur analogue Cys
    aln = (aligner or make_aligner()).align(seq_a.replace("U", "C"), seq_b.replace("U", "C"))[0]
    m = [None] * (len(seq_a) + 1)
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(int(a1 - a0)):
            m[int(a0) + k + 1] = int(b0) + k + 1
    return m


def local_identity(worm, human, wmap, half=5):
    """Fraction of identical aligned pairs in a +/-half window; unaligned residues count as 0."""
    n = len(worm)
    ident = [0] * (n + 1)
    for i in range(1, n + 1):
        h = wmap[i]
        ident[i] = 1 if h and worm[i - 1] == human[h - 1] else 0
    out = [math.nan] * (n + 1)
    for i in range(1, n + 1):
        lo, hi = max(1, i - half), min(n, i + half)
        out[i] = sum(ident[lo:hi + 1]) / (hi - lo + 1)
    return out


# ---------------------------------------------------------------- per-ortholog

def process_ortholog(symbol, worm_seq, cache_dir, window, min_identity, aligner, min_ac=1):
    """Fetch, check and project one human ortholog; returns a dict of per-worm-residue arrays."""
    res = {"symbol": symbol, "status": "ok"}
    g38 = gnomad_query(GRCH38_QUERY % symbol, cache_dir / f"gnomad38_{symbol}.json")["gene"]
    g37 = gnomad_query(GRCH37_QUERY % symbol, cache_dir / f"gnomad37_{symbol}.json")["gene"]
    if not g38 or not g38.get("mane_select_transcript"):
        res["status"] = "no_mane_transcript"
        return res
    mane = g38["mane_select_transcript"]["ensembl_id"]
    prot = ensembl_sequence(mane, "protein", "GRCh38", cache_dir)
    cds = ensembl_sequence(mane, "cds", "GRCh38", cache_dir)
    res.update(mane=mane, human_length=len(prot), gene_constraint=g38.get("gnomad_constraint"))

    # CDS must translate to the protein, or per-codon possible counts are misaligned
    tr = str(Seq(cds[:3 * len(prot)]).translate())
    # Sec (U) reads as a TGA stop; a non-AUG start (TXNRD3 CTG) still encodes Met
    cds_ok = all(a == b or (b == "U" and a == "*") or (k == 0 and b == "M")
                 for k, (a, b) in enumerate(zip(tr, prot)))
    res["cds_translates"] = cds_ok

    mis, syn, check = count_observed(g38.get("variants") or [], mane, prot, min_ac)
    res["hgvsp_check"] = check
    res["n_missense_alleles"], res["n_synonymous_alleles"] = sum(mis), sum(syn)
    if check["mismatch"] or not cds_ok:
        res["status"] = "numbering_mismatch"
        return res
    pm, ps = possible_counts(cds, prot)
    mtr = window_mtr(mis, syn, pm, ps, window)
    plp, blb = clinvar_positions(g38.get("clinvar_variants") or [], mane, prot)
    res["n_clinvar_plp_residues"], res["n_clinvar_blb_residues"] = len(plp), len(blb)

    # RMC lives on the GRCh37 canonical transcript; map it to MANE when the proteins differ
    rmc_oe = [math.nan] * (len(prot) + 1)
    rmc_regions_h38 = []
    rmc = (g37 or {}).get("gnomad_v2_regional_missense_constraint")
    res["rmc_status"] = "absent"
    if g37 and rmc and rmc.get("passed_qc"):
        prot37 = ensembl_sequence(g37["canonical_transcript_id"], "protein", "GRCh37", cache_dir)
        res["grch37_protein_identical"] = prot37 == prot
        if rmc.get("has_no_rmc_evidence"):
            # no regional split: the whole protein carries the v2 gene-level o/e
            oe = (g37.get("gnomad_constraint") or {}).get("oe_mis")
            regions37, rcheck = [(1, len(prot37), oe, None)], {"match": 1, "mismatch": 0}
            res["rmc_status"] = "single_region"
        else:
            regions37, rcheck = parse_rmc(rmc, prot37)
            res["rmc_status"] = "regional"
        res["rmc_check"] = rcheck
        h37_to_h38 = (list(range(len(prot) + 1)) if prot37 == prot
                      else residue_map(prot37, prot, aligner))
        for s, e, oe, chisq in regions37:
            mapped = [h37_to_h38[p] for p in range(s, e + 1) if h37_to_h38[p]]
            if not mapped or oe is None:
                continue
            rmc_regions_h38.append((min(mapped), max(mapped), oe, chisq))
            for p in mapped:
                rmc_oe[p] = oe
    res["rmc_regions_human"] = rmc_regions_h38

    # project human tracks to worm through the alignment
    wmap = residue_map(worm_seq, prot, aligner)
    ident = local_identity(worm_seq, prot, wmap)
    n = len(worm_seq)
    res["worm_to_human"] = wmap
    res["identity"] = ident
    res["mtr"] = [math.nan] * (n + 1)
    res["rmc_oe"] = [math.nan] * (n + 1)
    res["plp"] = [False] * (n + 1)
    res["blb"] = [False] * (n + 1)
    for i in range(1, n + 1):
        h = wmap[i]
        if not h:
            continue
        # loose (RMC) projection keeps every aligned residue; tight (MTR) needs local identity
        res["rmc_oe"][i] = rmc_oe[h]
        if ident[i] >= min_identity:
            res["mtr"][i] = mtr[h]
        res["plp"][i] = h in plp
        res["blb"][i] = h in blb
    aligned = sum(1 for i in range(1, n + 1) if wmap[i])
    res["fraction_worm_aligned"] = aligned / n
    res["fraction_worm_tight"] = sum(1 for i in range(1, n + 1)
                                     if not math.isnan(res["mtr"][i])) / n
    res["global_identity"] = (sum(1 for i in range(1, n + 1)
                                  if wmap[i] and worm_seq[i - 1] == prot[wmap[i] - 1])
                              / max(aligned, 1))
    return res


def project_regions(wmap, regions_human, n):
    """Worm ranges covering the aligned partners of each human RMC region."""
    out = []
    for s, e, oe, chisq in regions_human:
        hits = [i for i in range(1, n + 1) if wmap[i] and s <= wmap[i] <= e]
        if hits:
            out.append({"start": min(hits), "end": max(hits), "oe": oe, "chisq": chisq,
                        "human": [s, e],
                        "aligned_fraction": len(hits) / (max(hits) - min(hits) + 1)})
    return out


# ---------------------------------------------------------------- evaluation

def auc(pos_scores, neg_scores):
    """Mann-Whitney AUC: P(score of a positive > score of a negative), ties count half."""
    if not pos_scores or not neg_scores:
        return None
    wins = 0.0
    for p in pos_scores:
        for q in neg_scores:
            wins += 1.0 if p > q else 0.5 if p == q else 0.0
    return wins / (len(pos_scores) * len(neg_scores))


def classify(value, constrained, tolerant):
    """Three-way class from a depletion value (lower = more constrained)."""
    if value is None or math.isnan(value):
        return "unmapped"
    if value <= constrained:
        return "constrained"
    return "tolerant" if value >= tolerant else "intermediate"


def nanmedian(vals):
    """Median of the non-NaN values, or NaN."""
    v = [x for x in vals if not math.isnan(x)]
    return statistics.median(v) if v else math.nan


def evaluate(cons_mtr, cons_rmc, cons_ident, plp, n, args):
    """AUCs of each signal for ClinVar P/LP partner residues, plus class enrichment."""
    out = {}
    # AUCs compare P/LP residues against all other residues scored by that signal;
    # signals are negated where low = constrained so higher always predicts P/LP
    for name, vals, sign in (("mtr", cons_mtr, -1), ("rmc_oe", cons_rmc, -1),
                             ("worm_human_identity", cons_ident, 1)):
        pos = [sign * vals[i] for i in range(1, n + 1) if plp[i] and not math.isnan(vals[i])]
        neg = [sign * vals[i] for i in range(1, n + 1)
               if not plp[i] and not math.isnan(vals[i])]
        a = auc(pos, neg)
        out[f"auc_{name}"] = round(a, 3) if a is not None else None
        out[f"n_plp_{name}"] = len(pos)
    # fraction of P/LP residues landing in each MTR class vs the background
    for name, vals, c, t in (("mtr", cons_mtr, args.mtr_constrained, args.mtr_tolerant),
                             ("rmc", cons_rmc, args.rmc_constrained, args.rmc_tolerant)):
        classes = [classify(vals[i], c, t) for i in range(1, n + 1)]
        bg = {k: classes.count(k) / n for k in ("constrained", "intermediate", "tolerant",
                                                 "unmapped")}
        pl = [classes[i - 1] for i in range(1, n + 1) if plp[i]]
        out[f"{name}_class_fraction"] = {k: round(v, 3) for k, v in bg.items()}
        out[f"{name}_plp_constrained_fraction"] = (round(pl.count("constrained") / len(pl), 3)
                                                   if pl else None)
        # enrichment among scored residues only, so unmapped regions do not inflate it;
        # robust where AUC is not (RMC has few distinct values, so near-ties decide AUC)
        scored = [i for i in range(1, n + 1) if classes[i - 1] != "unmapped"]
        pl_s = [i for i in scored if plp[i]]
        bg_c = sum(1 for i in scored if classes[i - 1] == "constrained") / max(len(scored), 1)
        pl_c = (sum(1 for i in pl_s if classes[i - 1] == "constrained") / len(pl_s)
                if pl_s else None)
        out[f"{name}_plp_enrichment"] = round(pl_c / bg_c, 2) if pl_c is not None and bg_c \
            else None
        # tag-site relevance: class of the terminal 15 residues
        out[f"{name}_nterm15"] = _majority(classes[:15])
        out[f"{name}_cterm15"] = _majority(classes[-15:])
    return out


def _majority(classes):
    """Most common class in a list (ties broken by first appearance)."""
    return max(dict.fromkeys(classes), key=classes.count) if classes else None


# ---------------------------------------------------------------- main

def _fmt(x):
    """Tab-file number formatting with NaN passthrough."""
    return "NaN" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.4g}"


def main(workingdir, name, gene="", input_file="", wbgene="", orthology=str(DEFAULT_ORTHOLOGY),
         window=31, min_identity=0.3, all_orthologs=False, min_ac=1, mtr_quantile=0.25,
         mtr_constrained=0.7,
         mtr_tolerant=0.95, rmc_constrained=0.4, rmc_tolerant=0.8):
    """Build the gnomAD constraint tracks for one worm protein; returns the summary dict."""
    args = argparse.Namespace(mtr_constrained=mtr_constrained, mtr_tolerant=mtr_tolerant,
                              rmc_constrained=rmc_constrained, rmc_tolerant=rmc_tolerant)
    wd = Path(workingdir)
    wd.mkdir(parents=True, exist_ok=True)
    cache_dir = wd / "gnomad_cache"
    cache_dir.mkdir(exist_ok=True)
    t0 = time.monotonic()

    # worm sequence: explicit FASTA wins; otherwise the reviewed UniProt entry for the symbol
    accession = None
    if input_file:
        worm_seq = str(next(SeqIO.parse(input_file, "fasta")).seq).rstrip("*")
    else:
        accession, worm_seq, wb = fetch_worm_protein(gene, cache_dir)
        wbgene = wbgene or wb
    if not wbgene:
        raise ValueError("no WBGene id: pass --wbgene or a --gene with a UniProt WormBase xref")
    n = len(worm_seq)
    summary = {"gene": gene, "wbgene": wbgene, "uniprot": accession, "worm_length": n,
               "window": window, "min_identity": min_identity, "orthology_file": str(orthology)}

    rows = load_orthologs(orthology, wbgene)
    chosen = select_orthologs(rows, all_orthologs)
    summary["alliance_rows"] = rows
    summary["orthologs_used"] = [r["symbol"] for r in chosen]
    if not chosen:
        summary["status"] = "no_human_ortholog"
        (wd / f"{name}.GNOMAD_map.json").write_text(json.dumps(summary, indent=1))
        return summary

    aligner = make_aligner()
    results = []
    # one ortholog failing (missing MANE, numbering mismatch) must not sink the others
    for r in chosen:
        try:
            results.append(process_ortholog(r["symbol"], worm_seq, cache_dir, window,
                                            min_identity, aligner, min_ac))
        except Exception as exc:  # noqa: BLE001 - recorded per ortholog, run continues
            results.append({"symbol": r["symbol"], "status": f"error: {exc}"})
    good = [r for r in results if r["status"] == "ok"]

    # consensus across orthologs: per-residue median; P/LP if any ortholog's partner is P/LP
    cons_mtr = [math.nan] + [nanmedian([r["mtr"][i] for r in good]) for i in range(1, n + 1)]
    cons_rmc = [math.nan] + [nanmedian([r["rmc_oe"][i] for r in good]) for i in range(1, n + 1)]
    cons_ident = [math.nan] + [nanmedian([r["identity"][i] for r in good])
                               for i in range(1, n + 1)]
    plp = [False] + [any(r["plp"][i] for r in good) for i in range(1, n + 1)]

    # absolute MTR mostly separates genes (SMAD4 is low everywhere); for choosing sites
    # within one protein the within-gene rank matters, so classes default to quantiles
    scored = sorted(x for x in cons_mtr[1:] if not math.isnan(x))
    if mtr_quantile and scored:
        mtr_constrained = scored[int(mtr_quantile * (len(scored) - 1))]
        mtr_tolerant = scored[int((1 - mtr_quantile) * (len(scored) - 1))]
        args.mtr_constrained, args.mtr_tolerant = mtr_constrained, mtr_tolerant
    summary.update(min_ac=min_ac, mtr_quantile=mtr_quantile,
                   mtr_thresholds=[mtr_constrained, mtr_tolerant])
    n_constr = [0] + [sum(1 for r in good if classify(r["mtr"][i], mtr_constrained,
                                                       mtr_tolerant) == "constrained")
                      for i in range(1, n + 1)]

    with open(wd / f"{name}.GNOMAD_mtr.txt", "w") as fh:
        fh.writelines(f"{i}\t{_fmt(cons_mtr[i])}\n" for i in range(1, n + 1))

    with open(wd / f"{name}.GNOMAD_regions.txt", "w") as fh:
        for r in good:
            r["rmc_regions_worm"] = project_regions(r["worm_to_human"], r["rmc_regions_human"], n)
            for reg in r["rmc_regions_worm"]:
                cls = classify(reg["oe"], rmc_constrained, rmc_tolerant)
                fh.write(f"GnomAD_RMC\t{reg['start']}\t{reg['end']}\t{r['symbol']} "
                         f"{reg['human'][0]}-{reg['human'][1]} o/e={reg['oe']:.2f} ({cls})\n")

    with open(wd / f"{name}.GNOMAD_residues.tsv", "w") as fh:
        cols = ["pos", "aa", "consensus_mtr", "consensus_rmc_oe", "mtr_class", "rmc_class",
                "n_orthologs_constrained", "clinvar_plp"]
        for r in good:
            s = r["symbol"]
            cols += [f"{s}_pos", f"{s}_aa", f"{s}_identity", f"{s}_mtr", f"{s}_rmc_oe"]
        fh.write("\t".join(cols) + "\n")
        for i in range(1, n + 1):
            row = [str(i), worm_seq[i - 1], _fmt(cons_mtr[i]), _fmt(cons_rmc[i]),
                   classify(cons_mtr[i], mtr_constrained, mtr_tolerant),
                   classify(cons_rmc[i], rmc_constrained, rmc_tolerant),
                   str(n_constr[i]), str(int(plp[i]))]
            for r in good:
                h = r["worm_to_human"][i]
                row += [str(h or ""), "", _fmt(r["identity"][i]), _fmt(r["mtr"][i]),
                        _fmt(r["rmc_oe"][i])]
            fh.write("\t".join(row) + "\n")

    # keep the JSON small: drop the per-residue arrays, keep the per-ortholog scalars
    drop = {"worm_to_human", "identity", "mtr", "rmc_oe", "plp", "blb"}
    summary["orthologs"] = [{k: v for k, v in r.items() if k not in drop} for r in results]
    summary["status"] = "ok" if good else "no_usable_ortholog"
    if good:
        summary["evaluation"] = evaluate(cons_mtr, cons_rmc, cons_ident, plp, n, args)
    summary["elapsed_s"] = round(time.monotonic() - t0, 1)
    (wd / f"{name}.GNOMAD_map.json").write_text(json.dumps(summary, indent=1))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--gene", default="", help="C. elegans gene symbol (e.g. sma-4)")
    parser.add_argument("-f", "--input_file", default="", help="worm protein FASTA (optional)")
    parser.add_argument("--wbgene", default="", help="WBGene id; required with --input_file")
    parser.add_argument("--orthology", default=str(DEFAULT_ORTHOLOGY),
                        help="Alliance ORTHOLOGY-ALLIANCE_TSV_COMBINED.tsv(.gz)")
    parser.add_argument("--workingdir", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--window", type=int, default=31, help="MTR window in codons")
    parser.add_argument("--min_identity", type=float, default=0.3,
                        help="local worm/human identity needed for a per-residue MTR value")
    parser.add_argument("--min_ac", type=int, default=1,
                        help="minimum gnomAD allele count for a variant to count (2 = no singletons)")
    parser.add_argument("--mtr_quantile", type=float, default=0.25,
                        help="MTR classes from within-gene quantiles (0 = absolute thresholds)")
    parser.add_argument("--all_orthologs", action="store_true",
                        help="use every Alliance human ortholog, not only reciprocal-best")
    a = parser.parse_args()
    s = main(a.workingdir, a.name, a.gene, a.input_file, a.wbgene, a.orthology, a.window,
             a.min_identity, a.all_orthologs, a.min_ac, a.mtr_quantile)
    print(json.dumps({k: s.get(k) for k in ("gene", "status", "orthologs_used", "evaluation",
                                            "elapsed_s")}, indent=1))
