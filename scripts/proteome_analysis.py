"""
proteome_analysis.py

Read-only analyses over proteome.sqlite3 (plus the modification files in the batch run folders):

  census    genes, products (canonical entries, isoform records, WormBase-only), products
            per gene, AlphaFold DB coverage
  wormtag   overlap with the WormTagDB export (genes already carrying a fluorescent tag)
  termini   classification of every product / gene as C-terminal, N-terminal or internal

Termini rules (per product, per end; a window of --window residues at that end):
  * no blocking UniProt feature, modification site or DeepTMHMM signal / TM helix overlaps it
  * mean conservation (JSD) < 0.5 and mean pLDDT < 50   (the app's own thresholds)
  * the terminal residue is cytosolic ("inside") or the protein is secreted
  * the window is present, uninterrupted, in every isoform listed for the protein
Order: C-terminus first, then N-terminus, else internal. A guide cut within 10 bp of the
terminal insertion point is only a warning and never changes the class. Missing data (no
structure, no topology, no conservation) passes the filter and is counted separately.
A gene takes a terminal class only if all of its products pass at that end.

Outputs go to <run_dir>/analysis/: termini_products.tsv, termini_genes.tsv, summary.txt.

Matt Rich, 2026
"""

import argparse
import collections
import sqlite3
import statistics
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import build_proteome_db as bdb  # noqa: E402

WORMTAG = Path(__file__).parent.parent.parent / "TAG_ORTHOLOGS" / "wormtagdb_tagged_genes.tsv"
CONS_MAX = 0.5      # conservation (JSD) must be below this, as in scores.config.json
PLDDT_MAX = 50.0    # pLDDT (0-100) must be below this; scores.config.json uses 0.5 of 0-1
GUIDE_BP = 10       # guide cut within this many bp of the insertion point
WINDOWS = [1, 5, 10, 20]
BAD_TOPO = {"signal", "TMhelix"}


def open_db(path):
    """Read-only, immutable connection (SQLite locking is unreliable on the network volume)."""
    return sqlite3.connect("file:{}?immutable=1".format(path), uri=True)


def gene_key(row):
    """Gene identifier: WBGene id when known, else the protein's own accession/id."""
    return row.wb_gene or row.wormbase_gene or "acc:" + (row.accession or row.id)


def load_proteins(conn):
    """Proteins table with a gene key column."""
    df = pd.read_sql("SELECT pid, id, kind, accession, gene_name, wb_gene, wormbase_gene, "
                     "wb_transcript, seq_length, has_structure FROM proteins", conn)
    df["gene"] = [gene_key(r) for r in df.itertuples()]
    return df


# ── census ────────────────────────────────────────────────────────────────────

def dist_text(values):
    """Mean / median / max and a 1,2,3,4-5,6+ histogram of a list of counts."""
    bins = collections.Counter("1" if v == 1 else "2" if v == 2 else "3" if v == 3
                               else "4-5" if v <= 5 else "6+" for v in values)
    return "mean {:.2f}, median {}, max {}; per-gene histogram {}".format(
        statistics.mean(values), statistics.median(values), max(values),
        ", ".join("{}: {}".format(k, bins[k]) for k in ("1", "2", "3", "4-5", "6+")))


def census(conn, prot):
    """Lines describing gene, product and isoform counts and AFDB coverage."""
    out = ["== census =="]
    out.append("products (rows in proteins): {}".format(len(prot)))
    for kind, n in prot.groupby("kind").size().items():
        out.append("  {}: {}".format(kind, n))
    genes = prot.groupby("gene")
    wb = prot[prot.gene.str.startswith("WBGene")]
    out.append("genes: {} distinct ({} WBGene ids, {} products without a WBGene id kept as "
               "their own entry)".format(prot.gene.nunique(), wb.gene.nunique(),
                                         int((~prot.gene.str.startswith("WBGene")).sum())))
    out.append("canonical (kind=uniprot) entries: {}, genes with one: {}".format(
        int((prot.kind == "uniprot").sum()), prot[prot.kind == "uniprot"].gene.nunique()))
    out.append("products per gene: " + dist_text(genes.size().tolist()))
    # isoforms per gene: every accession seen for the gene, in its own rows or in the
    # isoform lists the conservation step recorded for any of its canonical entries
    iso = pd.read_sql("SELECT pid, accession FROM isoforms WHERE accession IS NOT NULL", conn)
    iso = iso.merge(prot[["pid", "gene"]], on="pid")
    names = {}
    for r in prot.itertuples():
        names.setdefault(r.gene, set()).add(r.accession or r.id)
    for r in iso.itertuples():
        names[r.gene].add(r.accession)
    counts = [len(v) for v in names.values()]
    out.append("distinct isoform accessions per gene (own rows + isoform lists): "
               + dist_text(counts))
    n_listed = iso.pid.nunique()
    out.append("canonical entries that carry an isoform list: {} of {}".format(
        n_listed, int((prot.kind == "uniprot").sum())))
    # AFDB coverage
    out.append("")
    out.append("AlphaFold DB coverage (has_structure):")
    for kind, grp in prot.groupby("kind"):
        out.append("  {:<16}{:>6} / {:<6}({:.1f}%)".format(
            kind, int(grp.has_structure.sum()), len(grp), 100 * grp.has_structure.mean()))
    out.append("  all products    {:>6} / {:<6}({:.1f}%)".format(
        int(prot.has_structure.sum()), len(prot), 100 * prot.has_structure.mean()))
    res = (prot.has_structure * prot.seq_length).sum() / prot.seq_length.sum()
    out.append("  residues in modelled proteins: {:.1f}%".format(100 * res))
    g_any = genes.has_structure.max().sum()
    g_all = genes.has_structure.min().sum()
    out.append("  genes with a model for any product: {} / {} ({:.1f}%); for every product: "
               "{} ({:.1f}%)".format(int(g_any), len(genes), 100 * g_any / len(genes),
                                     int(g_all), 100 * g_all / len(genes)))
    canon = prot[prot.kind == "uniprot"].groupby("gene").has_structure.max()
    out.append("  genes by canonical entry with a model: {} / {}".format(int(canon.sum()),
                                                                       len(canon)))
    return out


# ── WormTagDB ─────────────────────────────────────────────────────────────────

def match_tags(prot, tag):
    """{db gene key: tag row indexes} matching on WBGene id, symbol or sequence name."""
    names = collections.defaultdict(set)
    for r in prot.itertuples():
        for tok in (r.wb_gene, r.wormbase_gene, r.gene_name, r.wb_transcript):
            if tok:
                names[tok.lower()].add(r.gene)
    hits = collections.defaultdict(set)
    for i, r in tag.iterrows():
        # WormTagDB symbols can be joined, e.g. "cdf-3/ZK185.5"; any part may match
        for tok in [r.wbid] + r.gene.split("/"):
            for key in names.get(tok.lower(), ()):
                hits[key].add(i)
    return hits


def wormtag(prot):
    """Lines comparing the database's genes to the WormTagDB export; returns FP gene keys."""
    out = ["", "== WormTagDB ({}) ==".format(WORMTAG.name)]
    tag = pd.read_csv(WORMTAG, sep="\t", keep_default_na=False)
    is_fp = tag.fluor_tags != ""
    is_cgc = tag.cgc_strains != ""
    hits = match_tags(prot, tag)
    genes = set(prot.gene)
    out.append("export rows (one per tagged gene): {} with any tag, {} with a fluorescent "
               "protein (FP), {} with a CGC strain".format(len(tag), int(is_fp.sum()),
                                                           int(is_cgc.sum())))
    # the symbol-list view: a plain gene list matched on the export's symbol column
    symbols = WORMTAG.parent / "data" / "celegans_proteome.tsv"
    if symbols.exists():
        listed = {ln.strip() for ln in open(symbols) if ln.strip()}
        any_n, fp_n = tag.gene.isin(listed).sum(), tag[is_fp].gene.isin(listed).sum()
        out.append("{:<20}{:>8}{:>9}{:>12}{:>8}{:>8}".format(
            "gene_list", "n_genes", "any_tag", "any_tag_pct", "FP_tag", "FP_pct"))
        out.append("{:<20}{:>8}{:>9}{:>12}{:>8}{:>8}".format(
            "celegans_proteome", len(listed), any_n, "{:.1f}%".format(100 * any_n / len(listed)),
            fp_n, "{:.1f}%".format(100 * fp_n / len(listed))))
    fp_keys = set()
    out.append("database gene entries ({}; a different universe from the list above), matched on "
               "WBGene id, symbol or sequence name:".format(len(genes)))
    for label, mask in (("any tag", pd.Series(True, index=tag.index)), ("FP tag", is_fp),
                        ("CGC strain", is_cgc)):
        keys = {k for k, idx in hits.items() if any(mask[i] for i in idx)}
        if label == "FP tag":
            fp_keys = keys
        out.append("  {:<11}{:>6} ({:.1f}%)".format(label, len(keys), 100 * len(keys) / len(genes)))
    matched = set().union(*hits.values()) if hits else set()
    out.append("  export rows with no database match: {} (FP: {})".format(
        len(tag) - len(matched), int(is_fp.sum() - is_fp[list(matched)].sum())))
    return out, fp_keys


# ── termini ───────────────────────────────────────────────────────────────────

def load_modifications(prot, status, run_dir):
    """{pid: [(start, stop)]} modification ranges read from the run folders (inclusive stop)."""
    mods = {}
    for r in prot.itertuples():
        if status.get(r.id, {}).get("tasks", {}).get("modifications") != "ok":
            continue
        path = Path(run_dir) / r.id / (r.id + bdb.SUFFIX["modifications"])
        if path.exists() and path.stat().st_size:
            # stop_shift -1: regex_sites writes match.end()+1 (see build_proteome_db)
            mods[r.pid] = [(s, e) for _, _, s, e, _ in bdb._range_features(
                str(path), "modification", -1)]
    return mods


def load_features(conn):
    """Per-pid blocking-UniProt, UniProt-site and topology intervals."""
    feats = pd.read_sql("SELECT pid, track, source, start, stop, description FROM features "
                        "WHERE track IN ('uniprot', 'topology')", conn)
    block, site, topo = (collections.defaultdict(list) for _ in range(3))
    for r in feats.itertuples():
        if r.track == "topology":
            topo[r.pid].append((r.start, r.stop, r.description))
        elif r.source == "UniProt":
            block[r.pid].append((r.start, r.stop, r.description))
        else:
            site[r.pid].append((r.start, r.stop, r.description))
    return block, site, topo


def load_terminal_residues(conn, width):
    """{pid: {pos: (conservation, pLDDT)}} for the first and last `width` residues."""
    df = pd.read_sql("SELECT r.pid, r.pos, r.conservation, r.plddt FROM residues r "
                     "JOIN proteins p USING (pid) WHERE r.pos <= ? OR r.pos > p.seq_length - ?",
                     conn, params=(width, width))
    cols = collections.defaultdict(dict)
    for r in df.itertuples():
        cons = None if pd.isna(r.conservation) else r.conservation / bdb.SCALE["conservation"]
        pl = None if pd.isna(r.plddt) else r.plddt / bdb.SCALE["plddt"]
        cols[r.pid][r.pos] = (cons, pl)
    return cols


def load_isoform_spans(conn):
    """{pid: {iso_idx: {"present": [...], "inserts": [after...]}}} for non-query isoforms."""
    query_idx = dict(conn.execute("SELECT pid, iso_idx FROM isoforms WHERE is_query = 1"))
    spans = collections.defaultdict(lambda: collections.defaultdict(
        lambda: {"present": [], "inserts": []}))
    for pid, idx, kind, start, stop in conn.execute(
            "SELECT pid, iso_idx, kind, start, stop FROM isoform_spans"):
        if idx == query_idx.get(pid):
            continue   # the query itself is not an alternative isoform
        slot = spans[pid][idx]
        if kind == "present":
            slot["present"].append((start, stop))
        elif kind == "insert":
            slot["inserts"].append(start)
    # an isoform with no rows at all (fully skipped) must still count, so list all of them
    for pid, idx in conn.execute("SELECT pid, iso_idx FROM isoforms WHERE is_query = 0"):
        spans[pid][idx]
    return spans


def window_of(end, length, width):
    """(first, last) residue of the terminal window; the N-terminus is residues 1..width."""
    width = min(width, length)
    return (length - width + 1, length) if end == "C" else (1, width)


def overlaps(intervals, lo, hi):
    """True when any (start, stop, ...) interval overlaps [lo, hi]."""
    return any(iv[0] <= hi and iv[1] >= lo for iv in intervals)


def isoform_ok(iso, lo, hi, end):
    """True when [lo, hi] is present and unbroken in every alternative isoform of a protein."""
    for slot in iso.values():
        covered = sum(min(e, hi) - max(s, lo) + 1 for s, e in slot["present"]
                      if s <= hi and e >= lo)
        if covered < hi - lo + 1:
            return False   # part of the window is skipped in this isoform
        # an insertion after residue p splits the window; after 0 extends the N-terminus
        first, last = (lo, hi) if end == "C" else (0, hi - 1)
        if any(first <= a <= last for a in slot["inserts"]):
            return False
    return True


def mean_or_none(values):
    """Mean of the non-None values, or None when there are none."""
    values = [v for v in values if v is not None]
    return sum(values) / len(values) if values else None


def classify_end(p, end, width, data, with_site, ignore=()):
    """Filter results for one terminus of one protein: {check: True/False/None}."""
    length = p.seq_length
    lo, hi = window_of(end, length, width)
    term_pos = length if end == "C" else 1
    block, site, topo, mods, resid, iso = data
    topo_iv = topo.get(p.pid, [])
    ivs = [(s, e) for s, e, _ in block.get(p.pid, [])]
    ivs += [(s, e) for s, e in mods.get(p.pid, [])]
    ivs += [(s, e) for s, e, d in topo_iv if d in BAD_TOPO]
    if with_site:
        ivs += [(s, e) for s, e, _ in site.get(p.pid, [])]
    res = {"feature_free": not overlaps(ivs, lo, hi)}
    vals = [resid.get(p.pid, {}).get(pos) for pos in range(lo, hi + 1)]
    cons = mean_or_none([v[0] for v in vals if v])
    plddt = mean_or_none([v[1] for v in vals if v])
    # None = no data for the window; it passes but is counted separately
    res["conservation"] = None if cons is None else cons < CONS_MAX
    res["plddt"] = None if plddt is None else plddt < PLDDT_MAX
    if topo_iv:
        label = next((d for s, e, d in topo_iv if s <= term_pos <= e), None)
        has_tm = any(d == "TMhelix" for _, _, d in topo_iv) or any(
            "Transmembrane" in d for _, _, d in block.get(p.pid, []))
        has_signal = any(d == "signal" for _, _, d in topo_iv) or any(
            "Signal" in d for _, _, d in block.get(p.pid, []))
        res["topology"] = label == "inside" or (has_signal and not has_tm)
    else:
        res["topology"] = None
    res["isoforms"] = isoform_ok(iso[p.pid], lo, hi, end) if p.pid in iso else None
    # ignored checks stay out of the pass/fail decision (the columns are dropped)
    return {k: v for k, v in res.items() if k not in ignore}


def passes(res):
    """A terminus passes when no check is False (None, i.e. no data, passes)."""
    return all(v is not False for v in res.values())


def guide_distances(conn, prot):
    """{(pid, 'N'|'C'): nearest guide cut distance} for designed proteins, else absent."""
    if conn is None or not bdb.table_in(conn, "site_guides"):
        return {}
    length = dict(zip(prot.pid, prot.seq_length))
    designed = {r[0] for r in conn.execute("SELECT pid FROM reagent_status WHERE status='success'")}
    out = {}
    for pid, res, dist in conn.execute("SELECT pid, residue_index, MIN(distance) FROM site_guides "
                                       "GROUP BY pid, residue_index"):
        if pid in length and res == length[pid]:
            out[(pid, "C")] = dist
        elif res == 1:
            out[(pid, "N")] = dist
    # a designed protein with no guide row at its terminus is recorded as infinite distance
    for pid in designed:
        out.setdefault((pid, "C"), float("inf"))
        out.setdefault((pid, "N"), float("inf"))
    return out


def classify_products(prot, data, guides, width, with_site=False, ignore=()):
    """One row per product: per-end checks, class (C/N/internal) and guide status."""
    rows = []
    for p in prot.itertuples():
        row = {"pid": p.pid, "id": p.id, "kind": p.kind, "gene": p.gene, "length": p.seq_length}
        cls = "internal"
        for end in ("C", "N"):
            res = classify_end(p, end, width, data, with_site, ignore)
            row.update({"{}_{}".format(end, k): v for k, v in res.items()})
            row[end + "_pass"] = passes(res)
        # C first, then N; the guide check is a warning on the chosen terminus only
        if row["C_pass"]:
            cls = "C"
        elif row["N_pass"]:
            cls = "N"
        row["class"] = cls
        d = guides.get((p.pid, cls)) if cls != "internal" else None
        row["guide_bp"] = d
        row["guide_status"] = ("n/a" if cls == "internal" else "not designed" if d is None
                               else "ok" if d <= GUIDE_BP else "WARN no guide <=%dbp" % GUIDE_BP)
        rows.append(row)
    return pd.DataFrame(rows)


def classify_genes(prods):
    """One row per gene: C if every product passes at C, else N if every product passes at N."""
    rows = []
    for gene, grp in prods.groupby("gene"):
        cls = "C" if grp.C_pass.all() else "N" if grp.N_pass.all() else "internal"
        status = "n/a"
        if cls != "internal":
            d = [x for x in grp.guide_bp if x is not None and not pd.isna(x)]
            # the gene is fine when any designed product has a guide in range
            status = ("not designed" if not d else "ok" if min(d) <= GUIDE_BP
                      else "WARN no guide <=%dbp" % GUIDE_BP)
        rows.append({"gene": gene, "n_products": len(grp), "class": cls, "guide_status": status,
                     "ids": ";".join(grp.id)})
    return pd.DataFrame(rows)


def class_lines(label, df):
    """Counts of each class and guide status for a products/genes table."""
    out = ["{} (n={}):".format(label, len(df))]
    for cls in ("C", "N", "internal"):
        sub = df[df["class"] == cls]
        extra = ""
        if cls != "internal":
            extra = "  guide: " + ", ".join("{} {}".format(k, v) for k, v in
                                             sub.guide_status.value_counts().items())
        out.append("  {:<9}{:>6} ({:.1f}%){}".format(cls, len(sub), 100 * len(sub) / len(df),
                                                    extra))
    return out


def filter_lines(prods):
    """Per-filter pass / fail / no-data counts for each terminus (canonical entries)."""
    out = ["per-filter results, canonical entries (pass / fail / no data):"]
    canon = prods[prods.kind == "uniprot"]
    for end in ("C", "N"):
        for chk in ("feature_free", "conservation", "plddt", "topology", "isoforms"):
            col = canon["{}_{}".format(end, chk)]
            out.append("  {}-term {:<13}{:>6} / {:>6} / {:>6}".format(
                end, chk, int((col == True).sum()), int((col == False).sum()),  # noqa: E712
                int(col.isna().sum())))
        out.append("  {}-term all pass     {:>6}".format(end, int(canon[end + "_pass"].sum())))
    return out


def window_means(prot, end, width, resid):
    """(conservation, pLDDT) window means per protein as float arrays, NaN where no data."""
    cons, plddt = [], []
    for p in prot.itertuples():
        lo, hi = window_of(end, p.seq_length, width)
        vals = [resid.get(p.pid, {}).get(pos) for pos in range(lo, hi + 1)]
        c = mean_or_none([v[0] for v in vals if v])
        d = mean_or_none([v[1] for v in vals if v])
        cons.append(float("nan") if c is None else c)
        plddt.append(float("nan") if d is None else d)
    return np.array(cons), np.array(plddt)


def threshold_sweep(prot, data, width, steps, outdir):
    """Gene and canonical-entry C/N/internal counts over a grid of conservation x pLDDT limits."""
    # every filter except the two swept ones, evaluated once
    base = classify_products(prot, data, {}, width, ignore=("plddt", "conservation"))
    resid = data[4]
    means = {end: window_means(prot, end, width, resid) for end in ("C", "N")}
    canon = (prot.kind == "uniprot").to_numpy()
    rows = []
    for tc in steps:
        for tp in steps:
            ok = {}
            for end in ("C", "N"):
                cons, plddt = means[end]
                # missing data passes; the limit is strict (<), as in the main run
                ok[end] = (base[end + "_pass"].to_numpy()
                           & (np.isnan(cons) | (cons < tc)) & (np.isnan(plddt) | (plddt < tp * 100)))
            frame = pd.DataFrame({"gene": base.gene.to_numpy(), "C": ok["C"], "N": ok["N"]})
            per_gene = frame.groupby("gene")[["C", "N"]].all()
            g_c = int(per_gene.C.sum())
            g_n = int((~per_gene.C & per_gene.N).sum())
            e_c = int(ok["C"][canon].sum())
            e_n = int((~ok["C"] & ok["N"])[canon].sum())
            rows.append({"conservation_lt": round(tc, 1), "plddt_lt": round(tp, 1),
                         "genes_C": g_c, "genes_N": g_n, "genes_internal": len(per_gene) - g_c - g_n,
                         "entries_C": e_c, "entries_N": e_n,
                         "entries_internal": int(canon.sum()) - e_c - e_n})
    sweep = pd.DataFrame(rows)
    sweep.to_csv(outdir / "termini_threshold_sweep.tsv", sep="\t", index=False)
    out = ["", "== threshold sweep (window {} aa; genes; rows = conservation (JSD) < x, columns = "
               "pLDDT/100 < y; missing data passes) ==".format(width)]
    for col, label in (("genes_C", "C-terminal"), ("genes_N", "N-terminal"),
                       ("genes_internal", "internal")):
        out.append("{} genes:".format(label))
        out.append("  cons\\pLDDT " + "".join("{:>7.1f}".format(t) for t in steps))
        for tc in steps:
            sub = sweep[sweep.conservation_lt == round(tc, 1)]
            out.append("  {:>10.1f} ".format(tc) + "".join("{:>7}".format(v) for v in sub[col]))
    return out


def termini(conn, prot, status, cfg, reagent_conn, fp_genes, outdir, width):
    """Classify products and genes at the main window plus a window sensitivity table."""
    block, site, topo = load_features(conn)
    mods = load_modifications(prot, status, cfg["run_dir"])
    resid = load_terminal_residues(conn, max(WINDOWS))
    iso = load_isoform_spans(conn)
    data = (block, site, topo, mods, resid, iso)
    guides = guide_distances(reagent_conn, prot)
    out = ["", "== termini (window {} aa; modifications read for {} proteins) ==".format(
        width, len(mods))]
    prods = classify_products(prot, data, guides, width)
    genes = classify_genes(prods)
    out += class_lines("products, all kinds", prods)
    out += class_lines("products, canonical (kind=uniprot)", prods[prods.kind == "uniprot"])
    out += class_lines("genes", genes)
    out += filter_lines(prods)
    # classes among genes with and without a fluorescent tag already
    genes["fp_tagged"] = genes.gene.isin(fp_genes)
    out.append("genes by class and existing WormTagDB FP tag:")
    for cls in ("C", "N", "internal"):
        sub = genes[genes["class"] == cls]
        out.append("  {:<9}untagged {:>6}   FP-tagged {:>6}".format(
            cls, int((~sub.fp_tagged).sum()), int(sub.fp_tagged.sum())))
    out.append("")
    out.append("window sensitivity (genes: C / N / internal):")
    for w in WINDOWS:
        g = classify_genes(classify_products(prot, data, guides, w))
        c = g["class"].value_counts()
        out.append("  window {:>2}: {} / {} / {}".format(w, c.get("C", 0), c.get("N", 0),
                                                         c.get("internal", 0)))
    g = classify_genes(classify_products(prot, data, guides, width, with_site=True))
    c = g["class"].value_counts()
    out.append("  window {:>2} also blocking on UniProt 'site' features: {} / {} / {}".format(
        width, c.get("C", 0), c.get("N", 0), c.get("internal", 0)))
    out += threshold_sweep(prot, data, width, [i / 10 for i in range(11)], outdir)
    # the same classification without the pLDDT and conservation filters
    skip = ("plddt", "conservation")
    prods_nf = classify_products(prot, data, guides, width, ignore=skip)
    genes_nf = classify_genes(prods_nf)
    genes_nf["fp_tagged"] = genes_nf.gene.isin(fp_genes)
    out += ["", "== termini ignoring pLDDT and conservation (window {} aa) ==".format(width)]
    out += class_lines("products, all kinds", prods_nf)
    out += class_lines("products, canonical (kind=uniprot)", prods_nf[prods_nf.kind == "uniprot"])
    out += class_lines("genes", genes_nf)
    out.append("genes by class and existing WormTagDB FP tag:")
    for cls in ("C", "N", "internal"):
        sub = genes_nf[genes_nf["class"] == cls]
        out.append("  {:<9}untagged {:>6}   FP-tagged {:>6}".format(
            cls, int((~sub.fp_tagged).sum()), int(sub.fp_tagged.sum())))
    prods_nf.to_csv(outdir / "termini_products_no_plddt_cons.tsv", sep="\t", index=False)
    genes_nf.to_csv(outdir / "termini_genes_no_plddt_cons.tsv", sep="\t", index=False)
    prods.to_csv(outdir / "termini_products.tsv", sep="\t", index=False)
    genes.to_csv(outdir / "termini_genes.tsv", sep="\t", index=False)
    return out


def main(argv=None):
    """CLI entry point."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--window", type=int, default=10, help="terminal window in residues")
    ap.add_argument("--db", default=None, help="database (default: batch.config.json path)")
    ap.add_argument("--reagents-db", default=None,
                    help="database holding reagent tables (default: --db); a snapshot of the "
                         "in-progress build file works")
    args = ap.parse_args(argv)
    cfg = bdb.db_config()
    conn = open_db(args.db or cfg["path"])
    reagent_conn = open_db(args.reagents_db) if args.reagents_db else conn
    prot = load_proteins(conn)
    status = bdb.read_last_status(cfg["run_dir"])
    outdir = Path(cfg["run_dir"]) / "analysis"
    outdir.mkdir(exist_ok=True)
    lines = census(conn, prot)
    tag_lines, fp = wormtag(prot)
    lines += tag_lines
    lines += termini(conn, prot, status, cfg, reagent_conn, fp, outdir, args.window)
    # the yeast comparison is computed elsewhere (yeast_tag_rescue_analysis/) and kept as text
    yeast = outdir / "yeast_summary.txt"
    if yeast.exists():
        lines.append(yeast.read_text().rstrip("\n"))
    text = "\n".join(lines)
    (outdir / "summary.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
