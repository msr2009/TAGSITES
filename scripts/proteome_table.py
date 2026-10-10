"""
proteome_table.py

One row per product (canonical entry, isoform record or WormBase-only protein) of
proteome.sqlite3, with every analysis in one TSV for exploration in pandas.

Columns: identity and gene, length and AFDB model, WormTagDB status (tagged any / FP-tagged and
the tag terminus), per-terminus constitutive flag, window-mean conservation and pLDDT at 1, 5,
10 and 20 residues, blocking features at the terminal window (UniProt features, modification
sites and DeepTMHMM signal / TM helix labels, by name), the per-filter checks and C / N /
internal class from proteome_analysis, and the guide distance when reagents are available.

Windowed checks (constitutive, blocking features, class) use --window residues (default 10).

Output: WORMPRO/proteome_isoforms.tsv (--out to change).

Matt Rich, 2026
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import build_proteome_db as bdb  # noqa: E402
import proteome_analysis as pa  # noqa: E402

STRAINS = pa.WORMTAG.parent / "filtered_strains_2026-10-06_v2.csv"
OUT = Path(__file__).parent.parent / "WORMPRO" / "proteome_isoforms.tsv"
TERM = {"C-terminal": "C", "N-terminal": "N", "Internal": "Internal"}


def load_named_modifications(prot, status, run_dir):
    """{pid: [(start, stop, name)]} modification sites read from the run folders."""
    mods = {}
    for r in prot.itertuples():
        if status.get(r.id, {}).get("tasks", {}).get("modifications") != "ok":
            continue
        path = Path(run_dir) / r.id / (r.id + bdb.SUFFIX["modifications"])
        if path.exists() and path.stat().st_size:
            # stop_shift -1: regex_sites writes match.end()+1 (see build_proteome_db)
            mods[r.pid] = [(s, e, d) for _, _, s, e, d in bdb._range_features(
                str(path), "modification", -1)]
    return mods


def tag_table():
    """WormTagDB per-gene frame (gene, wbid, fp, locations of FP-tagged strains)."""
    df = pd.read_csv(STRAINS, keep_default_na=False)
    df["gene"] = df.Gene.str.replace(r"<[^>]+>", "", regex=True).str.strip()
    df["wbid"] = df["WormBase ID"].str.replace(r"<[^>]+>", "", regex=True).str.strip()
    df["is_fp"] = ~df.Fluor.isin(["", "NA"])
    rows = []
    for (gene, wbid), grp in df.groupby(["gene", "wbid"]):
        locs = {TERM.get(t, "Other") for t in grp[grp.is_fp]["Tag Location"]}
        rows.append({"gene": gene, "wbid": wbid, "fp": bool(grp.is_fp.any()),
                     "terminus": ";".join(sorted(locs)) or None})
    return pd.DataFrame(rows)


def wormtag_columns(prot):
    """{gene key: (tagged, fp_tagged, terminus)} for every database gene with a tag record."""
    tag = tag_table()
    hits = pa.match_tags(prot, tag)
    out = {}
    for key, idx in hits.items():
        sub = tag.loc[sorted(idx)]
        locs = sorted({t for v in sub.terminus.dropna() for t in v.split(";")})
        out[key] = (True, bool(sub.fp.any()), ";".join(locs) or None)
    return out


def blocking_names(p, end, width, block, mods, topo):
    """Names of UniProt features, modification sites and signal / TM labels in the window."""
    lo, hi = pa.window_of(end, p.seq_length, width)
    names = ["UniProt: " + d for s, e, d in block.get(p.pid, []) if s <= hi and e >= lo]
    names += ["modification: " + d for s, e, d in mods.get(p.pid, []) if s <= hi and e >= lo]
    names += ["topology: " + d for s, e, d in topo.get(p.pid, [])
              if d in pa.BAD_TOPO and s <= hi and e >= lo]
    return "; ".join(sorted(set(names))) or None


def site_names(p, end, width, site):
    """Names of informational UniProt_site features (modified residue, mutagenesis, ...) in the window."""
    lo, hi = pa.window_of(end, p.seq_length, width)
    names = [d for s, e, d in site.get(p.pid, []) if s <= hi and e >= lo]
    return "; ".join(sorted(set(names))) or None


def build(conn, reagent_conn, cfg, width):
    """The full per-product table."""
    prot = pa.load_proteins(conn)
    status = bdb.read_last_status(cfg["run_dir"])
    block, site, topo = pa.load_features(conn)
    mods_named = load_named_modifications(prot, status, cfg["run_dir"])
    mods = {k: [(s, e) for s, e, _ in v] for k, v in mods_named.items()}
    resid = pa.load_terminal_residues(conn, max(pa.WINDOWS))
    iso = pa.load_isoform_spans(conn)
    data = (block, site, topo, mods, resid, iso)
    guides = pa.guide_distances(reagent_conn, prot)
    cls = pa.classify_products(prot, data, guides, width)
    tags = wormtag_columns(prot)
    n_alt = pd.read_sql("SELECT pid, COUNT(*) AS n_listed_isoforms FROM isoforms "
                        "WHERE is_query = 0 GROUP BY pid", conn)
    df = prot.merge(n_alt, on="pid", how="left").fillna({"n_listed_isoforms": 0})
    df["wormtagdb_tagged"] = [tags.get(g, (False, False, None))[0] for g in df.gene]
    df["wormtagdb_fp"] = [tags.get(g, (False, False, None))[1] for g in df.gene]
    df["wormtagdb_terminus"] = [tags.get(g, (False, False, None))[2] for g in df.gene]
    for end, name in (("N", "nterm"), ("C", "cterm")):
        # a product with no alternative isoform is constitutive by definition
        df[name + "_constitutive"] = [
            True if p.pid not in iso else pa.isoform_ok(iso[p.pid], *pa.window_of(
                end, p.seq_length, width), end) for p in df.itertuples()]
        for w in pa.WINDOWS:
            cons, plddt = pa.window_means(df, end, w, resid)
            df["{}_cons_{}aa".format(name, w)] = cons
            df["{}_plddt_{}aa".format(name, w)] = plddt
        df[name + "_blocking_mods"] = [blocking_names(p, end, width, block, mods_named, topo)
                                       for p in df.itertuples()]
        # informational only: never part of the pass / fail checks
        df[name + "_uniprot_sites"] = [site_names(p, end, width, site) for p in df.itertuples()]
    checks = [c for c in cls.columns if c[:2] in ("C_", "N_")] + ["class", "guide_bp",
                                                                   "guide_status"]
    cls = cls[["pid"] + checks].rename(columns={
        c: ("cterm_" if c[0] == "C" else "nterm_") + c[2:] if c[:2] in ("C_", "N_") else
        "terminus_" + c for c in checks})
    df = df.merge(cls, on="pid")
    df["n_products_in_gene"] = df.groupby("gene").pid.transform("size")
    return df


def main(argv=None):
    """CLI entry point."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--window", type=int, default=10, help="terminal window in residues")
    ap.add_argument("--db", default=None, help="database (default: batch.config.json path)")
    ap.add_argument("--reagents-db", default=None, help="database holding reagent tables")
    ap.add_argument("--out", default=str(OUT), help="output TSV")
    args = ap.parse_args(argv)
    cfg = bdb.db_config()
    conn = pa.open_db(args.db or cfg["path"])
    reagent_conn = pa.open_db(args.reagents_db) if args.reagents_db else conn
    df = build(conn, reagent_conn, cfg, args.window)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print("wrote {} rows x {} columns to {}".format(len(df), df.shape[1], out))


if __name__ == "__main__":
    main()
