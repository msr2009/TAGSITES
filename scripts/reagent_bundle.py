"""
reagent_bundle.py

Assembles the reagent download for chosen sites and guides, with no Shiny dependency: the
homology-arm, guide, oligo, SapTrap and genotyping TSVs plus the annotated GenBank records.
Shared by the Reagents tab (modules/reagents_server.py) and the batch app (app_batch.py), so
both produce the same files for the same choices.

Naming goes through a `names` dict of callables so the single-protein app keeps its
{run}_{site}{strand}{distance} names while the batch app uses {gene}_{uniprot}_{site}_{tag}
(see single_names / batch_names).

    python scripts/reagent_bundle.py --run_json data/x/x.run.json --sites 45,120 \\
        --insert FLAG --repair_type ssodn --out x_reagents.zip

Matt Rich, 2026
"""

import argparse
import io
import json
import os
import sys
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from genbank_export import (  # noqa: E402
    PRIMER_TYPES as GENBANK_PRIMER_TYPES,
    build_gdna_record,
    build_knockin_record,
    load_cds_strand,
    load_exons,
    load_region,
)
from plasmid_assembly import OLIGO_MAX, build_saptrap_site  # noqa: E402
from reagent_sequences import (  # noqa: E402
    build_ssodn,
    design_pcr_primers,
    load_tags,
    reconstruct_wt_arms,
    truncate_arms,
    truncate_arms_with_tolerance,
)

TAGS_PATH = Path(__file__).parent.parent / "tables" / "tags.tsv"

# arm-length defaults by repair type (the Reagents tab's defaults)
ARM_DEFAULTS = {"ha": 1000, "ssodn": 45, "amplicon": 45, "plasmid": 57}

# defaults for the repair-specific options the Reagents tab exposes as inputs
OPTION_DEFAULTS = {"pcr_tm": 60.0, "pcr_phos": False, "pcr_fwd_ext": "", "pcr_rev_ext": "",
                   "saptrap_tm": 60.0, "saptrap_tolerance": 0.2, "saptrap_genbank": False}


def load_tag_table(path=TAGS_PATH):
    """Insert library {name: dna} from tables/tags.tsv ({} when absent)."""
    return load_tags(path) if Path(path).exists() else {}


def add_guide_ids(df):
    """Add the guide_id column ({strand}_{pam_fwd_start}_{spacer}) to a reagents frame."""
    df["guide_id"] = (df["guide_strand"].astype(str) + "_" + df["pam_fwd_start"].astype(str)
                      + "_" + df["spacer"].astype(str))
    return df


def load_reagents_tsv(path):
    """Read a {run}.reagents.tsv with guide_id added."""
    return add_guide_ids(pd.read_csv(path, sep="\t"))


def default_guides(df, sites):
    """{residue_index: guide_id} picking the guide nearest the insertion for each site."""
    out = {}
    if df is None:
        return out
    sub = df[df["residue_index"].isin([int(s) for s in sites])]
    # rows sort by distance, so the first per residue is the nearest guide
    for rid, grp in sub.sort_values("distance", kind="stable").groupby("residue_index"):
        out[int(rid)] = str(grp.iloc[0]["guide_id"])
    return out


def wt_arms_for_row(row):
    """Return (wt_left, wt_right, wt_is_true) for a reagents row.

    Prefers the genuine left_arm_wt/right_arm_wt TSV columns; falls back to
    reconstructing WT arms from the recorded mutation for older TSVs. These
    feed truncate_arms/truncate_arms_with_tolerance so a recut-blocking
    mutation farther than arm_length from the junction still gets included
    rather than truncated away. wt_is_true is False when reconstruction
    failed and only the mutated arms are available.
    """
    if "left_arm_wt" in row and "right_arm_wt" in row:
        return str(row["left_arm_wt"]), str(row["right_arm_wt"]), True
    recon = reconstruct_wt_arms(
        str(row["left_arm"]), str(row["right_arm"]),
        int(row["insert_pos"]),
        str(row.get("recut_block_method", "") or ""),
        str(row.get("mutation_desc", "") or ""),
    )
    if recon is not None:
        return recon[0], recon[1], True
    return str(row["left_arm"]), str(row["right_arm"]), False


def select_rows(df, selection, arm_len, repair, tolerance=0.2):
    """List of (row, left_arm, right_arm) for each selected guide.

    selection is {residue_index: iterable of guide_id}. For the SapTrap long-arm
    (PCR-primer) path, arms are widened on the far/outer side per tolerance so
    sapi_amplification_primers has room to search; every other repair type gets the
    plain exact truncation.
    """
    out = []
    if df is None or not selection:
        return out
    use_widened = repair == "plasmid" and arm_len > OLIGO_MAX
    for rid, gid_set in selection.items():
        for gid in gid_set:
            matches = df[(df["residue_index"] == rid) & (df["guide_id"] == gid)]
            if matches.empty:
                continue
            row = matches.iloc[0]
            wt_left_src, wt_right_src, wt_is_true = wt_arms_for_row(row)
            wt_kwargs = dict(
                left_arm_wt=wt_left_src if wt_is_true else None,
                right_arm_wt=wt_right_src if wt_is_true else None,
            )
            try:
                if use_widened:
                    left, right = truncate_arms_with_tolerance(
                        str(row["left_arm"]), str(row["right_arm"]), arm_len, tolerance,
                        **wt_kwargs,
                    )
                else:
                    left, right = truncate_arms(
                        str(row["left_arm"]), str(row["right_arm"]), arm_len, **wt_kwargs,
                    )
            except ValueError:
                continue
            out.append((row, left, right))
    return out


def guide_label(row):
    """Compact guide label: {site}{strand}{distance} e.g. K45+3."""
    return "{}{}{}{}".format(
        row["amino_acid"], int(row["residue_index"]),
        str(row["guide_strand"]), int(row["distance"]),
    )


def site_name(row):
    """Site label {aa}{residue}, e.g. K45."""
    return "{}{}".format(row["amino_acid"], int(row["residue_index"]))


def tag_label(name):
    """Tag name for use in reagent names; 'custom' for a custom insert."""
    return "custom" if not name or name == "Custom" else name.replace(" ", "-")


def single_names(rn, tag):
    """Naming of the single-protein Reagents tab: {rn}_{site} arms, {rn}_{K45+3} guides."""
    return {
        "zip": rn,
        "record": rn,
        "arm": lambda row: "{}_{}".format(rn, site_name(row)),
        "guide": lambda row: "{}_{}".format(rn, guide_label(row)),
        "oligo": lambda row: "{}_{}_{}".format(rn, guide_label(row), tag),
        "saptrap": lambda row, name: name,
        "saptrap_gb": lambda row, rec: "{}_{}_saptrap.gb".format(rn, rec.id),
        "genbank": lambda row: "{}_{}".format(rn, guide_label(row)),
    }


def batch_names(prefix, tag, ranks=None):
    """Batch naming: every reagent of a site is {prefix}_{site}_{tag}, e.g. trxr-1_Q17745_K45_FLAG.

    With `ranks` ({residue_index: rank}) the site's rank joins the name as
    trxr-1_Q17745_K45_r1_FLAG. The batch's output is reagents rather than the site itself,
    so the rank says which of several equally-scoring sites this one was — r1 is the app's
    first choice. It is reproducible for the same scores (see
    utils.scoring.pick_suggested_regions); a hand-added site has no rank and is left bare.

    One guide per site in the batch app, so strand and distance stay out of the name (they
    remain columns of the guides TSV). SapTrap oligos, bare {site}_sgRNA_top from
    plasmid_assembly, get the same base so a combined file has no collisions.
    """
    def base(row):
        rank = (ranks or {}).get(int(row["residue_index"]))
        site = site_name(row) if rank is None else "{}_r{}".format(site_name(row), rank)
        return "{}_{}_{}".format(prefix, site, tag)

    return {
        "zip": "{}_{}".format(prefix, tag),
        "record": prefix,
        "arm": base,
        "guide": base,
        "oligo": base,
        "saptrap": lambda row, name: base(row) + name[len(site_name(row)):],
        "saptrap_gb": lambda row, rec: "{}_saptrap.gb".format(base(row)),
        "genbank": base,
    }


def ha_tsv(rows, names):
    """[name]_homology-arms.tsv: name, 5p_HA, 3p_HA, guide."""
    lines = ["name\t5p_HA\t3p_HA\tguide"]
    for row, left, right in rows:
        lines.append("{}\t{}\t{}\t{}".format(names["arm"](row), left, right, str(row["spacer"])))
    return "\n".join(lines) + "\n"


def guides_tsv(rows, names):
    """[name]_guides.tsv: name, guide sequence."""
    lines = ["name\tsequence"]
    for row, left, right in rows:
        lines.append("{}\t{}".format(names["guide"](row), str(row["spacer"])))
    return "\n".join(lines) + "\n"


def oligos_tsv(rows, names, repair, insert, opts):
    """[name]_oligos.tsv for ssodn/amplicon; None when not applicable."""
    if repair not in ("ssodn", "amplicon"):
        return None
    lines = ["name\tsequence\tnotes"]
    for row, left, right in rows:
        base = names["oligo"](row)
        cut = int(row["cut_pos"])
        ins = int(row["insert_pos"])
        if repair == "ssodn":
            seq, _, _ = build_ssodn(left, right, insert, cut, ins, True)
            lines.append("{}_repair\t{}\t".format(base, seq))
        elif repair == "amplicon" and insert:
            (_, fwd), (_, rev), meta = design_pcr_primers(
                left, right, insert, opts["pcr_tm"], opts["pcr_phos"], cut, ins
            )
            if opts["pcr_fwd_ext"]:
                fwd = opts["pcr_fwd_ext"] + fwd
            if opts["pcr_rev_ext"]:
                rev = opts["pcr_rev_ext"] + rev
            notes = "fallback=True" if meta["used_fallback"] else ""
            lines.append("{}_F\t{}\t{}".format(base, fwd, notes))
            lines.append("{}_R\t{}\t{}".format(base, rev, notes))
    return "\n".join(lines) + "\n" if len(lines) > 1 else None


def saptrap_tsv(rows, names, arm_len, opts):
    """[name]_saptrap-oligos.tsv plus (filename, SeqRecord) pairs when GenBank is on."""
    lines = ["name\tsequence\tkind\tnotes"]
    gb_pairs = []
    for row, left, right in rows:
        oligo_rows, records = build_saptrap_site(
            row, left, right, tm_target=opts["saptrap_tm"], genbank=opts["saptrap_genbank"],
            arm_length=arm_len, tolerance=opts["saptrap_tolerance"],
        )
        for orow in oligo_rows:
            lines.append("{}\t{}\t{}\t{}".format(
                names["saptrap"](row, orow["name"]), orow["sequence"],
                orow.get("kind", ""), orow.get("notes", "")))
        for rec in records:
            gb_pairs.append((names["saptrap_gb"](row, rec), rec))
    tsv_text = "\n".join(lines) + "\n" if len(lines) > 1 else None
    return tsv_text, gb_pairs


def genotyping_tsv(results, df):
    """[name]_genotyping-primers.tsv from {residue: {amplicon_type: primer dict}}, or None."""
    if not results:
        return None
    lines = ["residue_index\tamino_acid\tamplicon_type\tfwd_seq\tfwd_tm\trev_seq\trev_tm"
             "\tproduct_size"]
    for rid in sorted(results):
        primers = results[rid]
        if not primers:
            continue
        aa = ""
        if df is not None:
            match = df[df["residue_index"] == rid]
            if not match.empty:
                aa = str(match.iloc[0]["amino_acid"])
        for amplicon_type, p in primers.items():
            lines.append("{}\t{}\t{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}".format(
                rid, aa, amplicon_type,
                p["fwd_seq"], p["fwd_tm"], p["rev_seq"], p["rev_tm"],
                int(p["product_size"]),
            ))
    return "\n".join(lines) + "\n" if len(lines) > 1 else None


def load_genotyping(path, sites):
    """Precomputed {residue: {amplicon_type: primer dict}} from a .reagents.genotyping.tsv."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    gt = pd.read_csv(path, sep="\t")
    gt = gt[gt["residue_index"].isin([int(s) for s in sites])]
    # one primer pair per (residue, amplicon type)
    for r in gt.itertuples():
        out.setdefault(int(r.residue_index), {})[str(r.amplicon_type)] = {
            "fwd_seq": r.fwd_seq, "fwd_tm": float(r.fwd_tm), "rev_seq": r.rev_seq,
            "rev_tm": float(r.rev_tm), "product_size": int(r.product_size)}
    return out


def genbank_records(rows, names, region_fa, genewise, insert, tag, primers_by_site):
    """Two GenBank records per selected guide (WT region and knock-in) as (filename, record).

    The knock-in record is per-guide rather than per-site because the recut-blocking
    edits it annotates are guide-specific. Returns [] when the region FASTA is unavailable.
    """
    region = load_region(region_fa)
    if region is None:
        return []
    exons = load_exons(genewise)
    strand = load_cds_strand(genewise)
    rn = names["record"]
    pairs = []
    for row, left, right in rows:
        label = names["genbank"](row)
        primers = (primers_by_site or {}).get(int(row["residue_index"])) or {}
        wt_left, wt_right, wt_is_true = wt_arms_for_row(row)
        try:
            gdna = build_gdna_record(
                row, region, wt_left, wt_right, exons, primers,
                GENBANK_PRIMER_TYPES, run_name=rn, cds_strand=strand,
            )
            knockin = build_knockin_record(
                row, region, left, right, insert,
                left_wt=wt_left if wt_is_true else None,
                right_wt=wt_right if wt_is_true else None,
                exons=exons, primers=primers, types=GENBANK_PRIMER_TYPES,
                insert_name=tag, run_name=rn, cds_strand=strand,
            )
        except Exception:
            continue
        pairs.append(("{}_gDNA.gb".format(label), gdna))
        pairs.append(("{}_knockin.gb".format(label), knockin))
    return pairs


def _gb_text(rec):
    """A SeqRecord as GenBank text."""
    from Bio import SeqIO
    buf = io.StringIO()
    SeqIO.write(rec, buf, "genbank")
    return buf.getvalue()


def build_files(rows, names, repair, insert, tag, arm_len, opts=None, genotyping=None,
                df=None, region_fa="", genewise=""):
    """Ordered [(filename, text)] of the reagent download, as the Reagents tab ZIP holds it."""
    opts = {**OPTION_DEFAULTS, **(opts or {})}
    zn = names["zip"]
    files = [("{}_homology-arms.tsv".format(zn), ha_tsv(rows, names)),
             ("{}_guides.tsv".format(zn), guides_tsv(rows, names))]
    geno = genotyping_tsv(genotyping, df)
    if geno:
        files.append(("{}_genotyping-primers.tsv".format(zn), geno))
    # plasmid mode gets SapTrap oligos; ssODN/amplicon get repair oligos
    if repair == "plasmid":
        tsv_text, gb_pairs = saptrap_tsv(rows, names, arm_len, opts)
        if tsv_text:
            files.append(("{}_saptrap-oligos.tsv".format(zn), tsv_text))
        files += [(fname, _gb_text(rec)) for fname, rec in gb_pairs]
    else:
        oligos = oligos_tsv(rows, names, repair, insert, opts)
        if oligos:
            files.append(("{}_oligos.tsv".format(zn), oligos))
    # annotated GenBank records always ride along with the TSVs
    files += [(fname, _gb_text(rec)) for fname, rec in
              genbank_records(rows, names, region_fa, genewise, insert, tag, genotyping)]
    return files


def write_zip(files, dest=None):
    """Write [(name, text)] to a ZIP at dest (path or buffer); return bytes when dest is None."""
    buf = dest if dest is not None else io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, text in files:
            zf.writestr(name, text)
    if dest is None:
        buf.seek(0)
        return buf.read()
    return None


def run_paths(run_json):
    """(run_name, reagents TSV, region FASTA, Genewise output, genotyping TSV) of a run JSON."""
    from task_registry import companion_path
    with open(run_json) as fh:
        j = json.load(fh)
    gb = j.get("global", {})
    rn, wd = gb.get("run_name", ""), gb.get("working_dir", "")
    args = j.get("tasks", {}).get("REAGENTS_reagents", {}).get("args", {})
    tsv = args.get("output", "")
    # the runner's own *_genewise.* files win over the run JSON's inputs (as in the app)
    prefix = os.path.join(wd, rn + "_genewise") if wd and rn else ""
    fa = [prefix + ".genewise_genomic.fa"] if prefix else []
    gw = [prefix + ".genewise.out.txt"] if prefix else []
    fa.append(args.get("genomic_fasta", ""))
    gw.append(args.get("genewise", ""))
    region_fa = next((p for p in fa if p and os.path.exists(p)), "")
    genewise = next((p for p in gw if p and os.path.exists(p)), "")
    geno = companion_path(tsv, "reagents", "genotyping") if tsv else ""
    return rn, tsv, region_fa, genewise, geno


def main(run_json, sites, insert="", repair_type="ha", arm_length=None, out=None,
         prefix=None, opts=None):
    """Build the reagent ZIP for chosen sites of a finished run; nearest guide per site."""
    tags = load_tag_table()
    rn, tsv, region_fa, genewise, geno_path = run_paths(run_json)
    df = load_reagents_tsv(tsv)
    insert_name = insert if insert in tags else ("Custom" if insert else "")
    insert_seq = tags[insert].upper() if insert in tags else insert.upper()
    tag = tag_label(insert_name)
    arm_len = arm_length or ARM_DEFAULTS.get(repair_type, 500)
    selection = {r: {g} for r, g in default_guides(df, sites).items()}
    rows = select_rows(df, selection, arm_len, repair_type,
                       (opts or {}).get("saptrap_tolerance", 0.2))
    names = batch_names(prefix, tag) if prefix else single_names(rn or "run", tag)
    files = build_files(rows, names, repair_type, insert_seq, tag, arm_len, opts,
                        genotyping=load_genotyping(geno_path, sites), df=df,
                        region_fa=region_fa, genewise=genewise)
    out = out or "{}_reagents.zip".format(names["zip"])
    write_zip(files, out)
    print("wrote {} ({} files, {} guides)".format(out, len(files), len(rows)))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--run_json", required=True)
    ap.add_argument("--sites", required=True, help="comma-separated residue numbers")
    ap.add_argument("--insert", default="", help="tags.tsv name, or a raw DNA sequence")
    ap.add_argument("--repair_type", default="ha", choices=sorted(ARM_DEFAULTS))
    ap.add_argument("--arm_length", type=int, default=None)
    ap.add_argument("--prefix", default=None, help="batch naming prefix, e.g. trxr-1_Q17745")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    main(a.run_json, [int(s) for s in a.sites.split(",") if s.strip()], a.insert,
         a.repair_type, a.arm_length, a.out, a.prefix)
