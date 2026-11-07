"""
batch_output.py

The batch app's Output Reagents: one ZIP holding a folder per accepted protein (the same files
the single-protein Reagents tab writes, named {gene}_{uniprot}_{site}_{tag}) plus combined
TSVs across the batch, each with a leading `protein` column, and a summary.

Arms are regenerated here, only for the accepted sites' chosen guides (rows_with_arms): the
batch's run folders hold arm-less reagent TSVs.

    python scripts/batch_output.py <batch.json> [--out batch_reagents.zip]

Matt Rich, 2026
"""

import argparse
import io
import os
import sys
import time
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import batch_state as BS  # noqa: E402
import reagent_bundle as RB  # noqa: E402

# per-protein files that are merged into one batch-wide file (by name suffix)
COMBINED = ("_homology-arms.tsv", "_guides.tsv", "_oligos.tsv", "_saptrap-oligos.tsv",
            "_genotyping-primers.tsv")


def protein_inserts(protein, settings, tags):
    """[(tag label, DNA)] for a protein: its own list insert, else the batch's inserts."""
    names = [protein["insert"]] if protein.get("insert") else list(settings.get("inserts") or [])
    custom = settings.get("custom_insert", "")
    out = []
    for name in names:
        # a tags.tsv name, the batch's custom insert, or raw DNA from the list column
        if name in tags:
            out.append((RB.tag_label(name), tags[name].upper()))
        elif name == "Custom":
            out.append(("custom", custom.strip().upper()))
        else:
            out.append(("custom", name.strip().upper()))
    return out or [("none", "")]


def protein_rows(protein):
    """Reagent frame (with arms) for a protein's chosen guides, plus its run's companion paths."""
    rn, tsv, region_fa, genewise, geno = RB.run_paths(protein["run_json"])
    guides = {int(r): g for r, g in (protein.get("guides") or {}).items()}
    if protein.get("source") == "db":
        import proteome_export as X
        df = RB.add_guide_ids(X.rows_with_arms(protein["pid"], list(guides.items())))
    else:
        df = RB.load_reagents_tsv(tsv)
    return df, guides, {"tsv": tsv, "region_fa": region_fa, "genewise": genewise, "geno": geno}


def _model_check(tsv):
    """Gene-model check status and message for the summary ('' when absent)."""
    from cds_check import read_sidecar
    chk = read_sidecar(str(Path(tsv).with_suffix("")) + ".model_check.json") if tsv else None
    return (chk or {}).get("status", ""), (chk or {}).get("message", "") or ""


def protein_files(protein, settings, tags):
    """[(filename, text)] for one accepted protein across its inserts, and its summary row."""
    df, guides, paths = protein_rows(protein)
    repair = settings.get("repair_type", "ha")
    arm_len = int(settings.get("arm_length") or RB.ARM_DEFAULTS.get(repair, 500))
    sites = sorted(guides)
    selection = {r: {g} for r, g in guides.items()}
    geno = RB.load_genotyping(paths["geno"], sites)
    prefix = BS.name_prefix(protein)
    # rank of the region each site came from, recorded by Review; empty for older batches
    ranks = {int(r): int(v) for r, v in (protein.get("site_ranks") or {}).items()}
    files = []
    inserts = protein_inserts(protein, settings, tags)
    for tag, dna in inserts:
        rows = RB.select_rows(df, selection, arm_len, repair,
                              settings.get("options", {}).get("saptrap_tolerance", 0.2))
        files += RB.build_files(rows, RB.batch_names(prefix, tag, ranks), repair, dna, tag, arm_len,
                                settings.get("options"), genotyping=geno, df=df,
                                region_fa=paths["region_fa"], genewise=paths["genewise"])
    status, message = _model_check(paths["tsv"])
    site_labels = []
    for r in sites:
        hit = df[df["residue_index"] == r]
        site_labels.append("{}{}".format(hit.iloc[0]["amino_acid"], r) if len(hit) else str(r))
    summary = {"protein": prefix, "input": protein["input"], "product": protein["product_id"],
               "gene": protein.get("gene") or "", "sites": ",".join(site_labels),
               "inserts": ",".join(t for t, _ in inserts), "repair_type": repair,
               "arm_length": arm_len, "gene_model": status, "gene_model_note": message,
               "accepted_at": protein.get("accepted_at") or ""}
    return files, summary


def _merge(texts_by_protein):
    """One TSV from per-protein TSV texts: header once, a leading protein column."""
    header, lines = None, []
    for prefix, text in texts_by_protein:
        rows = text.rstrip("\n").split("\n")
        if header is None:
            header = "protein\t" + rows[0]
        lines += ["{}\t{}".format(prefix, r) for r in rows[1:]]
    return "\n".join([header] + lines) + "\n" if header else ""


def build_batch_zip(batch, dest):
    """Write the batch ZIP to dest (path or buffer); return (n proteins, n files, seconds)."""
    t0 = time.time()
    tags = RB.load_tag_table()
    settings = batch.get("settings", {})
    name = batch["name"]
    combined = {suffix: [] for suffix in COMBINED}
    summaries, n_files, errors = [], 0, []
    with zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in batch["proteins"]:
            if p.get("status") != "accepted":
                continue
            prefix = BS.name_prefix(p)
            try:
                files, summary = protein_files(p, settings, tags)
            except Exception as e:   # one bad protein is reported, not fatal to the batch
                errors.append("{}\t{}: {}".format(prefix, type(e).__name__, e))
                continue
            # several inserts give several files of one kind; merge them per protein
            per_kind = {}
            for fname, text in files:
                suffix = next((s for s in COMBINED if fname.endswith(s)), None)
                if suffix:
                    per_kind.setdefault(suffix, []).append(text)
                else:
                    zf.writestr("{}/{}".format(prefix, fname), text)
                    n_files += 1
            for suffix, texts in per_kind.items():
                # insert-independent rows (genotyping primers) repeat per insert; keep one
                body = [ln for t in texts for ln in t.rstrip("\n").split("\n")[1:]]
                merged = "\n".join([texts[0].split("\n", 1)[0]] + list(dict.fromkeys(body)))
                merged += "\n"
                zf.writestr("{}/{}{}".format(prefix, prefix, suffix), merged)
                combined[suffix].append((prefix, merged))
                n_files += 1
            summaries.append(summary)
        for suffix, texts in combined.items():
            if texts:
                zf.writestr("{}{}".format(name, suffix), _merge(texts))
        if summaries:
            cols = list(summaries[0])
            body = ["\t".join(cols)] + ["\t".join(str(s[c]).replace("\t", " ") for c in cols)
                                        for s in summaries]
            zf.writestr("{}_summary.tsv".format(name), "\n".join(body) + "\n")
        if errors:
            zf.writestr("{}_errors.tsv".format(name), "protein\terror\n" + "\n".join(errors) + "\n")
    return len(summaries), n_files, time.time() - t0, errors


def zip_bytes(batch):
    """The batch ZIP as bytes."""
    buf = io.BytesIO()
    build_batch_zip(batch, buf)
    return buf.getvalue()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Write the reagent ZIP of a batch's accepted proteins.")
    ap.add_argument("batch_json")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    b = BS.load_batch(a.batch_json)
    out = a.out or os.path.join(os.path.dirname(a.batch_json), b["name"] + "_reagents.zip")
    n, nf, secs, errs = build_batch_zip(b, out)
    print("wrote {}: {} protein(s), {} per-protein files, {:.1f} s".format(out, n, nf, secs))
    for e in errs:
        print("  error:", e)
