"""
batch_inputs.py

Turns a pasted or uploaded protein list into batch entries resolved against the proteome
database. A list is one term per line (gene name, UniProt accession or isoform id, WBGene id,
WormBase transcript), optionally tab/comma separated with an insert name in the second column
(a header row naming an `insert` column is also accepted). FASTA records are matched by exact
sequence instead.

A gene can resolve to several products (UniProt canonical rows, UniProt isoforms, WormBase-only
transcripts; trxr-1 has two canonical rows), so each entry keeps every product as a choice and
preselects one: a product whose id equals the term, else kind uniprot over the others, then
reagents designed, then a structure, then the longest.

    python scripts/batch_inputs.py genes.txt

Matt Rich, 2026
"""

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import proteome_db as P  # noqa: E402

_KIND_RANK = {"uniprot": 0, "uniprot_isoform": 1, "wormbase_only": 2}


def parse_list(text, tag_names=()):
    """[{input, insert, sequence}] from list text (terms, optional insert column, or FASTA).

    Cells split on tabs, commas and spaces. A cell is an insert when a header row names an
    `insert` column, or when it is the second cell of a line and names a known tag (or is raw
    DNA); every other cell is a protein term, so "trxr-1, unc-54" on one line is two proteins.
    """
    text = (text or "").strip()
    if not text:
        return []
    if text.startswith(">"):
        return _parse_fasta(text)
    tags = {t.lower() for t in tag_names}
    lines = [re.split(r"[\t,; ]+", ln.strip()) for ln in text.splitlines()]
    lines = [[c for c in ln if c] for ln in lines]
    insert_col = None
    # a header row is recognized by an 'insert' column name
    if lines and any(c.lower() == "insert" for c in lines[0]):
        insert_col = [c.lower() for c in lines[0]].index("insert")
        lines = lines[1:]
    out, seen = [], set()
    for cells in lines:
        if not cells or cells[0].startswith("#"):
            continue
        ins = ""
        if insert_col is not None:
            ins = cells[insert_col] if insert_col < len(cells) else ""
            terms = [c for i, c in enumerate(cells) if i != insert_col][:1]
        elif len(cells) == 2 and (cells[1].lower() in tags or _is_dna(cells[1])):
            terms, ins = [cells[0]], cells[1]
        else:
            terms = cells
        # one entry per term; a repeated term keeps its first insert
        for term in terms:
            if term.lower() in seen:
                continue
            seen.add(term.lower())
            out.append({"input": term, "insert": ins, "sequence": ""})
    return out


def _is_dna(cell):
    """True for a raw DNA insert (>= 12 nt of ACGT)."""
    return len(cell) >= 12 and set(cell.upper()) <= set("ACGT")


def _parse_fasta(text):
    """FASTA records as entries keyed by their header's first word."""
    out, name, seq = [], None, []
    for line in text.splitlines() + [">"]:
        if line.startswith(">"):
            if name is not None and seq:
                out.append({"input": name, "insert": "", "sequence": "".join(seq)})
            name, seq = (line[1:].split() or ["seq{}".format(len(out) + 1)])[0], []
        elif line.strip():
            seq.append(line.strip())
    return out


def _choice_rows(pids, conn):
    """Product rows (with reagent status) for a list of pids."""
    if not pids:
        return []
    marks = ",".join("?" * len(pids))
    df = P.query("SELECT p.pid, p.id, p.kind, p.accession, p.gene_name, p.wb_transcript, "
                 "p.seq_length, p.has_structure, p.sequence_differs, r.status AS reagent_status "
                 "FROM proteins p LEFT JOIN reagent_status r USING (pid) "
                 "WHERE p.pid IN ({})".format(marks), list(pids), conn)
    return df.to_dict("records")


def _rank(row, term):
    """Sort key: exact id match, kind, reagents designed, structure, length (longest first)."""
    return (0 if str(row["id"]).lower() == term.lower() else 1,
            _KIND_RANK.get(row["kind"], 9),
            0 if row.get("reagent_status") == "success" else 1,
            0 if row.get("has_structure") else 1,
            -int(row.get("seq_length") or 0))


def choice_label(row):
    """Dropdown text for one product: id · kind · length · reagents · structure."""
    reag = "reagents ok" if row.get("reagent_status") == "success" else "no reagents"
    struct = "" if row.get("has_structure") else " · no structure"
    return "{} · {} · {} aa · {}{}".format(row["id"], row["kind"].replace("_", " "),
                                          row.get("seq_length"), reag, struct)


def _clean(row):
    """A product row with JSON-safe plain values."""
    return {k: (None if v != v else (v.item() if hasattr(v, "item") else v))
            for k, v in row.items()}


def resolve(entries, conn=None):
    """Batch protein entries: every matching product as a choice, the preferred one selected."""
    own = conn is None
    conn = conn or P.open_db()
    try:
        out = []
        for e in entries:
            if e.get("sequence"):
                hits = P.match_sequence(e["sequence"], conn=conn)
                pids = hits["pid"].tolist() if hits is not None and len(hits) else []
            else:
                hits = P.find_protein(e["input"], conn=conn)
                pids = hits["pid"].tolist() if hits is not None and len(hits) else []
            choices = sorted((_clean(r) for r in _choice_rows(pids, conn)),
                             key=lambda r: _rank(r, e["input"]))
            out.append(entry_for(e, choices, choices[0] if choices else None))
        mark_duplicates(out)
        return out
    finally:
        if own:
            conn.close()


def entry_for(e, choices, chosen):
    """The batch protein entry for an input with its chosen product (None = not in database)."""
    base = {"input": e["input"], "insert": e.get("insert", ""), "sequence": e.get("sequence", ""),
            "choices": choices}
    if chosen is None:
        return {**base, "key": e["input"], "pid": None, "product_id": None, "gene": None,
                "kind": None, "source": "missing", "reagent_status": None}
    return {**base, "key": chosen["id"], "pid": int(chosen["pid"]), "product_id": chosen["id"],
            "gene": chosen.get("gene_name"), "kind": chosen["kind"], "source": "db",
            "reagent_status": chosen.get("reagent_status")}


def mark_duplicates(entries):
    """Set duplicate_of on entries whose product an earlier entry already uses."""
    first = {}
    for e in entries:
        e["duplicate_of"] = None
        if e["pid"] is None:
            continue
        if e["key"] in first:
            e["duplicate_of"] = first[e["key"]]
        else:
            first[e["key"]] = e["input"]
    return entries


def choose(entry, product_id):
    """Re-point an entry at another of its choices (the Setup isoform dropdown)."""
    chosen = next((c for c in entry["choices"] if c["id"] == product_id), None)
    if chosen is None:
        return entry
    e = {"input": entry["input"], "insert": entry.get("insert", ""),
         "sequence": entry.get("sequence", "")}
    return entry_for(e, entry["choices"], chosen)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Resolve a protein list against the proteome DB.")
    ap.add_argument("list_file")
    a = ap.parse_args()
    from reagent_bundle import load_tag_table
    for ent in resolve(parse_list(Path(a.list_file).read_text(), load_tag_table())):
        alts = len(ent["choices"])
        print("{:<16} -> {:<16} {:<16} {:<10} {} choice(s){}".format(
            ent["input"], str(ent["product_id"]), str(ent["gene"]), ent["source"], alts,
            ("  insert=" + ent["insert"] if ent["insert"] else "")
            + ("  DUPLICATE of " + ent["duplicate_of"] if ent["duplicate_of"] else "")))
