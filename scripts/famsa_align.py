"""
famsa_align.py

In-process multiple sequence alignment with pyfamsa (FAMSA), used by
conservation_remote.py in place of an EBI Clustal Omega job; that job stays as
the fallback when pyfamsa is missing or raises. No binary, no data files, so it
runs on shinyapps.io.

Measured against the stored EBI Clustal Omega alignments (6 ortholog sets, 17-100
sequences): JSD conservation Spearman 0.92-0.97 vs Clustal, the same agreement
local MAFFT reaches, in 0.05-3 s.

FAMSA rejects U, O and J, so those residues are aligned as C, K and L and
restored afterwards. Records are aligned under their index, not their header,
because headers may repeat; output keeps input order and headers, wrapped at 60
columns like Clustal's aln-fasta.

Pinned to pyfamsa 0.5.3.post1: its Linux wheel is manylinux2014, while 0.7.0
ships only manylinux_2_34 (glibc 2.34+).

CLI:  python scripts/famsa_align.py in.fasta out.aln
"""

import argparse

# residues FAMSA cannot read -> stand-ins aligned in their place
FAMSA_MASK = {"U": "C", "O": "K", "J": "L"}
LINE_WIDTH = 60


def mask_residues(fasta_text, mapping):
    """Replace mapped residues in sequence lines of a FASTA string; headers are untouched."""
    table = str.maketrans(mapping)
    return "\n".join(ln if ln.startswith(">") else ln.translate(table)
                     for ln in fasta_text.split("\n"))


def restore_residues(aln_path, fasta_str_list, mapping):
    """Put masked residues back into an alignment file, matching records by input order."""
    # ungapped residue index -> original letter, per record
    sites = [{i: c for i, c in enumerate(rec.split("\n", 1)[1].replace("\n", ""))
              if c in mapping}
             for rec in fasta_str_list]
    if not any(sites):
        return
    with open(aln_path) as f:
        lines = f.read().split("\n")
    out, rec_idx, resi = [], -1, 0
    for ln in lines:
        if ln.startswith(">"):
            rec_idx, resi = rec_idx + 1, 0
            out.append(ln)
            continue
        chars = []
        for c in ln:
            # gaps do not advance the residue index; only real residues are checked
            if c != "-":
                if rec_idx < len(sites) and resi in sites[rec_idx]:
                    c = sites[rec_idx][resi]
                resi += 1
            chars.append(c)
        out.append("".join(chars))
    with open(aln_path, "w") as f:
        f.write("\n".join(out))


def _parse_records(fasta_str_list):
    """(header, sequence) per FASTA record string."""
    recs = []
    for rec in fasta_str_list:
        head, _, body = rec.lstrip().partition("\n")
        recs.append((head.lstrip(">").strip(), body.replace("\n", "").strip()))
    return recs


def align(fasta_str_list, aln_path):
    """Align FASTA record strings with pyfamsa and write a wrapped FASTA alignment."""
    from pyfamsa import Aligner, Sequence

    recs = _parse_records([mask_residues(r, FAMSA_MASK) for r in fasta_str_list])
    seqs = [Sequence(str(i).encode(), s.encode()) for i, (_, s) in enumerate(recs)]
    msa = Aligner(guide_tree="sl", keep_duplicates=True).align(seqs)
    by_idx = {int(s.id.decode()): s.sequence.decode() for s in msa}
    if len(by_idx) != len(recs):
        raise ValueError(f"FAMSA returned {len(by_idx)} of {len(recs)} sequences")

    # write in input order under the original headers
    with open(aln_path, "w") as f:
        for i, (head, _) in enumerate(recs):
            row = by_idx[i]
            f.write(f">{head}\n")
            f.writelines(row[j:j + LINE_WIDTH] + "\n" for j in range(0, len(row), LINE_WIDTH))
    restore_residues(aln_path, fasta_str_list, FAMSA_MASK)


def _read_records(path):
    """FASTA file -> list of record strings."""
    with open(path) as f:
        text = f.read()
    return [">" + r for r in text.split(">")[1:]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align a protein FASTA with pyfamsa.")
    parser.add_argument("fasta_in", help="input FASTA")
    parser.add_argument("aln_out", help="output FASTA alignment")
    args = parser.parse_args()
    align(_read_records(args.fasta_in), args.aln_out)
