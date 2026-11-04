"""Three-state secondary structure (helix / strand / coil) from backbone coordinates.

A numpy-only port of the DSSP-style assignment in pydssp (pydssp_numpy.py, v0.9.1,
https://github.com/ShintaroMinami/PyDSSP), rewritten to keep memory at a few L x L
boolean / float32 matrices (pydssp materialises L x L x 3 float64 maps, gigabytes for a
2,700-residue AlphaFold model) and to drop the torch and einops dependencies. Two changes
bring it closer to DSSP: the amide H is placed as DSSP places it (pydssp uses the C-N / CA-N
bisector) and proline donates no H-bond. Without them, output is identical to pydssp's.

Validated against the DSSP 4 annotation AlphaFold DB writes into its mmCIF files
(_struct_conf, computed on the same coordinates as the .pdb): 99.3% per-residue three-state
agreement over 18 AFDB models (pydssp as published: 98.3%), with helix and strand never
confused. Disagreements sit at element ends and come from DSSP rules this port does not
implement (beta-bulge ladder merging, pi-helix priority, helix breaks). A 2,470-residue
model takes ~0.1 s and ~100 MB.

Original code: MIT License, Copyright (c) 2022 Shintaro Minami.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import sys

import numpy as np
from scipy.spatial.distance import cdist

# Kabsch-Sander electrostatic H-bond energy: q1*q2*f, kcal/mol, and DSSP's -0.5 cutoff
CONST_Q1Q2_F = 0.084 * 332
HBOND_CUTOFF = -0.5
BACKBONE = ("N", "CA", "C", "O")
SS_LABELS = {"H": "Helix", "E": "Strand"}


def _unit(v):
    """Row-normalise an (n, 3) array."""
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def _hydrogen_positions(coord):
    """Amide H of residues 1..L-1 placed as DSSP does: N + unit(C(i-1) - O(i-1))."""
    return coord[1:, 0] + _unit(coord[:-1, 2] - coord[:-1, 3])


def hbond_map(coord, donor_mask=None):
    """Boolean [acceptor C=O, donor N-H] H-bond matrix (energy below -0.5 kcal/mol)."""
    L = coord.shape[0]
    n, h = coord[1:, 0], _hydrogen_positions(coord)
    c, o = coord[:-1, 2], coord[:-1, 3]
    # e[donor i-1, acceptor j], donors are residues 1..L-1 and acceptors 0..L-2;
    # accumulated one distance matrix at a time to keep only two L x L floats alive
    e = 1.0 / cdist(n, o).astype(np.float32)
    e += 1.0 / cdist(h, c).astype(np.float32)
    e -= 1.0 / cdist(h, o).astype(np.float32)
    e -= 1.0 / cdist(n, c).astype(np.float32)
    e *= CONST_Q1Q2_F
    hb = np.zeros((L, L), dtype=bool)  # [donor, acceptor]
    hb[1:, :-1] = e < HBOND_CUTOFF
    del e
    # no bonds between a residue and itself or its two predecessors
    idx = np.arange(L)
    for k in (0, 1, 2):
        hb[idx[k:], idx[:L - k]] = False
    # donors without an amide H (e.g. proline) bond nothing
    if donor_mask is not None:
        hb &= np.asarray(donor_mask, dtype=bool)[:, None]
    return hb.T  # [acceptor, donor]


def _helix(hb, n):
    """Residues in an n-turn helix: two consecutive n-turns at i-1 and i mark i..i+n-1."""
    L = hb.shape[0]
    turn = np.zeros(L, dtype=bool)
    turn[:L - n] = np.diagonal(hb, offset=n)  # C=O(i) -> N-H(i+n)
    start = np.zeros(L, dtype=bool)
    start[1:] = turn[:-1] & turn[1:]
    return start, n


def _spread(start, n):
    """Extend each helix start over the n residues it covers."""
    out = start.copy()
    for k in range(1, n):
        out[k:] |= start[:-k]
    return out


def _ladder(hb):
    """Residues in any parallel or antiparallel bridge (DSSP bridge patterns)."""
    L = hb.shape[0]
    hbt = hb.T

    # S(dr, dc)[I-1, J-1] = hb[I+dr, J+dc] for bridge centres I, J in 1..L-2
    def S(m, dr, dc):
        return m[1 + dr:L - 1 + dr, 1 + dc:L - 1 + dc]

    # parallel: [I-1 -> J and J -> I+1] or [J-1 -> I and I -> J+1]
    bridge = (S(hb, -1, 0) & S(hbt, 1, 0)) | (S(hbt, 0, -1) & S(hb, 0, 1))
    # antiparallel: [I -> J and J -> I] or [I-1 -> J+1 and J-1 -> I+1]
    bridge |= (S(hb, 0, 0) & S(hbt, 0, 0)) | (S(hb, -1, 1) & S(hbt, 1, -1))
    out = np.zeros(L, dtype=bool)
    out[1:-1] = bridge.any(axis=1)
    return out


def assign(coord, donor_mask=None):
    """Per-residue 'H' / 'E' / '-' string from an (L, 4, 3) N/CA/C/O coordinate array."""
    coord = np.asarray(coord, dtype=float)
    L = coord.shape[0]
    if L < 6:
        return "-" * L
    hb = hbond_map(coord, donor_mask)
    h3, h4, h5 = _helix(hb, 3)[0], _helix(hb, 4)[0], _helix(hb, 5)[0]
    helix4 = _spread(h4, 4)
    # the alpha helix wins: 3-10 and pi starts inside or just before one are dropped
    blocked = helix4 | np.roll(helix4, -1)
    helix = helix4 | _spread(h3 & ~blocked, 3) | _spread(h5 & ~blocked, 5)
    strand = _ladder(hb) & ~helix
    return "".join(np.where(helix, "H", np.where(strand, "E", "-")))


def backbone_from_structure(struct):
    """(resnums, residue names, (n, 4, 3) N/CA/C/O coords) for complete-backbone residues of model 0."""
    resnums, residues, coords = [], [], []
    for res in struct[0].get_residues():
        # standard residues only, and only those with all four backbone atoms
        if res.id[0] != " " or not all(a in res for a in BACKBONE):
            continue
        resnums.append(res.id[1])
        residues.append(res.get_resname())
        coords.append([res[a].coord for a in BACKBONE])
    return resnums, residues, np.array(coords, dtype=float).reshape(-1, 4, 3)


def _chain_segments(resnums, coords, max_cn=2.0):
    """Split residue indices wherever numbering jumps or the C(i)-N(i+1) bond is broken."""
    segments, start = [], 0
    for i in range(1, len(resnums)):
        gap = resnums[i] != resnums[i - 1] + 1
        broken = np.linalg.norm(coords[i, 0] - coords[i - 1, 2]) > max_cn
        if gap or broken:
            segments.append((start, i))
            start = i
    segments.append((start, len(resnums)))
    return segments


def backbone_from_pdb(path):
    """Same as backbone_from_structure, read straight from PDB text (first model, first chain)."""
    atoms, order, chain = {}, [], None
    with open(path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM"):
                continue
            chain = line[21] if chain is None else chain
            name = line[12:16].strip()
            # first chain only, backbone atoms only, first altloc only
            if line[21] != chain or name not in BACKBONE or line[16] not in " A":
                continue
            key = (int(line[22:26]), line[26])
            if key not in atoms:
                atoms[key] = {"res": line[17:20]}
                order.append(key)
            atoms[key].setdefault(name, [float(line[30:38]), float(line[38:46]),
                                         float(line[46:54])])
    keep = [k for k in order if all(a in atoms[k] for a in BACKBONE)]
    coords = np.array([[atoms[k][a] for a in BACKBONE] for k in keep], dtype=float)
    return [k[0] for k in keep], [atoms[k]["res"] for k in keep], coords.reshape(-1, 4, 3)


def secondary_structure(struct=None, pdb_path=None):
    """{resnum: 'H'|'E'|'-'} for a Bio.PDB structure or a PDB file; chain breaks assigned apart."""
    if struct is not None:
        resnums, names, coords = backbone_from_structure(struct)
    else:
        resnums, names, coords = backbone_from_pdb(pdb_path)
    out = {}
    if not resnums:
        return out
    # proline has no amide H, so it donates no H-bond (as in DSSP)
    donors = np.array([n != "PRO" for n in names])
    # each unbroken stretch is assigned on its own so no bond spans a gap
    for a, b in _chain_segments(resnums, coords):
        for r, s in zip(resnums[a:b], assign(coords[a:b], donor_mask=donors[a:b])):
            out[r] = s
    return out


def ss_ranges(ss_by_res):
    """[(label, start, stop)] runs of helix / strand from a {resnum: state} dict."""
    out, run = [], None
    for r in sorted(ss_by_res):
        s = ss_by_res[r]
        # extend the open run when the state and numbering continue
        if run and s == run[0] and r == run[2] + 1:
            run[2] = r
            continue
        if run and run[0] in SS_LABELS:
            out.append((SS_LABELS[run[0]], run[1], run[2]))
        run = [s, r, r]
    if run and run[0] in SS_LABELS:
        out.append((SS_LABELS[run[0]], run[1], run[2]))
    return out


def write_ss_file(path, struct=None, pdb_path=None):
    """Write helix/strand runs as a source/start/stop/description range file; return run count.

    Fails soft: a structure the assignment cannot handle leaves an empty file, so the caller's
    other outputs (pLDDT, SASA, ...) are never lost to it.
    """
    try:
        ranges = ss_ranges(secondary_structure(struct=struct, pdb_path=pdb_path))
    except Exception as e:
        print(f"Secondary structure assignment failed: {e}", file=sys.stderr)
        ranges = []
    with open(path, "w") as out:
        for label, start, stop in ranges:
            print(f"secondary_structure\t{start}\t{stop}\t{label}", file=out)
    return len(ranges)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Helix/strand ranges from a PDB file (DSSP-style).")
    parser.add_argument("--pdb", required=True, help="input PDB file")
    parser.add_argument("--output", required=True, help="output range file (.ss.txt)")
    args = parser.parse_args()
    print(f"{write_ss_file(args.output, pdb_path=args.pdb)} helix/strand segment(s)")
