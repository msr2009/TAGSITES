"""paired_features.py

UniProt gives an intrachain disulfide bond (and an intrachain cross-link) as a location
whose start and end are the two bonded residues, not a stretch of sequence. Shown as a
start..end range, a bond between C45 and C112 painted all 68 residues as a blocking
UniProt feature. This splits each such feature into two single-residue rows that share
one description naming the pair (e.g. "Disulfide bond: C45-C112 (PubMed:...)"), so both
partners get the same color and legend entry and nothing in between is marked.

Interchain bonds already arrive as one residue and pass through unchanged, as do rows
that have been split before, so applying this twice is harmless.

Dependency-free so uniprot_features.py, build_proteome_db.py and utils/results.py share it.

Matt Rich, 2026
"""

# UniProt feature labels (tables/uniprot_features.txt) whose start/end are two bonded residues
PAIRED_LABELS = ("Disulfide bond", "Cross-link")


def _pair_name(start, stop, sequence=None):
    """'C45-C112' when the sequence covers both residues, else '45-112'."""
    if sequence and stop <= len(sequence):
        return f"{sequence[start - 1]}{start}-{sequence[stop - 1]}{stop}"
    return f"{start}-{stop}"


def pair_description(description, start, stop, sequence=None):
    """Rewrite 'Label: detail (evidence)' as 'Label: C45-C112[, detail] (evidence)'."""
    label, _, rest = description.partition(": ")
    # UniProt often leaves the detail blank, so it falls back to the type name itself
    if rest.startswith(label):
        rest = rest[len(label):].lstrip(" ,;")
    pair = _pair_name(start, stop, sequence)
    if not rest:
        return f"{label}: {pair}"
    sep = " " if rest.startswith("(") else ", "
    return f"{label}: {pair}{sep}{rest}"


def is_paired(source, start, stop, description):
    """True for an unsplit UniProt bond row spanning two residues."""
    return (source == "UniProt" and int(stop) > int(start)
            and str(description).startswith(tuple(f"{lab}: " for lab in PAIRED_LABELS)))


def split_paired_rows(rows, sequence=None):
    """Replace each (source, start, stop, description, ...) bond row with its two residues.

    Extra trailing fields are carried onto both rows unchanged; every other row passes through.
    """
    out = []
    for row in rows:
        source, start, stop, desc = row[0], int(row[1]), int(row[2]), row[3]
        # anything that is not a two-residue bond keeps its span
        if not is_paired(source, start, stop, desc):
            out.append(row)
            continue
        new_desc = pair_description(desc, start, stop, sequence)
        for pos in (start, stop):
            out.append((source, pos, pos, new_desc, *row[4:]))
    return out
