"""utils/tag_filters.py — user-driven isoform / topology restriction of taggable positions

Turns the Results-tab "Advanced / Rescoring" choices (which isoforms are visible, which
ones the tag should land in, constitutive-only mode, allowed Phobius topology labels) into
the set of positions that remain available for tagging. The complement of that set is fed
to utils.scoring.position_mask_criterion() so the restriction behaves exactly like any
other config mask criterion.

Kept separate from utils.scoring: this module is pure set arithmetic over isoform records
and range annotations, with no notion of scoring criteria, so it can also drive the
Reagents-tab isoform-specificity annotation and the CLI without pulling in the scoring
engine. Isoform record shape is the one produced by scripts/derive_isoforms.py and surfaced
by utils.results._build_isoform_pane: "present"/"skipped" spans in query coordinates, plus
"accession"/"name". Phobius descriptions are assumed already normalized to short labels by
utils.results._translate_phobius_desc ("Cytoplasmic", "Extracellular", "Transmembrane",
"Signal peptide") — callers outside the Results-tab load path must translate first.
"""


def isoform_key(iso, index):
    """Stable per-run identifier for an isoform record: accession, else name, else index."""
    return iso.get("accession") or iso.get("name") or f"idx{index}"


def isoform_positions(iso, seq_len):
    """Set of 1-based query positions covered by an isoform's `present` spans."""
    positions = set()
    for start, stop in iso.get("present", []):
        positions.update(range(max(1, int(start)), min(int(stop), seq_len) + 1))
    return positions


def isoform_allowed_positions(isoforms, visible_keys, tagged_keys, constitutive, seq_len):
    """Positions taggable under the isoform rules; None means "no isoform restriction".

    constitutive: intersection over every visible isoform (sequence common to all of them).
    otherwise: the tag must land in exactly the checked ("tagged") isoforms — present in at
    least one tagged isoform AND absent from every visible untagged one. Checking every
    visible isoform is equivalent to no restriction (nothing to discriminate against).
    """
    visible = [(isoform_key(iso, i), iso) for i, iso in enumerate(isoforms)
               if isoform_key(iso, i) in visible_keys]
    if len(visible) <= 1:
        return None  # nothing to discriminate between

    pos_by_key = {k: isoform_positions(iso, seq_len) for k, iso in visible}

    if constitutive:
        return set.intersection(*pos_by_key.values())

    tagged = [pos_by_key[k] for k in pos_by_key if k in tagged_keys]
    untagged = [pos_by_key[k] for k in pos_by_key if k not in tagged_keys]
    if not tagged:
        return set()  # nothing checked -> nothing taggable
    if not untagged:
        return None  # everything checked -> unrestricted
    return set.union(*tagged) - set.union(*untagged)


def topology_allowed_positions(range_df, allowed_labels, seq_len):
    """Positions inside Phobius rows whose (already-translated) label is in allowed_labels.

    Returns None ("no restriction") when allowed_labels is empty/falsy or range_df has no
    Phobius rows at all.
    """
    if not allowed_labels or range_df is None or range_df.empty:
        return None
    phobius = range_df[range_df["source"] == "Phobius"]
    if phobius.empty:
        return None
    allowed = set(allowed_labels)
    positions = set()
    for _, row in phobius[phobius["description"].isin(allowed)].iterrows():
        start, stop = int(row["start"]), int(row["stop"])
        positions.update(range(max(1, start), min(stop, seq_len) + 1))
    return positions


def topology_tag_allowed_positions(range_df, all_labels, checked_labels, seq_len):
    """Positions taggable under a Tag-checkbox matrix over Phobius topology labels.

    Mirrors isoform_allowed_positions' matrix semantics: all labels checked (or no
    Phobius data at all) -> None (unrestricted); none checked -> empty set (nothing
    taggable); otherwise restricted to the checked labels' spans.
    """
    if not all_labels:
        return None
    if not checked_labels:
        return set()
    if set(checked_labels) >= set(all_labels):
        return None
    return topology_allowed_positions(range_df, list(checked_labels), seq_len) or set()


def allowed_tag_positions(isoforms, visible_keys, tagged_keys, constitutive,
                          range_df, allowed_labels, seq_len):
    """Intersection of the isoform and topology restrictions; None when neither applies."""
    parts = [p for p in (
        isoform_allowed_positions(isoforms, visible_keys, tagged_keys, constitutive, seq_len),
        topology_allowed_positions(range_df, allowed_labels, seq_len),
    ) if p is not None]
    if not parts:
        return None
    return set.intersection(*parts)


def position_isoform_labels(isoforms, seq_len):
    """Map each 1-based position to the list of isoform labels (name/accession) covering it.

    Used to annotate reagent sites with their isoform specificity — every isoform is
    considered here regardless of any Results-tab show/hide or tag selection, since reagent
    design is never filtered by those controls.
    """
    labels_by_pos = {p: [] for p in range(1, seq_len + 1)}
    for iso in isoforms:
        label = iso.get("name") or iso.get("accession") or ""
        for p in isoform_positions(iso, seq_len):
            if p in labels_by_pos:
                labels_by_pos[p].append(label)
    return labels_by_pos


def describe_position_isoforms(pos, labels_by_pos, all_labels):
    """Human-readable isoform specificity for a position: "constitutive", "only in X", or a list.

    all_labels is the full set of isoform labels in the run, used to detect the
    constitutive case (present in every one of them).
    """
    present_in = labels_by_pos.get(pos, [])
    if not present_in:
        return ""
    if len(all_labels) > 1 and set(present_in) == set(all_labels):
        return "constitutive"
    if len(present_in) == 1:
        return f"only in {present_in[0]}"
    return ", ".join(present_in)
