"""
test_design_reagents.py — offline integration test for scripts/design_tag_reagents.py

Runs design_reagents() on the snb-1 fixtures (109 residues, no network needed)
and performs structural + correctness checks on the output DataFrame.
"""

import re
import sys
from pathlib import Path
import pytest
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from design_tag_reagents import design_reagents
from crispr_util import iupac_to_regex, reverse_complement


EXPECTED_COLS = [
    "residue_index", "amino_acid", "insert_pos", "exon_index",
    "is_split_codon", "dist_to_5p_splice", "dist_to_3p_splice",
    "guide_strand", "spacer", "pam_seq", "pam_fwd_start", "cut_pos",
    "distance", "pam_in_arm", "recut_block_method", "mutation_desc",
    "left_arm", "right_arm",
]

PAM     = "NGG"
GUIDE_L = 20
ARM_L   = 200   # smaller arm to keep test fast
N_GUIDES = 3


@pytest.fixture(scope="module")
def snb1_df(DATA):
    df, _genotyping_df = design_reagents(
        genewise_out=DATA / "snb-1.genewise.out.txt",
        genomic_fasta=DATA / "snb-1_genomic.fa",
        n_guides=N_GUIDES,
        arm_length=ARM_L,
        pam=PAM,
        guide_length=GUIDE_L,
        cut_offset=3,
    )
    return df


@pytest.fixture(scope="module")
def snb1_dna(DATA):
    return str(list(SeqIO.parse(DATA / "snb-1_genomic.fa", "fasta"))[0].seq)


# ── structure ─────────────────────────────────────────────────────────────────

def test_columns_present(snb1_df):
    for col in EXPECTED_COLS:
        assert col in snb1_df.columns, f"missing column: {col}"


def test_not_empty(snb1_df):
    assert len(snb1_df) > 0


def test_residue_index_range(snb1_df):
    assert snb1_df["residue_index"].min() == 1
    assert snb1_df["residue_index"].max() == 109


def test_guides_per_residue_at_most_n(snb1_df):
    counts = snb1_df.groupby("residue_index").size()
    assert (counts <= N_GUIDES).all()


# ── guide properties ──────────────────────────────────────────────────────────

def test_spacer_length(snb1_df):
    assert (snb1_df["spacer"].str.len() == GUIDE_L).all()


def test_pam_matches_ngg(snb1_df):
    pattern = re.compile(iupac_to_regex(PAM))
    bad = snb1_df[~snb1_df["pam_seq"].str.match(pattern.pattern)]
    assert len(bad) == 0, f"non-NGG PAMs found:\n{bad[['guide_strand','pam_seq']].head()}"


def test_arm_length_at_most_arm_l(snb1_df):
    assert (snb1_df["left_arm"].str.len() <= ARM_L).all()
    assert (snb1_df["right_arm"].str.len() <= ARM_L).all()


# ── recut blocking ────────────────────────────────────────────────────────────

def test_no_none_recut_method(snb1_df):
    """Every row should have a recut-blocking strategy."""
    assert (snb1_df["recut_block_method"] != "none").all()


def test_valid_recut_methods(snb1_df):
    valid = {"syn_1", "mut_1", "insertion"}
    bad = set(snb1_df["recut_block_method"].unique()) - valid
    assert not bad, f"unexpected recut methods: {bad}"


# ── spacer-in-genomic verification ───────────────────────────────────────────

def test_plus_strand_spacer_in_genome(snb1_df, snb1_dna):
    plus = snb1_df[snb1_df["guide_strand"] == "+"].head(5)
    for _, row in plus.iterrows():
        assert row["spacer"] in snb1_dna, (
            f"+ strand spacer {row['spacer']} not found in genomic sequence"
        )


def test_minus_strand_spacer_rc_in_genome(snb1_df, snb1_dna):
    minus = snb1_df[snb1_df["guide_strand"] == "-"].head(5)
    for _, row in minus.iterrows():
        rc = reverse_complement(row["spacer"])
        assert rc in snb1_dna, (
            f"- strand spacer RC({row['spacer']}) not found in genomic sequence"
        )


# ── PAM disruption sanity ─────────────────────────────────────────────────────

def test_mutated_arms_differ_where_expected(snb1_df, snb1_dna):
    """
    For rows where pam_in_arm is 'left' or 'right' and recut_block_method is not
    'insertion', the arm should differ from the raw genomic sequence at least once.
    """
    subset = snb1_df[
        snb1_df["pam_in_arm"].isin(["left", "right"]) &
        (snb1_df["recut_block_method"] != "insertion") &
        (snb1_df["recut_block_method"] != "none")
    ].head(10)

    for _, row in subset.iterrows():
        insert_pos = int(row["insert_pos"])
        arm_l = len(row["left_arm"])
        arm_r = len(row["right_arm"])
        raw_left  = snb1_dna[insert_pos - arm_l: insert_pos]
        raw_right = snb1_dna[insert_pos: insert_pos + arm_r]
        if row["pam_in_arm"] == "left":
            assert row["left_arm"] != raw_left, (
                f"left_arm unchanged for residue {row['residue_index']}"
            )
        else:
            assert row["right_arm"] != raw_right, (
                f"right_arm unchanged for residue {row['residue_index']}"
            )


# ── Bad-input validation (unc119 longpro/shortdna fixtures) ──────────────────

DATA_DIR = Path(__file__).parent / "data"
LONGPRO_SHORTDNA = DATA_DIR / "unc119_longpro_shortdna"
LONGPRO_LONGDNA  = DATA_DIR / "unc119_longpro_longdna"
UNC119_PROTEIN_LEN = 244   # G5EGP9, the long isoform

_MISSING_FIXTURES = pytest.mark.skip(
    reason="fixture data never committed: needs a real UNC119 (G5EGP9) long-isoform "
           "Genewise run against short-DNA and long-DNA genomic regions"
)


@_MISSING_FIXTURES
class TestShortDnaIsoformMismatch:
    """
    longpro_shortdna: 244 aa long-isoform protein aligned to the 1335 nt
    short-isoform genomic (3 exons).  Genewise score is high (428 bits) but
    CDS coverage is only 82%, which is below the 90% threshold.  The task
    should raise ValueError rather than silently producing truncated reagents.
    """

    @pytest.fixture(scope="class", autouse=False)
    @classmethod
    def paths(cls):
        return {
            "gw": LONGPRO_SHORTDNA / "unc119_longpro_shortdna_genewise.genewise.out.txt",
            "fa": LONGPRO_SHORTDNA / "unc119_longpro_shortdna_genewise.genewise_genomic.fa",
        }

    def test_raises_on_low_coverage(self, paths):
        with pytest.raises(ValueError, match="Mismatch between DNA and protein"):
            design_reagents(
                genewise_out=paths["gw"],
                genomic_fasta=paths["fa"],
                protein_length=UNC119_PROTEIN_LEN,
            )

    def test_error_mentions_coverage(self, paths):
        with pytest.raises(ValueError, match="coverage"):
            design_reagents(
                genewise_out=paths["gw"],
                genomic_fasta=paths["fa"],
                protein_length=UNC119_PROTEIN_LEN,
            )

    def test_high_score_does_not_mask_bad_coverage(self, paths):
        """Score 428 bits is high; the error must come from coverage, not score."""
        import scripts.parse_genewise as pg
        score = pg.parse_genewise_gff_score(str(paths["gw"]))
        assert score > 50.0, "fixture score should be above score threshold"
        # Now confirm design_reagents still raises despite the good score
        with pytest.raises(ValueError):
            design_reagents(
                genewise_out=paths["gw"],
                genomic_fasta=paths["fa"],
                protein_length=UNC119_PROTEIN_LEN,
            )


@_MISSING_FIXTURES
class TestLongDnaCorrectCase:
    """longpro_longdna: correct protein + correct genomic — must pass all checks."""

    @pytest.fixture(scope="class")
    @classmethod
    def result(cls):
        gw = LONGPRO_LONGDNA / "unc119_longpro_longdna_genewise.genewise.out.txt"
        fa = LONGPRO_LONGDNA / "unc119_longpro_longdna_genewise.genewise_genomic.fa"
        df, _genotyping_df = design_reagents(
            genewise_out=gw,
            genomic_fasta=fa,
            protein_length=UNC119_PROTEIN_LEN,
        )
        return df

    def test_produces_reagents(self, result):
        assert len(result) > 0

    def test_all_arms_full_length(self, result):
        """With adequate flanking, every arm should be exactly arm_length=1000."""
        assert result["left_arm"].str.len().min() == 1000
        assert result["right_arm"].str.len().min() == 1000

    def test_residue_count_matches_protein(self, result):
        """One row per residue × guide; unique residues should equal protein length."""
        assert result["residue_index"].nunique() == UNC119_PROTEIN_LEN
