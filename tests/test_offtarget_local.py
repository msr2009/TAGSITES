"""
Tests for scripts/offtarget_local.py against a tiny synthetic BLAST database.

The genome is two random sequences: chr1 holds the query locus and chr2 holds a mutated
copy of it plus planted near-copies of one guide. That is enough to confirm the searches
still find what they must (a duplicated segment, 1-3 mismatch guide off-targets, a
primer's second amplicon) and still ignore what they must (the locus itself, a site with
a 3' mismatch). Skipped when BLAST+ is not installed.
"""

import random
import shutil
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import offtarget_local as ol  # noqa: E402
import offtarget_screen as ots  # noqa: E402
from crispr_util import reverse_complement  # noqa: E402

pytestmark = pytest.mark.skipif(
    not all(shutil.which(t) for t in ("blastn", "makeblastdb", "blastdbcmd")),
    reason="BLAST+ not installed")

TAXID = "9999"
REGION = (2000, 2600)   # 0-based half-open slice of chr1 used as the query locus


def _mutate(seq, positions):
    """Substitute each listed position with a different base."""
    out = list(seq)
    for i in positions:
        out[i] = next(b for b in "ACGT" if b != out[i])
    return "".join(out)


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """Build the toy genome, its BLAST database, a config pointing at it, and the screens' inputs."""
    rng = random.Random(1)
    chr1 = "".join(rng.choice("ACGT") for _ in range(8000))
    chr2 = [rng.choice("ACGT") for _ in range(16000)]
    region = chr1[REGION[0]:REGION[1]]

    # a guide with an NGG PAM inside the locus, away from the primers at offsets 50 and 450
    off = next(i for i in range(100, 400) if region[i + 21:i + 23] == "GG")
    guide = region[off:off + 23]
    near = _mutate(guide, [2, 8, 13])             # 3 mismatches, all 5' of the last 5 nt
    bad_3p = _mutate(guide, [17])                 # 3' end broken: must not be reported

    # chr2: a mutated copy of the locus (mutations kept clear of the primer windows),
    # two identical copies of the 3-mismatch guide site, and the 3'-mismatch decoy
    copy_start = 3000
    chr2[copy_start:copy_start + len(region)] = _mutate(
        region, [i for i in range(120, 400, 20) if not off - 1 <= i <= off + 24])
    for pos in (8000, 9000):
        chr2[pos:pos + 23] = near
    chr2[12000:12023] = bad_3p
    chr2 = "".join(chr2)

    tmp = tmp_path_factory.mktemp("blastdb")
    fasta = tmp / "toy.fa"
    fasta.write_text(">chr1\n{}\n>chr2\n{}\n".format(chr1, chr2))
    cfg = ots.load_config()
    cfg["local"] = {"blastdb_dir": str(tmp), "threads": 1, "spacer_evalue": 100000,
                    "spacer_min_len": 7, "primer_tm_min": 0,
                    "genomes": {TAXID: {"db": "toy", "fasta": str(fasta)}}}
    ol.build_blastdb(TAXID, cfg)

    return {"cfg": cfg, "region": region, "guide": guide, "near": near, "bad_3p": bad_3p,
            "fwd": region[50:70], "rev": reverse_complement(region[450:470]),
            "copy_start": copy_start, "guide_off": off}


@pytest.fixture(scope="module")
def region_result(world):
    """Screen A on the locus, as the reagents pipeline would run it."""
    return ol.run_region_screen(world["region"], "", TAXID, cfg=world["cfg"])


def test_available_needs_a_built_database(world, tmp_path):
    """available() is true for the built database and false for an unbuilt one."""
    assert ol.available(TAXID, world["cfg"])
    empty = dict(world["cfg"], local=dict(world["cfg"]["local"], blastdb_dir=str(tmp_path)))
    assert not ol.available(TAXID, empty)
    assert not ol.available("1234", world["cfg"])


def test_region_screen_finds_the_duplicate_and_the_locus_itself(world, region_result):
    """The mutated chr2 copy is a duplicate; the chr1 query locus is self, keyed by chromosome."""
    dups = region_result["duplicates"]
    assert any(h["chrom"] == "chr2" and h["identity"] >= 90 for h in dups)
    assert not any(h["chrom"] == "chr1" for h in dups)
    lo, hi = region_result["self_spans"]["chr1"][0]
    assert lo <= REGION[0] + 1 and hi >= REGION[1] - 5


def test_spacer_screen_recovers_planted_mismatch_sites_at_every_copy(world, region_result):
    """Both identical 3-mismatch copies are reported with their PAM, not collapsed into one."""
    res = ol.run_spacer_screen([world["guide"]], "", TAXID, cfg=world["cfg"],
                               self_spans=region_result["self_spans"])
    sites = res["spacer_hits"][0]
    by_start = {s["subject_span"][0]: s for s in sites if s["acc"] == "chr2"}
    # the guide sits inside the duplicated segment, so its perfect copy is reported too
    in_copy = world["copy_start"] + world["guide_off"] + 1
    assert sorted(by_start) == sorted([in_copy, 8001, 9001])
    assert by_start[in_copy]["perfect"]
    assert all(by_start[p]["mismatches"] == 3 and by_start[p]["pam_ok"] for p in (8001, 9001))


def test_spacer_screen_ignores_a_broken_3prime_end_and_excludes_the_locus(world, region_result):
    """A 3' mismatch is not an off-target, and the guide's own locus is excluded by self span."""
    res = ol.run_spacer_screen([world["guide"]], "", TAXID, cfg=world["cfg"],
                               self_spans=region_result["self_spans"])
    sites = res["spacer_hits"][0]
    assert not any(s["acc"] == "chr1" for s in sites)
    assert 12001 not in [s["subject_span"][0] for s in sites]
    # without the exclusion the on-target site comes back perfect, on the right strand
    raw = ol.run_spacer_screen([world["guide"]], "", TAXID, cfg=world["cfg"])
    own = [s for s in raw["spacer_hits"][0] if s["acc"] == "chr1"]
    assert len(own) == 1 and own[0]["perfect"] and own[0]["pam_ok"]
    assert own[0]["h_strand"] == "+"


def test_spacer_screen_reads_the_pam_of_a_minus_strand_site(world):
    """A guide whose protospacer lies on the minus strand is found there with its PAM."""
    region = world["region"]
    i = next(i for i in range(100, 400) if region[i:i + 2] == "CC")
    guide = reverse_complement(region[i:i + 23])   # ends in NGG on the minus strand
    res = ol.run_spacer_screen([guide], "", TAXID, cfg=world["cfg"])
    own = [s for s in res["spacer_hits"][0] if s["acc"] == "chr1"]
    assert len(own) == 1
    assert own[0]["h_strand"] == "-" and own[0]["perfect"] and own[0]["pam_ok"]
    assert own[0]["pam"] == guide[20:]


def test_primer_screen_reports_the_second_amplicon_only(world, region_result):
    """The intended chr1 product is excluded; the chr2 copy's product is reported."""
    pairs = [{"id": "p", "fwd_seq": world["fwd"], "rev_seq": world["rev"]}]
    res = ol.run_primer_screen(pairs, TAXID, world["cfg"],
                               self_spans=region_result["self_spans"])["p"]
    amps = res["amplicons"]
    assert len(amps) == 1
    assert amps[0]["acc"].startswith("chr2:") and amps[0]["product_size"] == 420
    # with no self spans the on-target product is reported too
    both = ol.run_primer_screen(pairs, TAXID, world["cfg"], self_spans={})["p"]["amplicons"]
    assert sorted(a["acc"].split(":")[0] for a in both) == ["chr1", "chr2"]


def test_primer_screen_drops_sites_with_a_3prime_mismatch(world, region_result):
    """A primer whose 3' base is wrong has no binding site, so it forms no products."""
    broken = world["fwd"][:-1] + next(b for b in "ACGT" if b != world["fwd"][-1])
    res = ol.run_primer_screen([{"id": "p", "fwd_seq": broken, "rev_seq": world["rev"]}],
                               TAXID, world["cfg"], self_spans={})["p"]
    assert res["n_sites"]["fwd"] == 0 and res["amplicons"] == []


def test_shipped_config_keeps_the_spacer_evalue_loose():
    """A toy genome cannot show it, but at genome scale blastn drops 3-mismatch sites below 1e5."""
    assert ots.load_config()["local"]["spacer_evalue"] >= 1e5
