"""
pfam_scan.py

Shared local Pfam-A scanning engine, used by scripts/build_pfam_cache.py
(bulk pre-scan over a whole proteome) and scripts/domains_scan.py (per-protein
backend, on-demand fallback for cache misses). Runs pyhmmer.hmmsearch() —
Pfam-A HMMs as queries, protein sequences as the target database — which is
HMMER's efficient direction: model setup cost is paid once per profile and
amortized across every sequence scanned in the same call, unlike hmmscan
(sequence as query, HMMs as the searched database), which pays that cost
per sequence. See ~/.claude/plans/i-m-considering-a-large-humble-sun.md.

Thresholds use Pfam's own "gathering" bit-score cutoffs (--cut_ga in the
HMMER CLI) by default, matching what InterProScan's Pfam member analysis
uses, rather than an arbitrary e-value cutoff.

Coordinates: Pfam/InterProScan report envelope coordinates (the region
HMMER's forward algorithm assigns to the domain, not just the aligned core),
so this module reads Domain.env_from/env_to, not the narrower alignment
coordinates.
"""

import json
import sys
from pathlib import Path

import pyhmmer

sys.path.insert(0, str(Path(__file__).parent))
from reference_data import _load_config as _load_reference_data_config

_REPO_ROOT = Path(__file__).parent.parent
_CONFIG_PATH = _REPO_ROOT / "batch.config.json"

_DEFAULTS = {
    "bit_cutoffs": "gathering",
    "clan_filter": True,
    "cache_profiles": True,
    "cpus": 0,
}

_profile_cache = {}  # cfg-path-independent module cache: {hmm_path: [profiles]}


def _load_scan_config():
    """Merge batch.config.json's "domains_scan" block (if any) over _DEFAULTS."""
    cfg = dict(_DEFAULTS)
    if _CONFIG_PATH.exists():
        with open(_CONFIG_PATH) as f:
            user_cfg = json.load(f).get("domains_scan", {})
        cfg.update(user_cfg)
    return cfg


def pfam_hmm_path():
    """Resolve Pfam-A.hmm's path from reference_data's out_dir; raise a clear
    error naming the fetch command if it hasn't been downloaded yet.
    """
    ref_cfg = _load_reference_data_config()
    path = Path(ref_cfg["out_dir"]) / "Pfam-A.hmm"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — run `python scripts/reference_data.py --only pfam` first."
        )
    return path


def pfam_clans_path():
    """Resolve Pfam-A.clans.tsv's path (see pfam_hmm_path); may not exist if
    clan filtering was never enabled — callers should check existence.
    """
    ref_cfg = _load_reference_data_config()
    return Path(ref_cfg["out_dir"]) / "Pfam-A.clans.tsv"


def load_clans():
    """Parse Pfam-A.clans.tsv into {pfam_accession: clan_id}; entries with no
    clan (blank 2nd column) are omitted, since they can't overlap by clan.
    """
    path = pfam_clans_path()
    if not path.exists():
        return {}
    clans = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1]:
                clans[parts[0]] = parts[1]
    return clans


def load_profiles(scan_cfg=None):
    """Load and cache Pfam-A's ~21k HMM profiles from disk, keyed by the
    resolved hmm path so a second call in the same process is instant. This
    is the expensive one-time cost hmmsearch's direction is chosen to pay
    only once (see module docstring).
    """
    hmm_path = pfam_hmm_path()
    scan_cfg = scan_cfg or _load_scan_config()
    if not scan_cfg.get("cache_profiles", True):
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            return list(hmm_file)
    if hmm_path not in _profile_cache:
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            _profile_cache[hmm_path] = list(hmm_file)
    return _profile_cache[hmm_path]


def _bit_cutoffs_kwarg(scan_cfg):
    """Translate the "bit_cutoffs" config value into pyhmmer's hmmsearch
    kwarg; None means use plain e-value thresholding instead of Pfam's
    built-in per-model cutoffs.
    """
    value = scan_cfg.get("bit_cutoffs", "gathering")
    return {"bit_cutoffs": value} if value else {}


def _rows_from_tophits(query_name_unused, top_hits):
    """Flatten one query HMM's TopHits into (seq_name, pfam_acc, description,
    start, stop, score) rows — only included hits/domains (those passing the
    configured threshold), using envelope coordinates (see module docstring).
    pyhmmer >=0.12 exposes HMM/Hit name/accession/description as plain str
    (older versions returned bytes) — no .decode() needed here.
    """
    pfam_acc = top_hits.query.accession.split(".")[0] if top_hits.query.accession else None
    description = top_hits.query.description or top_hits.query.name or pfam_acc
    rows = []
    for hit in top_hits:
        if not hit.included:
            continue
        seq_name = hit.name
        for domain in hit.domains.included:
            rows.append((seq_name, pfam_acc, description, domain.env_from, domain.env_to,
                         domain.score))
    return rows


def filter_clan_overlaps(rows, clans):
    """Drop lower-scoring hits that overlap (share any residue) with a
    higher-scoring hit from the same Pfam clan, per sequence — matches real
    Pfam/InterProScan behavior of reporting only the best clan member over a
    given region instead of every clan paralog that happens to match.
    rows: list of (seq_name, pfam_acc, description, start, stop, score).
    """
    by_seq = {}
    for row in rows:
        by_seq.setdefault(row[0], []).append(row)

    kept = []
    for seq_rows in by_seq.values():
        seq_rows.sort(key=lambda r: -r[5])  # best score first
        accepted = []
        for row in seq_rows:
            _, pfam_acc, _, start, stop, _ = row
            clan = clans.get(pfam_acc)
            overlaps_better = False
            if clan is not None:
                for accepted_row in accepted:
                    if clans.get(accepted_row[1]) != clan:
                        continue
                    if start <= accepted_row[4] and accepted_row[3] <= stop:
                        overlaps_better = True
                        break
            if not overlaps_better:
                accepted.append(row)
        kept.extend(accepted)
    return kept


def scan_sequences(records, cfg=None):
    """Run hmmsearch (Pfam-A profiles as queries, `records` as the target
    database) and return {seq_name: [(pfam_acc, description, start, stop, score), ...]}
    for every sequence, including sequences with zero hits (empty list).
    records: iterable of (name, sequence_str) pairs.
    """
    scan_cfg = cfg or _load_scan_config()
    alphabet = pyhmmer.easel.Alphabet.amino()

    records = list(records)
    seq_block = pyhmmer.easel.DigitalSequenceBlock(alphabet, [
        pyhmmer.easel.TextSequence(name=name.encode(), sequence=str(seq)).digitize(alphabet)
        for name, seq in records
    ])

    profiles = load_profiles(scan_cfg)
    results = {name: [] for name, _ in records}

    for top_hits in pyhmmer.hmmsearch(
        profiles, seq_block, cpus=scan_cfg.get("cpus", 0), **_bit_cutoffs_kwarg(scan_cfg)
    ):
        for seq_name, pfam_acc, description, start, stop, score in _rows_from_tophits(None, top_hits):
            results[seq_name].append((pfam_acc, description, start, stop, score))

    if scan_cfg.get("clan_filter", True):
        clans = load_clans()
        if clans:
            for seq_name, hits in results.items():
                flat = [(seq_name, acc, desc, start, stop, score)
                        for acc, desc, start, stop, score in hits]
                filtered = filter_clan_overlaps(flat, clans)
                results[seq_name] = [(acc, desc, start, stop, score)
                                     for _, acc, desc, start, stop, score in filtered]

    return results
