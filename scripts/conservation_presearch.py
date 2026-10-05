"""
conservation_presearch.py

Batched DIAMOND search for the proteome run: searches thousands of proteins per
`diamond blastp` call (database setup dominates a per-protein call, ~55 s each
against Swiss-Prot + Rhabditida TrEMBL, versus ~0.2 s per protein batched) and
caches each sequence's raw hit list at
{reference_data.out_dir}/conservation_hits/{crc64}.json. conservation_local.main()
reads that cache instead of searching, so MAFFT + JSD scoring run per protein as
usual. Hits are identical to a single-query search (verified for snb-1 against
batches of 1, 50 and 500 queries).

Same databases, e-values, max-target-seqs and sensitivity as conservation_local
(batch.config.json's conservation_local.search_databases). Knobs, also in
conservation_local: threads (DIAMOND --threads) and chunk_size (queries per call).

Usage
-----
    python scripts/conservation_presearch.py                    # every local_store protein
    python scripts/conservation_presearch.py --fasta seqs.fa    # a FASTA instead
    python scripts/conservation_presearch.py --limit 500        # benchmark
    python scripts/conservation_presearch.py --threads 16 --chunk-size 2000
"""

import gzip
import json
import os
import platform
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

from Bio.SeqUtils.CheckSum import crc64

sys.path.insert(0, str(Path(__file__).parent))
from conservation_local import (
    _load_search_databases,
    _local_config,
    _run_diamond_blastp_by_query,
    hits_cache_dir,
)

DEFAULT_CHUNK_SIZE = 2000


def read_fasta_records(path):
    """Return [(name, sequence)] from a (optionally gzipped) FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    records, name, chunks = [], None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def local_store_records(limit=None):
    """Return [(accession, sequence)] from local_store's proteins table."""
    from local_store import open_index

    conn = open_index()
    try:
        query = "SELECT accession, sequence FROM proteins ORDER BY accession"
        if limit:
            query += f" LIMIT {int(limit)}"
        return conn.execute(query).fetchall()
    finally:
        conn.close()


def write_hits_file(path, hits):
    """Atomically write one sequence's hit list so a killed run leaves no partial file."""
    tmp_path = Path(str(path) + ".tmp")
    with open(tmp_path, "w") as f:
        json.dump(hits, f)
    os.replace(tmp_path, path)


def append_meta_run(cache_dir, run):
    """Append one invocation's record to conservation_hits/_meta.json."""
    meta_path = cache_dir / "_meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {"runs": []}
    meta["runs"].append(run)
    meta_path.write_text(json.dumps(meta, indent=2))


def run_presearch(records, threads=None, chunk_size=None):
    """Search every not-yet-cached unique sequence in records; return counts."""
    local_cfg = _local_config()
    threads = threads or local_cfg.get("threads")
    chunk_size = chunk_size or local_cfg.get("chunk_size") or DEFAULT_CHUNK_SIZE
    db_entries = _load_search_databases()

    cache_dir = hits_cache_dir()
    cache_dir.mkdir(parents=True, exist_ok=True)
    # stale .tmp files come from a killed run; the final names are never partial
    for stale in cache_dir.glob("*.tmp"):
        stale.unlink()

    # dedupe by CRC64 (the cache key); one directory listing, not a stat per sequence,
    # since out_dir may be a network mount
    already = {p.name for p in cache_dir.iterdir()}
    unique = {}
    for _name, seq in records:
        key = crc64(str(seq))
        if f"{key}.json" not in already:
            unique[key] = str(seq)
    # shortest-first keeps early chunks fast and the long tail together
    todo = sorted(unique.items(), key=lambda kv: len(kv[1]))
    print(
        f"[conservation_presearch] {len(records)} sequences, {len(unique)} not yet cached, "
        f"{len(db_entries)} databases, threads={threads or 'default'}, chunk={chunk_size}"
    )
    if not todo:
        print("[conservation_presearch] nothing to do")
        return {"searched": 0, "with_hits": 0}

    run = {
        "host": platform.node(),
        "started": datetime.now().isoformat(timespec="seconds"),
        "databases": [{k: v for k, v in e.items() if k != "path"} for e in db_entries],
        "threads": threads,
        "chunk_size": chunk_size,
    }
    start = time.time()
    n_done, n_with_hits = 0, 0
    try:
        for i in range(0, len(todo), chunk_size):
            chunk = todo[i : i + chunk_size]
            chunk_start = time.time()
            with tempfile.TemporaryDirectory() as tmpdir:
                query_fasta = Path(tmpdir) / "queries.fa"
                query_fasta.write_text("".join(f">{key}\n{seq}\n" for key, seq in chunk))
                by_query = _run_diamond_blastp_by_query(
                    query_fasta, db_entries, None, tmpdir, threads
                )
            # every query id is present in by_query, an empty list meaning "no hits"
            for key, _seq in chunk:
                hits = by_query[key]
                write_hits_file(cache_dir / f"{key}.json", hits)
                n_with_hits += bool(hits)
            n_done += len(chunk)
            elapsed = time.time() - chunk_start
            print(
                f"[conservation_presearch] {n_done}/{len(todo)}: chunk of {len(chunk)} in "
                f"{elapsed:.1f}s ({elapsed / len(chunk):.2f}s/seq)",
                flush=True,
            )
    finally:
        # recorded even on Ctrl-C so a partial run is still attributed
        run.update(
            {
                "finished": datetime.now().isoformat(timespec="seconds"),
                "n_searched": n_done,
                "n_with_hits": n_with_hits,
                "wall_seconds": round(time.time() - start, 1),
            }
        )
        append_meta_run(cache_dir, run)
    print(
        f"[conservation_presearch] done: {n_done} searched, {n_with_hits} with hits "
        f"-> {cache_dir} ({time.time() - start:.1f}s)"
    )
    return {"searched": n_done, "with_hits": n_with_hits}


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fasta", default=None, help="protein FASTA (.gz ok); default: local_store"
    )
    parser.add_argument("--limit", type=int, default=None, help="only the first N sequences")
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="DIAMOND --threads (default: conservation_local.threads, else all)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=None,
        help="queries per DIAMOND call (default: conservation_local.chunk_size)",
    )
    args = parser.parse_args()

    if args.fasta:
        recs = read_fasta_records(args.fasta)
        recs = recs[: args.limit] if args.limit else recs
    else:
        recs = local_store_records(args.limit)
    run_presearch(recs, threads=args.threads, chunk_size=args.chunk_size)
