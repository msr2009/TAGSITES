"""
offtarget_local.py

Local BLAST+ backend for the off-target / primer-specificity screens. It runs against a
BLAST database built from the genome FASTA already in data/reference/, so nothing here
touches the network.

Same contract as offtarget_remote.py / offtarget_blat.py: run_region_screen() is Screen A
(the whole region as one query), run_spacer_screen() is Screen B (all spacer+PAMs in one
query) with the PAM settled from the local database rather than ENA fetches, and
EBI-shaped HSP dicts feed the unchanged offtarget_screen.py core. run_primer_screen() is
new: a Primer-BLAST-style check of genotyping primers.

Why this is not in the Shiny deployment: it needs blastn/blastdbcmd and a ~100 MB
database, neither of which exist on shinyapps.io. available() reports whether both are
present, and callers fall back to BLAT/EBI when it is False.

Coordinate frame. Screen A HSPs get acc = "chrom:lo-hi", one per local alignment, the same
choice offtarget_blat.py makes: classify_hsps calls an accession "self" when its hits cover
the query, so keying by chromosome would label every hit on the gene's own chromosome as
self. Screen B HSPs use acc = chromosome, so the spans run_region_screen() returns are
keyed by chromosome and own_accessions is empty (a whole chromosome is not "our locus").

Primer screening. primer3 cannot search a genome, so BLAST only seeds candidate sites
(blastn-short, word size 7). Each is then extended to the primer's full length from the
database, counted for mismatches, checked for an exact 3' end, and given a primer3
heterodimer Tm. Sites are paired Primer-BLAST style: any two sites on opposite strands
facing each other within scoring.max_product. Gapped primer binding is not modelled.

Known limit: blastn-short at word size 7 can miss a 20 nt off-target whose 3+ mismatches
are spread evenly, since no 7-mer seed survives. Cas-OFFinder would be exhaustive.

Matt Rich, 2026
"""

import argparse
import collections
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import offtarget_screen as ots
from crispr_util import reverse_complement
from progress import report as _report, resolve_reporter

_REPO_ROOT = Path(__file__).parent.parent
# sequences go last because they are the long columns
_OUTFMT = "6 qseqid sseqid qstart qend sstart send pident length bitscore evalue sstrand qseq sseq"
_chrom_length_cache = {}


# ── Database ──────────────────────────────────────────────────────────────────

def _local_cfg(cfg):
    """The offtarget.config.json "local" block ({} when absent)."""
    return (cfg or {}).get("local") or {}


def _threads(cfg):
    """blastn -num_threads: TAGSITES_BLAST_THREADS if set (a worker pool sets it to 1), else config."""
    return int(os.environ.get("TAGSITES_BLAST_THREADS") or _local_cfg(cfg).get("threads", 1))


def _abs(path):
    """Resolve a config path against the repo root."""
    p = Path(path)
    return p if p.is_absolute() else _REPO_ROOT / p


def genome_entry(taxid, cfg):
    """The configured {db, fasta} for a taxid, or None."""
    return (_local_cfg(cfg).get("genomes") or {}).get(str(taxid).strip())


def db_path(taxid, cfg):
    """Path prefix of the taxid's BLAST database, or None when no genome is configured."""
    entry = genome_entry(taxid, cfg)
    if not entry:
        return None
    return _abs(_local_cfg(cfg).get("blastdb_dir", "data/reference/blastdb")) / entry["db"]


def available(taxid, cfg):
    """True when blastn and blastdbcmd are installed and the taxid's database is built."""
    path = db_path(taxid, cfg)
    if path is None or not all(shutil.which(t) for t in ("blastn", "blastdbcmd")):
        return False
    return Path(str(path) + ".nin").exists()


def build_blastdb(taxid, cfg=None, force=False):
    """Build the taxid's BLAST database from its genome FASTA; skip if it already exists."""
    cfg = cfg or ots.load_config()
    entry, out = genome_entry(taxid, cfg), db_path(taxid, cfg)
    if out is None:
        raise RuntimeError("no local genome configured for taxid {}".format(taxid))
    if Path(str(out) + ".nin").exists() and not force:
        print("[skip] {} already built".format(out))
        return out
    if not shutil.which("makeblastdb"):
        raise RuntimeError("makeblastdb not found; install BLAST+ (mamba install -c bioconda blast)")
    fasta = _abs(entry["fasta"])
    if not fasta.exists():
        raise FileNotFoundError("{} not found; run reference_data.py --only genome".format(fasta))
    out.parent.mkdir(parents=True, exist_ok=True)
    print("[offtarget_local] building {} from {}".format(out, fasta))
    opener = gzip.open if str(fasta).endswith(".gz") else open
    # streamed through stdin so a gzipped FASTA never has to be unpacked on disk
    with tempfile.TemporaryFile("w+") as log:
        proc = subprocess.Popen(["makeblastdb", "-in", "-", "-dbtype", "nucl", "-parse_seqids",
                                 "-blastdb_version", "4",   # v5 needs LMDB, which network/external volumes reject
                                 "-title", entry["db"], "-out", str(out)],
                                stdin=subprocess.PIPE, stdout=log, stderr=subprocess.STDOUT)
        try:
            with opener(fasta, "rb") as f:
                shutil.copyfileobj(f, proc.stdin)
            proc.stdin.close()
        except BrokenPipeError:
            pass   # makeblastdb quit early; its own log below says why
        code = proc.wait()
        log.seek(0)
        text = log.read()
    if code != 0:
        raise RuntimeError("makeblastdb failed for {}: {}".format(fasta, text.strip()[-500:]))
    _chrom_length_cache.pop(str(out), None)
    return out


def chrom_lengths(db):
    """{sequence id: length} for a database, cached per process."""
    key = str(db)
    if key not in _chrom_length_cache:
        res = subprocess.run(["blastdbcmd", "-db", key, "-entry", "all", "-outfmt", "%a %l"],
                             capture_output=True, text=True, check=True)
        _chrom_length_cache[key] = {p[0]: int(p[1]) for p in
                                    (ln.split() for ln in res.stdout.splitlines()) if len(p) == 2}
    return _chrom_length_cache[key]


def fetch_ranges(db, ranges):
    """Plus-strand bases for [(chrom, start, end)] (1-based inclusive); '' where out of range."""
    lengths = chrom_lengths(db)
    valid = [i for i, (c, s, e) in enumerate(ranges)
             if c in lengths and 1 <= s <= e <= lengths[c]]
    out = [""] * len(ranges)
    if not valid:
        return out
    with tempfile.NamedTemporaryFile("w", suffix=".txt") as batch:
        for i in valid:
            c, s, e = ranges[i]
            batch.write("{} {}-{} plus\n".format(c, s, e))
        batch.flush()
        res = subprocess.run(["blastdbcmd", "-db", str(db), "-entry_batch", batch.name,
                              "-outfmt", "%s"], capture_output=True, text=True, check=True)
    for i, seq in zip(valid, res.stdout.split()):
        out[i] = seq.upper()
    return out


# ── blastn ────────────────────────────────────────────────────────────────────

def run_blastn(queries, db, task, word_size, evalue, max_targets, threads=1, min_len=0):
    """Run blastn on a list of query strings; row dicts carry "q", the query's list index."""
    rows = []
    with tempfile.NamedTemporaryFile("w", suffix=".fa") as qf, \
            tempfile.TemporaryFile("w+") as errf:
        for i, seq in enumerate(queries):
            qf.write(">q{}\n{}\n".format(i, seq))
        qf.flush()
        cmd = ["blastn", "-query", qf.name, "-db", str(db), "-task", task,
               "-word_size", str(word_size), "-evalue", str(evalue),
               "-max_target_seqs", str(max_targets), "-dust", "no", "-soft_masking", "false",
               "-num_threads", str(threads), "-outfmt", _OUTFMT]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errf, text=True)
        # streamed, not collected: a word-size-7 genome search emits millions of tiny hits
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13 or int(f[7]) < min_len:
                continue
            rows.append({"q": int(f[0][1:]), "chrom": f[1], "q_from": int(f[2]),
                         "q_to": int(f[3]), "s_from": int(f[4]), "s_to": int(f[5]),
                         "identity": float(f[6]), "align_len": int(f[7]),
                         "bits": float(f[8]), "evalue": float(f[9]),
                         "strand": "-" if f[10].startswith("minus") else "+",
                         "qseq": f[11], "sseq": f[12]})
        if proc.wait() != 0:
            errf.seek(0)
            raise RuntimeError("blastn failed: {}".format(errf.read().strip()[-400:]))
    return rows


def to_hsp(row, acc, db_name):
    """An EBI-shaped HSP dict (see offtarget_screen.parse_blast_json) from a blastn row."""
    lo, hi = sorted((row["s_from"], row["s_to"]))
    return {
        "acc":       acc,
        "chrom":     row["chrom"],
        "desc":      "{} {}:{}-{} ({})".format(db_name, row["chrom"], lo, hi, row["strand"]),
        "os":        db_name,
        "q_from":    min(row["q_from"], row["q_to"]),
        "q_to":      max(row["q_from"], row["q_to"]),
        "h_from":    row["s_from"],
        "h_to":      row["s_to"],
        "h_strand":  row["strand"],
        "qseq":      row["qseq"],
        "hseq":      row["sseq"],
        "identity":  row["identity"],
        "bits":      row["bits"],
        "evalue":    row["evalue"],
        "align_len": row["align_len"],
    }


def _footprint(row, primer_len):
    """Full-length subject span (chrom-local start, end, strand) implied by a seeding HSP."""
    q_lo, q_hi = sorted((row["q_from"], row["q_to"]))
    lo_pad, hi_pad = q_lo - 1, primer_len - q_hi   # primer bases beyond each end of the HSP
    s_lo, s_hi = sorted((row["s_from"], row["s_to"]))
    # the primer's 5' end sits lower on the subject for a plus hit and higher for a minus hit
    if row["strand"] == "+":
        return s_lo - lo_pad, s_hi + hi_pad, "+"
    return s_lo - hi_pad, s_hi + lo_pad, "-"


def _require(taxid, cfg):
    """Database path for a screen, raising if the local backend cannot run."""
    if not available(taxid, cfg):
        raise RuntimeError("local BLAST database for taxid {} is unavailable (blastn missing "
                           "or `reference_data.py --only blastdb` not run)".format(taxid))
    return db_path(taxid, cfg)


# ── Screen A ──────────────────────────────────────────────────────────────────

def run_region_screen(region_seq, email, taxid, exons=None, cfg=None, report=None,
                      job_id_cb=None, resume_job_ids=None):
    """Screen A against the local genome; same contract as offtarget_remote's."""
    cfg = cfg or ots.load_config()
    reporter = resolve_reporter(report)
    db = _require(taxid, cfg)
    name, blast = db.name, cfg["blast"]
    _report(reporter, "Running local region off-target blastn ({} bp vs {})".format(
        len(region_seq), name), stage="offtarget_region")

    t0 = time.perf_counter()
    rows = run_blastn([region_seq], db, "blastn", blast["wordsize_region"],
                      blast["evalue_region"], blast["alignments"],
                      threads=_threads(cfg))
    hsps = []
    for r in rows:
        lo, hi = sorted((r["s_from"], r["s_to"]))
        hsps.append(to_hsp(r, "{}:{}-{}".format(r["chrom"], lo, hi), name))
    elapsed = time.perf_counter() - t0

    ots.classify_hsps(hsps, len(region_seq), exons or [], cfg)
    hsps, n_dropped = ots.apply_post_filter(hsps, cfg)
    dups = ots.duplicates(hsps)
    identical = ots.identical_segments(hsps, cfg)
    n_self = sum(1 for h in hsps if h["klass"] == "self")
    n_tx = sum(1 for h in hsps if h["klass"] == "transcript")
    _report(reporter, "Region screen ({} blastn hits in {:.1f}s): {} duplicated segment(s) "
                      "({} identical), {} self, {} transcript, {} below thresholds".format(
                          len(rows), elapsed, len(dups), len(identical), n_self, n_tx,
                          n_dropped), stage="offtarget_region")
    if identical:
        _report(reporter, "{} identical duplicated segment(s) >= {} bp — every oligo in "
                          "those windows will act on both copies".format(
                              len(identical), cfg["identical"]["identical_min_len"]),
                stage="offtarget_region", level="warning")

    # keyed by chromosome so they line up with Screen B's acc = chromosome
    selves = [dict(h, acc=h["chrom"]) for h in hsps if h["klass"] == "self"]
    return {
        "database":     name,
        "is_fallback":  False,
        "duplicates":   dups,
        "identical":    identical,
        "excluded":     [h for h in hsps if h["klass"] in ("self", "transcript")],
        "self_spans":   ots.self_spans(selves),
        "own_accessions": [],
        "n_below_threshold": n_dropped,
    }


# ── Screen B ──────────────────────────────────────────────────────────────────

def spacer_hsps(rows, spacer_list, block_len, sep_len, db):
    """Full-length ungapped HSPs for every seed row, one per distinct site.

    blastn at +1/-3 will not extend an alignment over a mismatch near the spacer's 5' end,
    and score_window ignores any HSP that does not cover the whole spacer. Each row is
    therefore only a seed: the spacer+PAM footprint it implies is read from the database
    and compared directly, which also puts the PAM bases in the alignment. Sites with a
    gap (a bulge) are scored as the ungapped footprint, so they usually fall out.
    """
    n = len(spacer_list)
    cand, seen = [], set()
    for r in rows:
        a = ots.block_of(r["q_from"], block_len, sep_len, n)
        b = ots.block_of(r["q_to"], block_len, sep_len, n)
        if a is None or b is None or a[0] != b[0]:
            continue   # a seed has to sit inside one spacer block to map back to a spacer
        idx = a[0]
        seed = dict(r, q_from=a[1] + 1, q_to=b[1] + 1)
        fp = _footprint(seed, len(spacer_list[idx]))
        key = (idx, r["chrom"]) + fp
        if key not in seen:   # several seeds can imply the same site
            seen.add(key)
            # 1-based position in the whole query where this spacer's block starts
            cand.append((idx, r["chrom"], fp, r["q_from"] - a[1]))
    bases = fetch_ranges(db, [(c, fp[0], fp[1]) for _, c, fp, _ in cand])
    hsps = []
    for (idx, chrom, (s0, s1, strand), q0), seq in zip(cand, bases):
        full = spacer_list[idx].upper()
        if len(seq) != len(full):
            continue   # footprint ran off the end of the sequence
        site = reverse_complement(seq) if strand == "-" else seq
        matches = sum(1 for x, y in zip(full, site) if x == y)
        hsps.append({
            "acc": chrom, "chrom": chrom,
            "desc": "{} {}:{}-{} ({})".format(Path(str(db)).name, chrom, s0, s1, strand),
            "os": Path(str(db)).name,
            "q_from": q0, "q_to": q0 + len(full) - 1,
            "h_from": s1 if strand == "-" else s0, "h_to": s0 if strand == "-" else s1,
            "h_strand": strand, "qseq": full, "hseq": site,
            "identity": 100.0 * matches / len(full), "bits": 0.0, "evalue": 0.0,
            "align_len": len(full),
        })
    return hsps


def run_spacer_screen(spacers, email, taxid, pam="NGG", cfg=None, report=None,
                      job_id_cb=None, resume_job_ids=None, self_spans=None,
                      own_accessions=None):
    """Screen B against the local genome: all spacers in one query, PAM read from the database."""
    cfg = cfg or ots.load_config()
    reporter = resolve_reporter(report)
    db = _require(taxid, cfg)
    sq, blast, local = cfg["spacer_query"], cfg["blast"], _local_cfg(cfg)
    spacer_list = list(spacers or [])[: int(sq["max_spacers"])]
    if not spacer_list:
        return {"spacer_hits": {}, "spacers": [], "block_len": 0}
    if len(spacers or []) > len(spacer_list):
        _report(reporter, "Screening the first {} of {} spacers (max_spacers)".format(
            len(spacer_list), len(spacers)), stage="offtarget_spacer", level="warning")

    query, block_len = ots.build_spacer_query(spacer_list, sq["separator_len"])
    _report(reporter, "Running local spacer off-target blastn ({} spacers, {} bp, wordsize {})"
            .format(len(spacer_list), len(query), blast["wordsize_spacer"]),
            stage="offtarget_spacer")
    t0 = time.perf_counter()
    rows = run_blastn([query], db, "blastn-short", blast["wordsize_spacer"],
                      local.get("spacer_evalue", blast["evalue_spacer"]),
                      blast.get("alignments_spacer", blast["alignments"]),
                      threads=_threads(cfg), min_len=int(local.get("spacer_min_len", 7)))
    raw = spacer_hsps(rows, spacer_list, block_len, sq["separator_len"], db)
    # screen_spacer_hits collapses sites with identical sequence and PAM, which is right for
    # ENA (one locus appears once per assembly) but wrong for a single-assembly database,
    # where identical sequences at different positions are separate off-targets. The Nth
    # copy of each (spacer, matched sequence) goes to round N, so copies never share a call.
    rounds, seen_seq = {}, collections.Counter()
    for h in raw:
        key = (ots.block_of(h["q_from"], block_len, sq["separator_len"], len(spacer_list))[0],
               h["hseq"])
        rounds.setdefault(seen_seq[key], []).append(h)
        seen_seq[key] += 1
    hits, n_self = {}, 0
    for k in sorted(rounds):
        part, n = ots.screen_spacer_hits(
            rounds[k], spacer_list, block_len, sq["separator_len"], pam, cfg,
            exclude_spans=self_spans, exclude_accessions=set(own_accessions or []))
        n_self += n
        for idx, sites in part.items():
            hits.setdefault(idx, []).extend(sites)
    if not self_spans:
        _report(reporter, "No self-locus coordinates available, so each guide's own "
                          "on-target site is counted as a hit — counts are inflated",
                stage="offtarget_spacer", level="warning")
    _report(reporter, "Spacer screen ({} seeds -> {} sites in {:.1f}s): {} of {} spacers have a "
                      "PAM-bearing off-target ({} on-target matches excluded)".format(
                          len(rows), len(raw), time.perf_counter() - t0, len(hits),
                          len(spacer_list), n_self), stage="offtarget_spacer")
    return {"spacer_hits": hits, "spacers": spacer_list, "block_len": block_len}


# ── Genotyping primers ────────────────────────────────────────────────────────

def _score_site(primer, site, three_prime_len):
    """(mismatches, three_prime_ok) for a primer against the same-length site read 5'->3'."""
    mm = sum(1 for a, b in zip(primer, site) if a != b)
    tail = max(1, three_prime_len)
    return mm, primer[-tail:] == site[-tail:]


def primer_sites(primers, db, cfg, report=None):
    """Priming-competent genomic sites per primer sequence: {primer: [site dicts]}.

    A site is kept when its 3' end is exact, it has at most scoring.max_mismatch
    mismatches, and its primer3 heterodimer Tm reaches primer_tm_min.
    """
    import primer3

    local, scoring = _local_cfg(cfg), cfg["scoring"]
    unique = sorted({p.upper() for p in primers if p})
    if not unique:
        return {}
    # every primer is a query in ONE blastn call; each row's "q" maps back to its primer
    rows = run_blastn(unique, db, "blastn-short", local.get("primer_wordsize", 7),
                      local.get("primer_evalue", 1000), cfg["blast"]["alignments"],
                      threads=_threads(cfg),
                      min_len=int(local.get("primer_min_len", 10)))
    candidates = []   # (primer, chrom, start, end, strand)
    seen = set()
    for r in rows:
        primer = unique[r["q"]]
        s0, s1, strand = _footprint(r, len(primer))
        key = (primer, r["chrom"], s0, s1, strand)
        if key not in seen:   # several HSPs can seed the same site
            seen.add(key)
            candidates.append((primer, r["chrom"], s0, s1, strand))

    # one blastdbcmd call reads every footprint; revcomp puts minus sites in primer sense
    bases = fetch_ranges(db, [(c, s0, s1) for _, c, s0, s1, _ in candidates])
    out = {p: [] for p in unique}
    tm_min = float(local.get("primer_tm_min", 0))
    for (primer, chrom, s0, s1, strand), seq in zip(candidates, bases):
        if len(seq) != len(primer):
            continue   # footprint ran off the end of the sequence
        site_seq = reverse_complement(seq) if strand == "-" else seq
        mm, tail_ok = _score_site(primer, site_seq, scoring["three_prime_len"])
        if not tail_ok or mm > scoring["max_mismatch"]:
            continue
        # the primer anneals to the template strand, which is the site's complement
        tm = primer3.bindings.calc_heterodimer(primer, reverse_complement(site_seq)).tm
        if tm < tm_min:
            continue
        out[primer].append({"chrom": chrom, "start": s0, "end": s1, "strand": strand,
                            "mismatches": mm, "tm": round(tm, 1), "seq": site_seq})
    return out


def _in_self(site, self_spans):
    """True when a site lies wholly inside the query's own locus."""
    return any(site["start"] >= lo and site["end"] <= hi
               for lo, hi in (self_spans or {}).get(site["chrom"], []))


def pair_amplicons(fwd_sites, rev_sites, self_spans, max_product):
    """Predicted products between any two facing sites of one primer pair, intended one excluded.

    fwd_sites/rev_sites are the sites of the forward and reverse primer. The intended
    amplicon is the perfect forward/perfect reverse pair inside the query's own locus;
    everything else, including fwd/fwd and rev/rev products, is reported.
    """
    tagged = [("fwd", s) for s in fwd_sites] + [("rev", s) for s in rev_sites]
    plus = [t for t in tagged if t[1]["strand"] == "+"]
    minus = [t for t in tagged if t[1]["strand"] == "-"]
    amps, seen = [], set()
    for which_a, a in plus:
        for which_b, b in minus:
            # facing pair: plus site upstream, minus site downstream, on one chromosome
            if a["chrom"] != b["chrom"] or b["end"] < a["end"]:
                continue
            size = b["end"] - a["start"] + 1
            if size > max_product:
                continue
            intended = (which_a == "fwd" and which_b == "rev" and a["mismatches"] == 0
                        and b["mismatches"] == 0 and _in_self(a, self_spans)
                        and _in_self(b, self_spans))
            key = (a["chrom"], a["start"], b["end"])
            if intended or key in seen:
                continue
            seen.add(key)
            amps.append({
                "acc":          "{}:{}-{}".format(a["chrom"], a["start"], b["end"]),
                "desc":         "{} primers {}/{}".format(a["chrom"], which_a, which_b),
                "product_size": size,
                "fwd":          (a["start"], a["end"]),
                "rev":          (b["start"], b["end"]),
                "perfect":      a["mismatches"] == 0 and b["mismatches"] == 0,
                "fwd_mismatches": a["mismatches"],
                "rev_mismatches": b["mismatches"],
            })
    return sorted(amps, key=lambda x: (not x["perfect"], x["product_size"]))


def run_primer_screen(primers, taxid, cfg=None, self_spans=None, report=None):
    """Off-target amplicons per primer pair: {pair id: {"amplicons": [...], "n_sites": {...}}}.

    primers is [{"id", "fwd_seq", "rev_seq"}]. self_spans comes from run_region_screen()
    and is what marks the intended amplicon; without it every on-target product would be
    reported, so callers should only run this alongside a local region screen.
    """
    cfg = cfg or ots.load_config()
    reporter = resolve_reporter(report)
    db = _require(taxid, cfg)
    t0 = time.perf_counter()
    seqs = [p[k] for p in primers for k in ("fwd_seq", "rev_seq") if p.get(k)]
    sites = primer_sites(seqs, db, cfg)
    max_product = int(cfg["scoring"]["max_product"])
    out = {}
    for p in primers:
        fwd = sites.get((p.get("fwd_seq") or "").upper(), [])
        rev = sites.get((p.get("rev_seq") or "").upper(), [])
        out[p["id"]] = {
            "amplicons": pair_amplicons(fwd, rev, self_spans, max_product),
            "n_sites":   {"fwd": len(fwd), "rev": len(rev)},
        }
    _report(reporter, "Primer screen: {} pair(s) over {} distinct primer(s) in {:.1f}s".format(
        len(primers), len(set(s.upper() for s in seqs)), time.perf_counter() - t0),
        stage="genotyping_primers")
    return out


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(genomic_fasta, taxid, config=None, report=None):
    """CLI entry point: Screen A for one genomic FASTA, printing a summary."""
    cfg = ots.load_config(config)
    seq = ""
    with open(genomic_fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq += line.strip()
    result = run_region_screen(seq, "", taxid, cfg=cfg, report=report)
    print("database:        {}".format(result["database"]))
    print("duplicates:      {}".format(len(result["duplicates"])))
    print("identical:       {}".format(len(result["identical"])))
    print("excluded:        {}".format(len(result["excluded"])))
    print("below threshold: {}".format(result["n_below_threshold"]))
    return result


def _parse_args(argv=None):
    """Command-line arguments for standalone use."""
    p = argparse.ArgumentParser(description="Local BLAST+ off-target screens.")
    sub = p.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="build the BLAST database for a taxid")
    b.add_argument("--taxid", required=True)
    b.add_argument("--force", action="store_true")
    r = sub.add_parser("region", help="run Screen A on a genomic FASTA")
    r.add_argument("--genomic-fasta", required=True)
    r.add_argument("--taxid", required=True)
    for s in (b, r):
        s.add_argument("--config", default=None, help="alternate offtarget.config.json")
    return p.parse_args(argv)


if __name__ == "__main__":
    args = _parse_args()
    if args.cmd == "build":
        build_blastdb(args.taxid, ots.load_config(args.config), force=args.force)
    else:
        main(args.genomic_fasta, args.taxid, config=args.config)
