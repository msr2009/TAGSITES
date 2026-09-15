"""
conservation_local.py

Local DIAMOND + MAFFT backend for ortholog conservation scoring — searches
one or more local DIAMOND databases built by scripts/reference_data.py
(default: full Swiss-Prot plus all unreviewed Nematoda TrEMBL — see its
fetch_swissprot()/fetch_rhabditida_trembl()) instead of submitting an EBI BLAST
job, aligns with MAFFT instead of a Clustal Omega job, then feeds the
alignment to the same score_conservation_py3.py used by conservation_remote.py
— so downstream code (utils/results.py, utils/scoring.py) sees the exact same
.aln/.jsd/.isoforms.json files regardless of which backend produced them.

The set of databases to search — and each one's own evalue/max_target_seqs/
sensitivity — is config-driven (batch.config.json's conservation_local.
search_databases list), not hardcoded, so a database can be added, removed,
or retuned without touching this file. Every search runs with
--ultra-sensitive by default (DIAMOND's closest available approximation to
real NCBI BLASTP sensitivity — conservation_remote.py submits an actual
BLASTP job, so this backend shouldn't quietly trade sensitivity for DIAMOND's
faster default mode) unless a database entry overrides "sensitivity".

Reuses blast_orthologs.py's pure helpers (group_hits_by_species,
format_seq_label, ensure_query_in_alignment_set) and derive_isoforms.py
unchanged, by shaping DIAMOND's TSV output into the same "hits" list schema
EBI's BLAST JSON uses — no changes needed to either of those modules to
support a second conservation backend, or to support searching more than one
database. Hits from every configured database are merged and re-sorted by
ascending e-value before being handed off, since group_hits_by_species()/
derive_isoforms() both assume the hit list arrives best-first (e.g.
derive_isoforms()'s "first hit is the query's own species" heuristic) — a
naive concatenation across databases would silently violate that whenever a
later database's best hit actually outranks an earlier database's.

Score/rank differences from conservation_remote.py are still expected
(DIAMOND is not NCBI BLAST, MAFFT is not Clustal Omega) — see the plan's
Verification step 1 for the parity bar (rank correlation, not byte equality).
"""

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from site_selection_util import read_fasta

sys.path.insert(0, str(Path(__file__).parent))
from progress import report as _report, resolve_reporter
from providers import _load_config as _load_batch_config
from derive_isoforms import derive_isoforms

_STITLE_OS_RE = re.compile(r"OS=(.*?)\s+OX=")

# default search-database list when batch.config.json has no conservation_local
# key at all — keeps an unconfigured setup working (Swiss-Prot only) rather
# than erroring outright
_DEFAULT_SEARCH_DATABASES = [{"name": "swissprot", "evalue": 1e-10, "max_target_seqs": 50}]

_DEFAULT_SENSITIVITY = "ultra-sensitive"


def _reference_dir(cfg=None):
    cfg = cfg or _load_batch_config().get("reference_data", {})
    out_dir = Path(cfg.get("out_dir", "data/reference"))
    if not out_dir.is_absolute():
        out_dir = Path(__file__).parent.parent / out_dir
    return out_dir


def _load_search_databases():
    """Return the configured list of search-database dicts, each with a
    resolved 'path' key added (out_dir/{name}.dmnd). Raises FileNotFoundError
    naming any database file that doesn't exist yet and the fetch command to
    build it, so a missing/mistyped database fails loudly rather than
    silently searching fewer databases than configured.
    """
    ref_dir = _reference_dir()
    entries = _load_batch_config().get("conservation_local", {}).get(
        "search_databases", _DEFAULT_SEARCH_DATABASES
    )
    resolved = []
    for entry in entries:
        entry = dict(entry)
        path = ref_dir / f"{entry['name']}.dmnd"
        if not path.exists():
            raise FileNotFoundError(
                f"{path} not found — run `python scripts/reference_data.py "
                f"--only {entry['name']}` first."
            )
        entry["path"] = path
        entry.setdefault("sensitivity", _DEFAULT_SENSITIVITY)
        resolved.append(entry)
    return resolved


def _parse_sseqid(sseqid):
    """UniProt FASTA headers are 'sp|ACCESSION|ENTRYNAME' or 'tr|ACCESSION|...';
    DIAMOND's sseqid is the first whitespace-delimited token, i.e. that whole
    string. Return the bare accession (middle pipe field), or the raw sseqid
    if it doesn't look like that format.
    """
    parts = sseqid.split("|")
    return parts[1] if len(parts) >= 3 else sseqid


def _parse_species(stitle):
    """Extract the OS=... organism name out of a UniProt FASTA description line."""
    m = _STITLE_OS_RE.search(stitle)
    return m.group(1) if m else ""


def _run_diamond_blastp_one(query_fasta, db_entry, evalue, workdir):
    """Run `diamond blastp` against one configured database and return a list
    of hit dicts in the same raw shape blast_orthologs.hit_to_dict() expects
    (i.e. one EBI-BLAST-JSON hit record each), each tagged with a "source_db"
    key naming the database — a harmless extra field existing consumers
    (hit_to_dict(), derive_isoforms()) simply don't read.
    """
    out_tsv = Path(workdir) / f"diamond_{db_entry['name']}.tsv"
    fields = ["sseqid", "stitle", "pident", "evalue", "qstart", "qend",
              "sstart", "send", "full_sseq"]
    max_target_seqs = db_entry.get("max_target_seqs", 50)
    db_evalue = db_entry.get("evalue", evalue)
    sensitivity = db_entry.get("sensitivity", _DEFAULT_SENSITIVITY)
    subprocess.run([
        "diamond", "blastp",
        "--query", str(query_fasta),
        "--db", str(db_entry["path"]),
        "--out", str(out_tsv),
        "--outfmt", "6", *fields,
        "--max-target-seqs", str(max_target_seqs),
        "--evalue", str(db_evalue),
        f"--{sensitivity}",
        "--quiet",
    ], check=True)

    hits = []
    with open(out_tsv) as f:
        for line in f:
            sseqid, stitle, pident, ev, qstart, qend, sstart, send, full_sseq = line.rstrip("\n").split("\t")
            hits.append({
                "hit_acc": _parse_sseqid(sseqid),
                "hit_os": _parse_species(stitle),
                "source_db": db_entry["name"],
                "hit_hsps": [{
                    "hsp_expect": ev,
                    "hsp_identity": pident,
                    "hsp_hit_to": send,
                    "hsp_hit_from": sstart,
                    "hsp_hseq": full_sseq,
                }],
            })
    return hits


def _run_diamond_blastp(query_fasta, db_entries, n, evalue, workdir):
    """Run diamond blastp against every configured database, merge the hit
    lists, and re-sort by ascending e-value — group_hits_by_species()/
    derive_isoforms() both assume the incoming hit list is already best-first
    (see this module's docstring), which only holds automatically within a
    single database's own output, not across several concatenated together.
    """
    all_hits = []
    for db_entry in db_entries:
        all_hits.extend(_run_diamond_blastp_one(query_fasta, db_entry, evalue, workdir))
    all_hits.sort(key=lambda h: float(h["hit_hsps"][0]["hsp_expect"]))
    return all_hits


def main(fasta_in, email, workingdir, name, output,
         n, evalue, db, length_percent,
         align_full_seqs, taxid, clients_folder, exclude_paralogs,
         taxid_file=None, min_perfect_len=40,
         report=None, job_id_cb=None, resume_job_ids=None):
    """DIAMOND blastp -> filter hits -> MAFFT align -> JSD scoring.

    Same signature as conservation_remote.main(); email/db/taxid/taxid_file/
    job_id_cb/resume_job_ids are accepted but unused (no EBI job is submitted,
    so there's nothing to tag or resume, and no per-database/taxid restriction
    beyond the fixed reference set this backend searches). align_full_seqs is
    also accepted but unused — full subject sequences come from the local
    reference DB at no extra cost, so there's no cheap/expensive path to
    choose between here the way there is for a per-hit EBI dbfetch.
    """
    from blast_orthologs import group_hits_by_species, format_seq_label, ensure_query_in_alignment_set

    reporter = resolve_reporter(report)
    seq_name, seq = read_fasta(fasta_in)
    seq_len = float(len(seq))
    out_prefix = str(Path(output).with_suffix(""))
    db_entries = _load_search_databases()

    ###########################
    # DIAMOND SEARCH
    ###########################

    db_names = ", ".join(e["name"] for e in db_entries)
    _report(reporter, f"Searching local reference database(s) ({db_names}) with DIAMOND…",
            stage="blast_submit")
    with tempfile.TemporaryDirectory() as tmpdir:
        query_fasta = Path(tmpdir) / "query.fa"
        query_fasta.write_text(f">{seq_name}\n{seq}\n")
        raw_hits = _run_diamond_blastp(query_fasta, db_entries, n, evalue, tmpdir)

    blast_output = {"query_len": int(seq_len), "hits": raw_hits}
    blast_json_path = f"{out_prefix}.json.json"
    with open(blast_json_path, "w") as f:
        json.dump(blast_output, f)

    if not raw_hits:
        _report(reporter, "DIAMOND returned no hits — check the sequence and search parameters.",
                stage="blast", level="error")
        return
    query_species = raw_hits[0]["hit_os"]

    _report(reporter, "Detecting isoforms…", stage="isoforms")
    iso_result = derive_isoforms(blast_output, min_perfect_len=min_perfect_len)
    isoform_path = f"{out_prefix}.isoforms.json"
    with open(isoform_path, "w") as _f:
        json.dump(iso_result, _f)
    n_iso = len(iso_result["isoforms"])
    if n_iso:
        _report(reporter, f"Found {n_iso} isoform(s) (source: {iso_result['source']}).", stage="isoforms")
    else:
        _report(reporter, "No additional isoforms detected (single-isoform gene).", stage="isoforms")

    blast_hits = group_hits_by_species(raw_hits, query_species, len(seq),
                                       evalue, length_percent, exclude_paralogs, n)

    ###########################
    # BUILD ALIGNMENT SET
    ###########################

    input_match = ""
    fasta_str_list = []
    ordered_hits = [(s, h) for s in blast_hits for h in blast_hits[s]]
    for _, h in ordered_hits:
        label = format_seq_label(h["acc"], h["species"])
        hit_seq = h["hitseq"]
        if hit_seq == seq:
            _report(reporter, f"found match to input: {label}", stage="dbfetch")
            input_match = label
        fasta_str_list.append(">{}\n{}\n".format(label, hit_seq))

    if input_match == "":
        _report(reporter, "no exact match to input; using best hit as query", stage="dbfetch")
    fasta_str_list, input_match = ensure_query_in_alignment_set(
        fasta_str_list, input_match, seq_name, seq)

    ###########################
    # ALIGN WITH MAFFT (local)
    ###########################

    fasta_out_path = f"{out_prefix}.fasta"
    with open(fasta_out_path, "w") as fasta_out:
        fasta_out.write("".join(fasta_str_list))

    aln_path = f"{out_prefix}.aln"
    if len(fasta_str_list) <= 1:
        _report(reporter, "only one sequence; skipping MAFFT, copying input as .aln", stage="align")
        with open(aln_path, "w") as f:
            f.write("".join(fasta_str_list))
    else:
        _report(reporter, "Running MAFFT…", stage="align_submit")
        with open(aln_path, "w") as aln_out:
            subprocess.run(["mafft", "--auto", "--quiet", fasta_out_path],
                           check=True, stdout=aln_out)
        _report(reporter, f"Alignment written → {aln_path}", stage="align")

    ###########################
    # RENDER ALIGNMENT IMAGE
    ###########################

    try:
        _aln_seqs = sum(1 for ln in open(aln_path) if ln.startswith(">"))
        _aln_len = next((len(ln.rstrip()) for ln in open(aln_path) if not ln.startswith(">")), 0)
        _report(reporter,
                f"Rendering alignment image ({_aln_seqs} sequences × {_aln_len} positions)…",
                stage="align_img")
    except Exception:
        _report(reporter, "Rendering alignment image…", stage="align_img")

    try:
        import build_heatmap_reportlab
        build_heatmap_reportlab.plot_alignment_reportlab(aln_path)
        _report(reporter, "Alignment image saved.", stage="align_img")
    except Exception as e:
        _report(reporter, f"alignment image generation failed: {e}", stage="align_img", level="warning")

    ###########################
    # CALCULATE JSD
    ###########################

    import score_conservation_py3 as sc
    best_hit_name = input_match
    jsd_path = f"{out_prefix}.jsd"
    _scripts = str(Path(__file__).parent) + "/"
    blosum_bg = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025,
                 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]
    _report(reporter, f"Scoring conservation: {aln_path} → {jsd_path}", stage="score")
    with open(jsd_path, "w") as jsd_out:
        sc.main(
            align_file=aln_path,
            window_size=3,
            win_lam=0.5,
            outfile_name=jsd_out,
            s_matrix_file=f"{_scripts}matrix/blosum62.bla",
            bg_distribution=blosum_bg,
            scoring_function=sc.js_divergence,
            use_seq_weights=True,
            gap_cutoff=0.75,
            use_gap_penalty=True,
            seq_specific_output=best_hit_name.rstrip(),
            normalize_scores=False,
        )
    return 0
