"""
conservation_remote.py

EBI BLAST + dbfetch + Clustal Omega backend for ortholog conservation scoring
— moved out of blast_orthologs.py verbatim as part of the local/remote backend
split (see scripts/providers.py). This is the default backend; a future
conservation_local.py will search a local DIAMOND-indexed ortholog reference
set and align with MAFFT instead of submitting queued EBI jobs, emitting the
same .aln/.jsd/.isoforms.json files consumed by the rest of the pipeline.

Matt Rich, 9/2024 / updated 2026 — EBI REST calls via ebi_rest.py; backend-split 2026
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from site_selection_util import read_fasta, resolve_taxids

# ensure scripts/ is importable when called from the app
sys.path.insert(0, str(Path(__file__).parent))
import ebi_rest
from progress import report as _report, resolve_reporter, timed_poll_adapter
from derive_isoforms import derive_isoforms


def main(fasta_in, email, workingdir, name, output,
         n, evalue, db, length_percent,
         align_full_seqs, taxid, clients_folder, exclude_paralogs,
         taxid_file=None, min_perfect_len=40,
         report=None, job_id_cb=None, resume_job_ids=None):
    """Run BLAST → filter hits → fetch full seqs → clustalo → JSD scoring.

    job_id_cb(index, jid), when given, is called with the EBI job ID for the
    (single) BLAST submission this task makes, tagged as index 0.
    resume_job_ids, when given, is the job ID list persisted from a previous
    attempt; if index 0 holds an ID, reattach to it via ebi_rest.resume_job()
    instead of resubmitting. Returns {"ebi_status": "pending"|"expired", ...}
    if the resumed job hasn't finished, instead of raising.
    taxid_file, when given, supplies additional taxids (one per line) merged
    with the manually-entered taxid string.
    """
    # imported lazily (not at module top) to avoid a load-time circular import —
    # blast_orthologs.py imports providers.py, which imports this module on demand
    from blast_orthologs import (
        group_hits_by_species, parse_dbfetch_fasta_record,
        format_seq_label, ensure_query_in_alignment_set,
    )

    taxid = resolve_taxids(taxid, taxid_file)
    reporter = resolve_reporter(report)

    seq_name, seq = read_fasta(fasta_in)
    seq_len = float(len(seq))

    out_prefix = str(Path(output).with_suffix(""))

    ############################
    # BLAST TO FIND ORTHOLOGS
    ###########################

    resume_id = (resume_job_ids or [None])[0]
    if resume_id:
        _report(reporter, "Checking previously-submitted NCBI BLAST job…", stage="blast_submit")
        state, payload = ebi_rest.resume_job(ebi_rest.NCBIBLAST, resume_id, "json")
        if state == "pending":
            return {"ebi_status": "pending", "detail": payload}
        if state == "expired":
            return {"ebi_status": "expired", "detail": payload}
        blast_json_bytes = payload
    else:
        blast_params = {
            "email":      email,
            "program":    "blastp",
            "stype":      "protein",
            "sequence":   str(seq),  # Biopython Seq objects must be coerced to str
            "database":   db,
            "outformat":  "json",
            "alignments": n,         # must match scores; 0 produces empty hits list
            "scores":     n,
            "exp":        ebi_rest.fmt_exp(evalue),
        }
        if taxid:
            blast_params["taxids"] = taxid

        _report(reporter, "Submitting NCBI BLAST job…", stage="blast_submit")
        poll_cb = ebi_rest.combined_poll_cb(
            ebi_rest.indexed_job_id_cb(job_id_cb, 0),
            timed_poll_adapter(reporter, stage="blast_submit"),
        )
        blast_job_id = ebi_rest.run_job(ebi_rest.NCBIBLAST, blast_params, poll_cb=poll_cb)

        # fetch JSON result and save to the path the rest of the script expects
        blast_json_bytes = ebi_rest.fetch_result(ebi_rest.NCBIBLAST, blast_job_id, "json")

    blast_json_path = f"{out_prefix}.json.json"
    with open(blast_json_path, "wb") as f:
        f.write(blast_json_bytes)

    ###########################
    # PROCESS BLAST OUTPUT
    ###########################

    with open(blast_json_path) as f:
        blast_output = json.load(f)

    qlen = blast_output["query_len"]
    if not blast_output.get("hits"):
        _report(reporter, "BLAST returned no hits — check the sequence and search parameters.",
                stage="blast", level="error")
        return
    query_species = blast_output["hits"][0]["hit_os"]

    # derive per-isoform coordinate mappings (UniProt-first, BLAST fallback) while the
    # full BLAST JSON is in memory
    _report(reporter, "Detecting isoforms…", stage="isoforms")
    iso_result = derive_isoforms(blast_output, min_perfect_len=min_perfect_len)
    isoform_path = f"{out_prefix}.isoforms.json"
    with open(isoform_path, "w") as _f:
        json.dump(iso_result, _f)
    n_iso = len(iso_result["isoforms"])
    if n_iso:
        _report(reporter, f"Found {n_iso} isoform(s) (source: {iso_result['source']}).",
                stage="isoforms")
    else:
        _report(reporter, "No additional isoforms detected (single-isoform gene).",
                stage="isoforms")

    blast_hits = group_hits_by_species(blast_output["hits"], query_species, len(seq),
                                       evalue, length_percent, exclude_paralogs, n)

    ###########################
    # GET FULL SEQS OF HITS
    ###########################

    input_match = ""
    fasta_str_list = []
    total_hits = sum(len(v) for v in blast_hits.values())
    # flat ordered list of (species, hit) preserving blast_hits iteration order
    ordered_hits = [(s, h) for s in blast_hits for h in blast_hits[s]]

    if not align_full_seqs:
        # sequences come directly from the BLAST result — no network call needed
        for _, h in ordered_hits:
            fasta_str_list.append(">{}\n{}\n".format(
                format_seq_label(h["acc"], h["species"]), h["hitseq"]))
    else:
        accs = [h["acc"] for _, h in ordered_hits]
        _WORKERS = 10
        _report(reporter,
                f"fetching {total_hits} sequences from UniProt "
                f"({min(_WORKERS, total_hits)} concurrent requests)…",
                stage="dbfetch")

        def _fetch_one(acc):
            """Fetch one UniProt FASTA; return (acc, header_word, seq_str) or (acc, None, None)."""
            try:
                raw = ebi_rest.dbfetch("uniprotkb", acc, "fasta", "raw").decode("utf-8")
                name, seq = parse_dbfetch_fasta_record(raw)
                return acc, name, seq
            except Exception as e:
                return acc, None, str(e)

        fasta_dict = {}   # acc → (header_word, seq_str)
        with ThreadPoolExecutor(max_workers=_WORKERS) as pool:
            futures = {pool.submit(_fetch_one, acc): acc for acc in accs}
            done = 0
            for future in as_completed(futures):
                done += 1
                acc, name, result = future.result()
                if name is None:
                    _report(reporter, f"dbfetch failed for {acc}: {result}; skipping",
                            stage="dbfetch", level="warning")
                else:
                    fasta_dict[acc] = (name, result)
                if done % 10 == 0 or done == total_hits:
                    _report(reporter, f"dbfetch {done}/{total_hits}", stage="dbfetch")

        for _, h in ordered_hits:
            acc = h["acc"]
            if acc not in fasta_dict:
                continue
            _, acc_fasta_seq = fasta_dict[acc]
            label = format_seq_label(acc, h["species"])

            if acc_fasta_seq == seq:
                _report(reporter, f"found match to input: {label}", stage="dbfetch")
                input_match = label

            hit_len = float(len(acc_fasta_seq))
            if length_percent != 0:
                if hit_len / seq_len < length_percent or hit_len / seq_len > 1 / length_percent:
                    continue

            fasta_str_list.append(">{}\n{}\n".format(label, acc_fasta_seq.rstrip("\n")))

    # if no exact match to input, insert our query sequence at the front
    if input_match == "":
        _report(reporter, "no exact match to input; using best BLAST hit as query", stage="dbfetch")
    fasta_str_list, input_match = ensure_query_in_alignment_set(
        fasta_str_list, input_match, seq_name, seq)

    ###########################
    # ALIGN WITH CLUSTAL OMEGA (via EBI REST)
    ###########################

    # write the multi-sequence FASTA that will be submitted to clustalo
    fasta_out_path = f"{out_prefix}.fasta"
    with open(fasta_out_path, "w") as fasta_out:
        fasta_out.write("".join(fasta_str_list))

    aln_path = f"{out_prefix}.aln"

    if len(fasta_str_list) <= 1:
        # clustalo requires ≥2 sequences; fall back to copying input as alignment
        _report(reporter, "only one sequence; skipping Clustal, copying input as .aln", stage="align")
        with open(aln_path, "w") as f:
            f.write("".join(fasta_str_list))
    else:
        with open(fasta_out_path) as f:
            fasta_content = f.read()

        clustalo_params = {
            "email":    email,
            "sequence": fasta_content,
            "stype":    "protein",
            "outfmt":   "fa",
        }
        _report(reporter, "Submitting Clustal Omega job…", stage="align_submit")
        clustalo_job_id = ebi_rest.run_job(ebi_rest.CLUSTALO, clustalo_params,
                                           poll_cb=timed_poll_adapter(reporter, stage="align_submit"))
        aln_bytes = ebi_rest.fetch_result(ebi_rest.CLUSTALO, clustalo_job_id, "aln-fasta")
        with open(aln_path, "wb") as f:
            f.write(aln_bytes)
        _report(reporter, f"Alignment written → {aln_path}", stage="align")

    ###########################
    # RENDER ALIGNMENT IMAGE
    ###########################

    # count sequences for progress context before the slow render step
    try:
        _aln_seqs = sum(1 for ln in open(aln_path) if ln.startswith(">"))
        _aln_len  = next((len(ln.rstrip()) for ln in open(aln_path) if not ln.startswith(">")), 0)
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
            align_file         = aln_path,
            window_size        = 3,
            win_lam            = 0.5,
            outfile_name       = jsd_out,
            s_matrix_file      = f"{_scripts}matrix/blosum62.bla",
            bg_distribution    = blosum_bg,
            scoring_function   = sc.js_divergence,
            use_seq_weights    = True,
            gap_cutoff         = 0.75,
            use_gap_penalty    = True,
            seq_specific_output= best_hit_name.rstrip(),
            normalize_scores   = False,
        )
    return 0
