"""
blast_orthologs.py

Matt Rich, 9/2024 / updated 2026 — EBI REST calls via ebi_rest.py
"""

import json
import os
import sys

from pathlib import Path
from site_selection_util import read_fasta

# ensure scripts/ is importable when called from the app
sys.path.insert(0, str(Path(__file__).parent))
import ebi_rest
from progress import report as _report, resolve_reporter, poll_adapter
from derive_isoforms import derive_isoform_segments


def hit_to_dict(j):
    """Convert one JSON BLAST hit to a compact dict."""
    return {
        "acc":     j["hit_acc"],
        "species": j["hit_os"],
        "evalue":  float(j["hit_hsps"][0]["hsp_expect"]),
        "identity": float(j["hit_hsps"][0]["hsp_identity"]),
        "length":  int(j["hit_hsps"][0]["hsp_hit_to"]) - int(j["hit_hsps"][0]["hsp_hit_from"]),
        "hitseq":  j["hit_hsps"][0]["hsp_hseq"].replace("-", ""),
    }


def main(fasta_in, email, workingdir, name, output,
         n, evalue, db, length_percent,
         align_full_seqs, taxid, clients_folder, exclude_paralogs,
         report=None):
    """Run BLAST → filter hits → fetch full seqs → clustalo → JSD scoring."""
    reporter = resolve_reporter(report)

    seq_name, seq = read_fasta(fasta_in)
    seq_len = float(len(seq))

    out_prefix = str(Path(output).with_suffix(""))

    ############################
    # BLAST TO FIND ORTHOLOGS
    ###########################

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
    if str(taxid) not in ("", "1", "1.0"):
        blast_params["taxids"] = str(taxid)

    _report(reporter, "Submitting NCBI BLAST job…", stage="blast_submit")
    blast_job_id = ebi_rest.run_job(ebi_rest.NCBIBLAST, blast_params, poll_cb=poll_adapter(reporter))

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

    # derive isoform coverage segments from same-organism hits while the full JSON is in memory
    _report(reporter, "Detecting isoforms from same-organism hits…", stage="isoforms")
    iso_segs = derive_isoform_segments(blast_output)
    isoform_path = f"{out_prefix}.isoforms.tsv"
    with open(isoform_path, "w") as _f:
        for _start, _stop, _desc in iso_segs:
            _f.write(f"isoforms\t{_start}\t{_stop}\t{_desc}\n")
    if iso_segs:
        _report(reporter, f"Found {len(iso_segs)} isoform segment(s) across "
                f"{iso_segs[-1][1]} residues.", stage="isoforms")
    else:
        _report(reporter, "No additional isoforms detected (single-isoform gene).",
                stage="isoforms")

    blast_hits = {query_species: []}
    for h in blast_output["hits"]:
        d_hit = hit_to_dict(h)

        # apply evalue and length filters
        if length_percent != 0:
            if d_hit["evalue"] > evalue:
                continue
            if (d_hit["length"] / float(len(seq)) < length_percent
                    or d_hit["length"] / float(len(seq)) > 1 / length_percent):
                continue

        if exclude_paralogs:
            # store all hits from query species; best hit per other species
            if d_hit["species"] == query_species:
                blast_hits[query_species].append(d_hit)
            else:
                if d_hit["species"] not in blast_hits:
                    blast_hits[d_hit["species"]] = [d_hit]
                elif blast_hits[d_hit["species"]][0]["evalue"] > d_hit["evalue"]:
                    blast_hits[d_hit["species"]] = [d_hit]
        else:
            if d_hit["species"] in blast_hits:
                blast_hits[d_hit["species"]].append(d_hit)
            else:
                blast_hits[d_hit["species"]] = [d_hit]

        # stop once we have enough species
        if len(blast_hits) >= n:
            break

    ###########################
    # GET FULL SEQS OF HITS
    ###########################

    input_match = ""
    fasta_str_list = []
    total_hits = sum(len(v) for v in blast_hits.values())
    fetched = 0
    _report(reporter, f"fetching {total_hits} sequences from UniProt…", stage="dbfetch")
    for s in blast_hits:
        for h in blast_hits[s]:
            if align_full_seqs:
                fetched += 1
                _report(reporter, f"dbfetch {fetched}/{total_hits}: {h['acc']}",
                        stage="dbfetch")
                try:
                    acc_fasta_bytes = ebi_rest.dbfetch("uniprotkb", h["acc"], "fasta", "raw")
                    acc_fasta = acc_fasta_bytes.decode("utf-8")
                except Exception as e:
                    _report(reporter, f"dbfetch failed for {h['acc']}: {e}; skipping",
                            stage="dbfetch", level="warning")
                    continue

                acc_fasta_name = acc_fasta.split("\n")[0].split()[0]
                acc_fasta_seq = "".join(acc_fasta.split("\n")[1:])

                if acc_fasta_seq == seq:
                    _report(reporter, f"found match to input: {acc_fasta.split(chr(10))[0]}",
                            stage="dbfetch")
                    input_match = acc_fasta_name[1:]

                # filter by total length
                hit_len = float(len("".join(acc_fasta.split("\n")[1:])))
                if length_percent != 0:
                    if hit_len / seq_len < length_percent or hit_len / seq_len > 1 / length_percent:
                        continue

                fasta_str_list.append(acc_fasta_name + "\n" + acc_fasta_seq.rstrip("\n") + "\n")
            else:
                fasta_str_list.append(">{}\n{}\n".format(h["acc"].split()[0], h["hitseq"]))

    # if no exact match to input, insert our query sequence at the front
    if input_match == "":
        _report(reporter, "no exact match to input; using best BLAST hit as query", stage="dbfetch")
        fasta_str_list = [">{}\n{}\n".format(seq_name, seq)] + fasta_str_list[1:]
        input_match = seq_name

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
                                           poll_cb=poll_adapter(reporter))
        aln_bytes = ebi_rest.fetch_result(ebi_rest.CLUSTALO, clustalo_job_id, "aln-fasta")
        with open(aln_path, "wb") as f:
            f.write(aln_bytes)

    ###########################
    # RENDER ALIGNMENT IMAGE
    ###########################

    # import in-process (no subprocess needed; build_heatmap_images.py is a local script)
    try:
        import build_heatmap_images
        build_heatmap_images.plot_alignment_matrix_matplotlib(aln_path, "svg", "")
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


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument("-f", "--fasta", "--input_file", action="store", type=str, dest="FASTA_IN",
        help="name of fasta file containing seq.", required=True)
    parser.add_argument("--email", action="store", type=str, dest="EMAIL",
        help="email address, required by EBI job submission.", required=True)
    parser.add_argument("--dir", "--working_dir", action="store", type=str, dest="WORKINGDIR",
        help="working directory for output", required=True)
    parser.add_argument("--name", "--run_name", action="store", type=str, dest="NAME",
        help="prefix name for output", required=True)
    parser.add_argument("--output", action="store", type=str, dest="OUTPUT",
        help="user-supplied output filename", default=None)

    parser.add_argument("--taxid", action="store", type=str, dest="TAXID",
        help="taxid to use for blast search", default=1)
    parser.add_argument("--taxid_file", action="store", type=str, dest="TAX_FILE",
        help="file containing taxids for BLAST, one per line", default=None)
    parser.add_argument("-e", "--evalue", action="store", type=float, dest="EVALUE",
        help="evalue threshold for keeping hits (1e-10)", default=1e-10)
    parser.add_argument("-n", "--max_hits", action="store", type=int, dest="MAX_HITS",
        help="max number of hits to keep (100)", default=100)
    parser.add_argument("--db", action="store", type=str, dest="DB",
        help="Uniprot database to search (uniprotkb)", default="uniprotkb")
    parser.add_argument("-l", "--length", action="store", type=float, dest="LENGTH",
        help="minimum match length as percent of input (0)", default=0)
    parser.add_argument("--align-blast-sequence", action="store_false", dest="FULLSEQS",
        help="only align BLAST hit sequences, not full UniProt seqs", default=True)
    parser.add_argument("--clients-folder", action="store", type=str, dest="CLIENTS_FOLDER",
        help="(unused; retained for CLI compatibility)", default="./scripts/")
    parser.add_argument("--exclude-paralogs", action="store_true", dest="EXCLUDE_PARALOGS",
        help="return only best match per non-subject species", default=False)

    args, unknowns = parser.parse_known_args()

    main(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.NAME, args.OUTPUT,
         args.MAX_HITS, args.EVALUE, args.DB, args.LENGTH, args.FULLSEQS,
         args.TAXID, args.CLIENTS_FOLDER, args.EXCLUDE_PARALOGS)
