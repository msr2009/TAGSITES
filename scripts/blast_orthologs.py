"""
blast_orthologs.py

CLI entry point for ortholog conservation scoring. The actual search +
alignment is delegated to a backend resolved by scripts/providers.py —
scripts/conservation_remote.py (EBI BLAST + dbfetch + Clustal Omega) by
default, or a local DIAMOND/MAFFT lookup when batch.config.json sets
backends.conservation to a local mode. With no config file present this
always resolves to "remote", so the Shiny app's behavior here is unchanged.

Matt Rich, 9/2024 / updated 2026 — EBI REST calls via ebi_rest.py; backend-split 2026
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from providers import resolve


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


def group_hits_by_species(raw_hits, query_species, seq_len, evalue, length_percent,
                          exclude_paralogs, n):
    """Filter raw BLAST hits by evalue/length, then group by species (best hit per species
    when exclude_paralogs=True; all matching hits otherwise). Stops once n species are seen.
    """
    blast_hits = {query_species: []}
    for h in raw_hits:
        d_hit = hit_to_dict(h)

        if length_percent != 0:
            if d_hit["evalue"] > evalue:
                continue
            if (d_hit["length"] / float(seq_len) < length_percent
                    or d_hit["length"] / float(seq_len) > 1 / length_percent):
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

        if len(blast_hits) >= n:
            break

    return blast_hits


def parse_dbfetch_fasta_record(raw_text):
    """Parse one raw dbfetch FASTA response into (header_word, sequence)."""
    name = raw_text.split("\n")[0].split()[0]
    seq  = "".join(ln for ln in raw_text.split("\n") if not ln.startswith(">"))
    return name, seq


def format_seq_label(acc, species):
    """Build a FASTA header/label combining UniProt accession and species,
    e.g. "P12345_H_sapiens" — always keeps the accession so sequences stay
    traceable even when the species string is missing or unusual. Genus is
    abbreviated and anything past the binomial (strain/collection codes,
    which can run to 100+ chars in UniProt's hit_os field) is dropped.
    """
    acc = acc.split()[0]
    parts = (species or "").split()
    if len(parts) >= 2:
        binomial = f"{parts[0][0]}_{parts[1]}"
    else:
        binomial = parts[0] if parts else ""
    species_clean = re.sub(r"[^A-Za-z0-9]+", "_", binomial).strip("_")
    return f"{acc}_{species_clean}" if species_clean else acc


def ensure_query_in_alignment_set(fasta_str_list, input_match, seq_name, seq):
    """If no hit exactly matched the query, insert the query itself at the front.

    Returns (fasta_str_list, input_match) — the alignment set Clustal Omega receives
    must contain the query sequence so conservation scoring has a reference row.
    """
    if input_match == "":
        fasta_str_list = [">{}\n{}\n".format(seq_name, seq)] + fasta_str_list[1:]
        input_match = seq_name
    return fasta_str_list, input_match


def main(fasta_in, email, workingdir, name, output,
         n, evalue, db, length_percent,
         align_full_seqs, taxid, clients_folder, exclude_paralogs,
         taxid_file=None, min_perfect_len=40,
         report=None, job_id_cb=None, resume_job_ids=None):
    """Run ortholog conservation scoring via the configured backend (EBI BLAST
    + Clustal Omega by default: BLAST → filter hits → fetch full seqs →
    clustalo → JSD scoring); same signature/return value as before the
    backend split.

    job_id_cb(index, jid), when given, is called with the EBI job ID for the
    (single) BLAST submission the remote backend makes, tagged as index 0.
    resume_job_ids, when given, is the job ID list persisted from a previous
    attempt; if index 0 holds an ID, reattach to it via ebi_rest.resume_job()
    instead of resubmitting. Returns {"ebi_status": "pending"|"expired", ...}
    if the resumed job hasn't finished, instead of raising.
    taxid_file, when given, supplies additional taxids (one per line) merged
    with the manually-entered taxid string.
    """
    backend_main = resolve("conservation")
    return backend_main(fasta_in, email, workingdir, name, output,
                         n, evalue, db, length_percent,
                         align_full_seqs, taxid, clients_folder, exclude_paralogs,
                         taxid_file=taxid_file, min_perfect_len=min_perfect_len,
                         report=report, job_id_cb=job_id_cb, resume_job_ids=resume_job_ids)


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
    parser.add_argument("--min_perfect_len", action="store", type=int, dest="MIN_PERFECT_LEN",
        help="min consecutive perfect-match run to count a BLAST hit as an isoform (40)",
        default=40)

    args, unknowns = parser.parse_known_args()

    main(args.FASTA_IN, args.EMAIL, args.WORKINGDIR, args.NAME, args.OUTPUT,
         args.MAX_HITS, args.EVALUE, args.DB, args.LENGTH, args.FULLSEQS,
         args.TAXID, args.CLIENTS_FOLDER, args.EXCLUDE_PARALOGS,
         taxid_file=args.TAX_FILE, min_perfect_len=args.MIN_PERFECT_LEN)
