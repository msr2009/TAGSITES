"""setup_server.py — Shiny module server for the Setup tab.

Thin reactive layer: reads user inputs, delegates data construction to
setup_logic.py, writes the run JSON, and publishes its path via shared_json.
"""

import asyncio
import json
import os
import requests
import shutil
import sys
from pathlib import Path

from shiny import module, reactive, render, ui

from config import (
    DEFAULT_SPECIES, GLOBAL_TOOLTIPS, TASK_PARAMETERS,
    task_hidden, task_defaults, task_choices, task_tooltips, task_output_suffix,
    GLOBAL_DEFAULTS, GLOBAL_KEYS,
)
from modules.setup_ui import label_with_tip, TASK_DESCRIPTIONS
from modules.ui_helpers import compact_file_input
from modules.setup_logic import (
    build_defaults_entry,
    build_global_block,
    build_reagents_entry,
    build_run_json,
    build_task_entry,
    make_task,
)
from scripts.site_selection_util import get_sequence, save_fasta, extract_bfactors_from_pdb

# ebi_rest / fetch_genomic_sequence live in scripts/; make them importable
# regardless of working dir
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
import ebi_rest
from fetch_genomic_sequence import fetch_genomic_sequence

_ROOT = Path(__file__).parent.parent
_PARAMS_DIR = _ROOT / "params"
_TABLES_DIR = _ROOT / "tables"


def _table_choices(ext="*.tsv"):
    """Return {absolute_path_str: display_name} for files in tables/ matching ext."""
    return {str(f): f.stem for f in sorted(_TABLES_DIR.glob(ext))}


# ── widget builders ───────────────────────────────────────────────────────────

def _make_param_widget(widget_id, param, value, tip, choices=None):
    """Return one input widget for a task parameter."""
    lbl = label_with_tip(param.replace("_", " "), tip)
    if param == "scores_file":
        file_choices = _table_choices("*.tsv")
        selected = value if value in file_choices else next(iter(file_choices), "")
        return ui.input_select(widget_id, label=lbl, choices=file_choices, selected=selected)
    if choices:
        selected = value if value in choices else choices[0]
        return ui.input_select(widget_id, label=lbl, choices=choices, selected=selected, size=1)
    return ui.input_text(widget_id, label=lbl, value=str(value) if value != "" else "")


def _build_param_inputs(task, task_values_snap):
    """Build the list of param widgets for one task card."""
    tips       = task.get("tooltips", {})
    all_choices = task.get("choices", {})
    hidden     = task_hidden(task["name"])
    widgets = []
    for p, default in task["params"].items():
        if p in hidden:
            continue
        widget_id = f"{task['id']}_{p}"
        # use the snapshotted value if present (preserves edits across re-renders)
        current = task_values_snap.get(task["id"], {}).get(p, default)
        widgets.append(_make_param_widget(widget_id, p, current, tips.get(p, ""), all_choices.get(p)))
    return widgets


def _format_fasta(hit):
    """Return a FASTA-formatted string for a UniProt hit (60 chars per line)."""
    acc  = hit["accession"]
    gene = hit["gene"]
    org  = hit["organism"]
    iso  = hit.get("isoform_name", "")
    desc = f"{gene} isoform {iso}" if iso else gene
    header = f">{acc} {desc} [{org}]"
    seq = hit["sequence"]
    lines = [header] + [seq[i:i+60] for i in range(0, len(seq), 60)]
    return "\n".join(lines)


def _fetch_isoforms(canonical_acc, base_hit):
    """Expand a canonical UniProt entry into one hit dict per isoform.

    Reads the ALTERNATIVE PRODUCTS comment for the isoform list, then fetches
    each non-canonical isoform sequence from the per-accession FASTA endpoint.
    Returns a list of hit dicts (canonical first), or [base_hit] on failure.
    """
    try:
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{canonical_acc}",
            params={"format": "json"},
            timeout=15,
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return [base_hit]

    alt = next(
        (c for c in data.get("comments", []) if c.get("commentType") == "ALTERNATIVE PRODUCTS"),
        None,
    )
    if not alt or not alt.get("isoforms"):
        return [base_hit]

    hits = []
    for iso in alt["isoforms"]:
        iso_ids = iso.get("isoformIds", [])
        if not iso_ids:
            continue
        iso_acc    = iso_ids[0]                               # e.g. "P04637-2"
        iso_num    = iso_acc.rsplit("-", 1)[-1] if "-" in iso_acc else "1"
        synonyms   = [s["value"] for s in iso.get("synonyms", [])]
        iso_name   = synonyms[0] if synonyms else iso.get("name", {}).get("value", "")
        status     = iso.get("isoformSequenceStatus", "")
        safe_id    = iso_acc.replace("-", "_")                # safe Shiny input ID

        if status == "Displayed":
            # canonical isoform — sequence already in base_hit
            hit = dict(base_hit)
            hit.update(accession=iso_acc, isoform=iso_num,
                       isoform_name=iso_name, safe_id=safe_id)
        else:
            try:
                r2 = requests.get(
                    f"https://rest.uniprot.org/uniprotkb/{iso_acc}.fasta",
                    timeout=15,
                )
                if not r2.ok:
                    continue
                fasta_lines = r2.text.strip().splitlines()
                seq = "".join(fasta_lines[1:])
            except Exception:
                continue
            hit = dict(base_hit)
            hit.update(accession=iso_acc, sequence=seq, length=str(len(seq)),
                       isoform=iso_num, isoform_name=iso_name,
                       has_afdb=False, safe_id=safe_id)  # AFDB covers canonical only

        hits.append(hit)

    return hits if hits else [base_hit]


# ── module server ─────────────────────────────────────────────────────────────

@module.server
def setup_server(input, output, session, shared_json):

    tasks      = reactive.Value([])   # list of in-memory task dicts
    task_snap  = reactive.Value({})   # {task_id: {param: value}} — preserved across re-renders

    # prevents the self-save from re-loading tasks back into the setup form
    _skip_next_load = [False]

    organism_taxid   = reactive.Value("")   # resolved taxid string from organism selection
    _pdb_plddt_valid = reactive.Value(True) # False when uploaded PDB lacks valid pLDDT B-factors
    _search_hits     = reactive.Value([])   # [(name, taxid_str)] from UniProt taxonomy search

    # UniProt protein-search state
    _uniprot_hits = reactive.Value([])   # list of hit dicts from UniProt gene search
    _selected_hit = reactive.Value(None) # chosen hit dict, or None
    _selected_pdb = reactive.Value("")   # AFDB PDB text for the chosen hit ("" = none/unavailable)

    # ── organism selection ────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.organism)
    def _organism_changed():
        """Resolve taxid from preset dropdown; pre-fill existing blast task widgets."""
        name = input.organism()
        if not name or name not in DEFAULT_SPECIES:
            return
        taxid = DEFAULT_SPECIES[name]
        if taxid is None:   # "Other (search...)" sentinel — wait for search result
            _search_hits.set([])
            organism_taxid.set("")   # clear stale taxid from a previously selected preset
            return
        organism_taxid.set(str(taxid))

    @render.ui
    def organism_search_ui():
        """Show a search row only when 'Other (search…)' is selected."""
        if input.organism() != "Other (search...)":
            return None
        hits = _search_hits()
        result_widget = None
        if hits:
            choices = {"": "— select —",
                       **{str(taxid): f"{name} ({taxid})" for name, taxid in hits}}
            result_widget = ui.div(
                ui.input_select("organism_search_result", "", choices=choices),
                style="margin-top: 0.2rem;",
            )
        return ui.div(
            ui.div(
                ui.input_text("organism_search_text", "",
                              placeholder="Type species name…"),
                ui.input_action_button("organism_search_btn", "Search",
                                       class_="btn-sm btn-outline-secondary"),
                style="display:flex; gap:0.4rem; align-items:flex-end;",
            ),
            *([result_widget] if result_widget else []),
            style="margin-top: 0.2rem;",
        )

    @reactive.effect
    @reactive.event(input.organism_search_btn)
    def _organism_search():
        """Query UniProt taxonomy API and populate the result select."""
        query = input.organism_search_text().strip() if "organism_search_text" in input else ""
        if not query:
            return
        try:
            resp = requests.get(
                "https://rest.uniprot.org/taxonomy/search",
                params={"query": query, "format": "json", "size": 5},
                timeout=10,
            )
            resp.raise_for_status()
            results = resp.json().get("results", [])
            hits = [(r["scientificName"], r["taxonId"])
                    for r in results if "taxonId" in r]
            _search_hits.set(hits)
        except Exception as e:
            ui.notification_show(f"Taxonomy search failed: {e}", type="error", duration=5)

    @reactive.effect
    def _organism_search_result_changed():
        """Set taxid when user picks from search results."""
        if "organism_search_result" not in input:
            return
        val = input.organism_search_result()
        if not val:
            return
        organism_taxid.set(val)

    # ── taxonomic lineage ─────────────────────────────────────────────────────

    # ranks worth showing as BLAST scope options (clades and "no rank" filtered out)
    _BLAST_RANKS = {
        "domain", "superkingdom", "kingdom", "subkingdom",
        "phylum", "subphylum", "class", "subclass",
        "order", "suborder", "infraorder", "superfamily",
        "family", "subfamily", "genus", "species",
    }

    _lineage = reactive.Value([])  # [{"rank": str, "name": str, "taxid": str}, ...]

    @reactive.effect
    def _fetch_lineage():
        """Fetch taxonomic lineage from UniProt whenever the resolved taxid changes."""
        tid = organism_taxid()
        if not tid:
            _lineage.set([])
            return
        try:
            resp = requests.get(
                f"https://rest.uniprot.org/taxonomy/{tid}",
                params={"format": "json"},
                timeout=10,
            )
            resp.raise_for_status()
            data = resp.json()
            entries = [
                {"rank": e["rank"], "name": e["scientificName"], "taxid": str(e["taxonId"])}
                for e in data.get("lineage", [])
                if e.get("rank") in _BLAST_RANKS
            ]
            # append the organism itself
            entries.append({
                "rank": data.get("rank", "species"),
                "name": data.get("scientificName", ""),
                "taxid": tid,
            })
            _lineage.set(entries)
        except Exception as e:
            ui.notification_show(f"Lineage fetch failed: {e}", type="warning", duration=4)
            _lineage.set([])

    @render.ui
    def lineage_ui():
        """Display taxonomic lineage as a reference table (rank, name, taxid)."""
        entries = _lineage()
        if not entries:
            return None
        rows = [
            ui.tags.tr(
                ui.tags.td(e["rank"],   style="padding:1px 6px; color:#6c757d;"),
                ui.tags.td(ui.tags.em(e["name"]), style="padding:1px 6px;"),
                ui.tags.td(e["taxid"],  style="padding:1px 8px; font-family:monospace;"),
            )
            for e in entries
        ]
        table = ui.tags.table(
            ui.tags.tbody(*rows),
            style="font-size:0.75rem; border-collapse:collapse;",
        )
        return ui.tags.details(
            ui.tags.summary(
                "Taxids for BLAST scoping",
                style=(
                    "font-size:0.8rem; font-weight:600; cursor:pointer; "
                    "padding:0.35rem 0.5rem; color:#212529; user-select:none;"
                ),
            ),
            ui.div(
                ui.tags.p("Use these taxids to limit BLAST searches to specific taxa.",
                          style="font-size:0.78rem; color:#495057; margin:0 0 0.3rem;"),
                table,
                style="padding:0.3rem 0.5rem 0.4rem;",
            ),
            style=(
                "border:1px solid #adb5bd; border-radius:4px; "
                "margin-top:0.4rem; background:#f8f9fa;"
            ),
        )

    # ── genomic sequence auto-fetch (Ensembl) ─────────────────────────────────

    _genomic_fetch_status = reactive.Value("")

    @reactive.effect
    def _prefill_genomic_gene_symbol():
        """Prefill the auto-fetch gene-symbol box from the selected UniProt hit's gene name."""
        hit = _selected_hit()
        if hit and hit.get("gene"):
            ui.update_text("genomic_gene_symbol", value=hit["gene"])

    @reactive.effect
    @reactive.event(input.fetch_genomic_btn)
    def _fetch_genomic():
        """Fetch genomic FASTA from Ensembl and populate the paste textarea.

        Never blocks saving on failure — manual upload/paste remain usable.
        """
        taxid  = organism_taxid()
        symbol = input.genomic_gene_symbol().strip()
        if not taxid:
            ui.notification_show("Select an organism first.", type="warning", duration=5)
            return
        if not symbol:
            ui.notification_show("Enter a gene symbol first.", type="warning", duration=5)
            return
        try:
            flank_bp = int(input.genomic_flank_bp())
        except Exception:
            flank_bp = 2000

        try:
            fasta_text, meta = fetch_genomic_sequence(taxid, symbol, flank_bp=flank_bp)
        except Exception as e:
            _genomic_fetch_status.set("")
            ui.notification_show(f"Genomic fetch failed: {e}", type="error", duration=8)
            return

        ui.update_text_area("genomic_seq_paste", value=fasta_text)
        status = f"Fetched {meta['gene_id']} ({meta['species']}) {meta['region']}"
        _genomic_fetch_status.set(status)
        ui.notification_show(status, type="message", duration=6)

        if meta["ambiguous"]:
            ui.notification_show(
                f"'{symbol}' matched {len(meta['candidates'])} genes; chose "
                f"{meta['gene_id']} (longest genomic span, most likely the "
                f"complete gene model). Verify this is the gene you meant — "
                f"other matches: {', '.join(meta['candidates'])}.",
                type="warning", duration=10,
            )

    @render.ui
    def fetch_genomic_status():
        """Show the last successful auto-fetch's summary line."""
        status = _genomic_fetch_status()
        if not status:
            return None
        return ui.div(status, class_="section-hint")

    # ── PDB pLDDT validity check ──────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.input_file)
    def _check_input_file():
        """Warn if the uploaded PDB has no valid AlphaFold pLDDT B-factors."""
        files = input.input_file()
        if not files:
            return
        if not files[0]["name"].lower().endswith(".pdb"):
            _pdb_plddt_valid.set(True)
            return
        try:
            bfs = extract_bfactors_from_pdb(files[0]["datapath"])
            valid = bool(bfs) and all(0 <= v <= 100 for v in bfs) and (max(bfs) - min(bfs)) > 5
        except Exception:
            valid = False
        _pdb_plddt_valid.set(valid)
        if not valid:
            ui.notification_show(
                "The uploaded PDB does not appear to contain AlphaFold pLDDT scores "
                "(B-factors outside 0–100 or uniform). "
                "The AlphaFold Database will be searched for a matching structure.",
                type="warning",
                duration=10,
            )

    # ── UniProt protein search ────────────────────────────────────────────────

    def _parse_tsv_hits(tsv_text):
        """Parse a UniProt TSV response into a list of hit dicts.

        The TSV must include the 'Reviewed' column (requested via fields=...reviewed...)
        so hits carry a reviewed=True/False flag for downstream filtering.
        """
        lines = tsv_text.strip().splitlines()
        if len(lines) < 2:
            return []
        header = lines[0].split("\t")
        col = {h: i for i, h in enumerate(header)}
        hits = []
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) < len(header):
                continue
            acc      = parts[col.get("Entry", 0)].strip()
            gene     = parts[col.get("Gene Names (primary)", 1)].strip()
            protein  = parts[col.get("Protein names", 2)].strip()
            length   = parts[col.get("Length", 3)].strip()
            organism = parts[col.get("Organism", 4)].strip()
            taxid_h  = parts[col.get("Organism (ID)", 5)].strip()
            sequence = parts[col.get("Sequence", 6)].strip()
            afdb_col = parts[col.get("AlphaFoldDB", 7)].strip() if "AlphaFoldDB" in col else ""
            reviewed = parts[col.get("Reviewed", 8)].strip().lower() == "reviewed" if "Reviewed" in col else False
            hits.append({
                "accession":    acc,
                "gene":         gene,
                "protein":      protein,
                "length":       length,
                "organism":     organism,
                "taxid":        taxid_h,
                "has_afdb":     bool(afdb_col),
                "sequence":     sequence,
                "reviewed":     reviewed,
                "isoform":      None,
                "isoform_name": "",
                "safe_id":      acc.replace("-", "_"),
            })
        return hits

    @reactive.effect
    @reactive.event(input.uniprot_search_btn)
    def _uniprot_search():
        """Search UniProt for proteins by gene name or accession, scoped by organism if set.

        Strategy:
        1. Use gene:q as the primary search — this finds all entries whose primary gene
           name is q, including computationally defined isoforms stored as separate
           TrEMBL entries (e.g. the 14 unc-44 isoforms in C. elegans).
        2. If the gene search returns nothing (q is likely an accession, not a gene name),
           fall back to a plain full-text search.
        3. For each hit, also expand curated isoforms via ALTERNATIVE PRODUCTS comment
           (e.g. TP53's 9 manually curated isoforms all live inside P04637).
        """
        q = input.uniprot_query().strip() if "uniprot_query" in input else ""
        if not q:
            return

        tid = organism_taxid()
        scope = f" AND taxonomy_id:{tid}" if tid else ""
        fields = "accession,gene_primary,protein_name,length,organism_name,organism_id,reviewed,sequence,xref_alphafolddb"

        def _fetch_tsv(query):
            resp = requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={"query": query, "format": "tsv", "fields": fields, "size": 25},
                timeout=15,
            )
            resp.raise_for_status()
            return resp.text

        try:
            # gene-field search is the primary path; returns all entries with gene name = q,
            # including computationally defined isoforms stored as separate TrEMBL entries
            all_hits = _parse_tsv_hits(_fetch_tsv(f"gene:{q}{scope}"))

            if not all_hits:
                # q is likely an accession — fetch it directly rather than doing a
                # free-text search (which returns many unrelated hits)
                resp = requests.get(
                    f"https://rest.uniprot.org/uniprotkb/{q}",
                    params={"format": "json"},
                    timeout=15,
                )
                if resp.ok:
                    entry = resp.json()
                    org   = entry.get("organism", {})
                    prot  = entry.get("proteinDescription", {})
                    # extract primary gene name from nested structure
                    genes = entry.get("genes", [])
                    gene  = genes[0].get("geneName", {}).get("value", "") if genes else ""
                    # extract protein name
                    rec_name = prot.get("recommendedName", {})
                    prot_name = rec_name.get("fullName", {}).get("value", "") if rec_name else ""
                    seq   = entry.get("sequence", {}).get("value", "")
                    afdb  = any(
                        x.get("database") == "AlphaFoldDB"
                        for x in entry.get("uniProtKBCrossReferences", [])
                    )
                    rev   = entry.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"
                    all_hits = [{
                        "accession":    entry.get("primaryAccession", q),
                        "gene":         gene,
                        "protein":      prot_name,
                        "length":       str(entry.get("sequence", {}).get("length", len(seq))),
                        "organism":     org.get("scientificName", ""),
                        "taxid":        str(org.get("taxonId", "")),
                        "has_afdb":     afdb,
                        "sequence":     seq,
                        "reviewed":     rev,
                        "isoform":      None,
                        "isoform_name": "",
                        "safe_id":      q.replace("-", "_"),
                    }]
        except Exception as e:
            ui.notification_show(f"UniProt search failed: {e}", type="error", duration=5)
            return

        if not all_hits:
            _uniprot_hits.set([])
            _selected_hit.set(None)
            _selected_pdb.set("")
            scope_msg = f" (scoped to taxonomy_id {tid})" if tid else ""
            ui.notification_show(f"No UniProt hits found{scope_msg}.", type="warning", duration=4)
            return

        # prefer reviewed (SwissProt) entries when available — they have curated isoform
        # annotations via ALTERNATIVE PRODUCTS; only include unreviewed (TrEMBL) entries
        # when no reviewed entry exists (e.g. unc-44 which is entirely unreviewed).
        # Exception: always keep TrEMBL entries whose primary gene name exactly matches
        # the query — otherwise e.g. searching "rab-5" drops the TrEMBL P91857 entry
        # because a reviewed RAB family member also appears in the results.
        reviewed_hits = [h for h in all_hits if h["reviewed"]]
        if reviewed_hits:
            exact_unreviewed = [h for h in all_hits
                                if not h["reviewed"] and h["gene"].lower() == q.lower()]
            base_hits = reviewed_hits + exact_unreviewed
        else:
            base_hits = all_hits

        # exact gene name matches first, then reviewed, then the rest
        def _hit_sort_key(h):
            exact = h["gene"].lower() != q.lower()   # False (0) sorts before True (1)
            return (exact, not h["reviewed"])

        base_hits.sort(key=_hit_sort_key)

        # expand curated isoforms (entries with ALTERNATIVE PRODUCTS comment, e.g. TP53)
        seen = {}
        for h in base_hits:
            for exp_h in _fetch_isoforms(h["accession"], h):
                if exp_h["accession"] not in seen:
                    seen[exp_h["accession"]] = exp_h

        _selected_hit.set(None)
        _selected_pdb.set("")
        _uniprot_hits.set(list(seen.values()))
        for h in seen.values():
            _register_hit_select(h)

    def _register_hit_select(hit):
        """Register a reactive effect that selects a hit when its Use button fires."""
        sid = hit["safe_id"]
        acc = hit["accession"]   # capture in closure

        @reactive.effect
        @reactive.event(lambda: input[f"use_hit_{sid}"]())
        def _handler():
            _selected_hit.set(hit)

            if not hit["has_afdb"]:
                _selected_pdb.set("")
                ui.notification_show(
                    f"{acc}: no AFDB structure — will be searched at run time.",
                    type="warning", duration=6,
                )
                return

            # AFDB indexes by canonical accession — strip isoform suffix (P30628-1 → P30628)
            parts = acc.rsplit("-", 1)
            afdb_acc = parts[0] if len(parts) == 2 and parts[1].isdigit() else acc
            try:
                pdb_bytes = ebi_rest.dbfetch("afdb", afdb_acc, "pdb", "raw")
                pdb_text = pdb_bytes.decode("utf-8", errors="replace")
            except Exception as e:
                ui.notification_show(f"AFDB fetch failed for {afdb_acc}: {e}", type="error", duration=5)
                _selected_pdb.set("")
                return

            if "ERROR" in pdb_text[:50]:
                _selected_pdb.set("")
                ui.notification_show(
                    f"{afdb_acc}: AFDB structure not available — will be searched at run time.",
                    type="warning", duration=6,
                )
                return

            _selected_pdb.set(pdb_text)
            ui.notification_show(
                f"Selected {acc} ({hit['gene'] or acc}) — AFDB structure loaded.",
                type="message", duration=4,
            )

    @render.ui
    def uniprot_results_ui():
        """Render UniProt hits as an accordion; each panel expands to show FASTA + select."""
        all_hits = _uniprot_hits()
        if not all_hits:
            return None

        # apply display filters
        afdb_only = input.uniprot_afdb_only()
        hits = [h for h in all_hits if not afdb_only or h["has_afdb"]]

        hidden = len(all_hits) - len(hits)
        if not hits:
            return ui.div(
                ui.tags.p(
                    f"No results match the current filters "
                    f"({hidden} result{'s' if hidden != 1 else ''} hidden).",
                    style="font-size:0.78rem; color:#6c757d; margin:0.3rem 0;",
                ),
            )

        panels = []
        for hit in hits:
            acc      = hit["accession"]
            gene     = hit["gene"] or acc
            iso_name = hit.get("isoform_name", "")
            iso_num  = hit.get("isoform", None)
            sid      = hit["safe_id"]

            afdb_badge = (
                ui.tags.span("✓ AFDB", class_="uc-afdb-yes")
                if hit["has_afdb"]
                else ui.tags.span("no AFDB", class_="uc-afdb-no")
            )
            # show isoform number for non-canonical isoforms; always show length
            iso_badge = (
                ui.tags.span(f"isoform {iso_num}", class_="uc-isoform-badge")
                if iso_num and iso_num != "1"
                else None
            )
            len_badge = ui.tags.span(f"{hit['length']} aa", class_="uc-canonical-badge")

            # show isoform common name (e.g. "p53alpha") alongside gene name
            title = f"{gene} · {iso_name}" if iso_name else gene
            header = ui.span(
                title,
                ui.tags.span(acc, class_="task-type-badge"),
                *([iso_badge] if iso_badge else []),
                len_badge,
                afdb_badge,
            )

            body = ui.div(
                # protein name + organism line
                ui.div(
                    ui.tags.span(hit["protein"],
                                 style="font-size:0.78rem; color:#495057;"),
                    ui.tags.span(f" · {hit['organism']}",
                                 style="font-size:0.75rem; color:#6c757d; font-style:italic;"),
                    style="margin-bottom:0.4rem;",
                ),
                # FASTA sequence in a scrollable code block
                ui.tags.pre(
                    _format_fasta(hit),
                    style=(
                        "font-size:0.7rem; background:#f8f9fa; border:1px solid #adb5bd; "
                        "border-radius:3px; padding:0.4rem 0.6rem; overflow-x:auto; "
                        "white-space:pre; margin-bottom:0.4rem; "
                        "max-height:120px; overflow-y:auto;"
                    ),
                ),
                ui.div(
                    ui.input_action_button(
                        f"use_hit_{sid}", "✓ Use this sequence",
                        class_="btn-sm btn-primary",
                    ),
                    style="margin-top:0.3rem;",
                ),
            )

            panels.append(ui.accordion_panel(header, body, value=sid))

        n_shown  = len(hits)
        hint = f"{n_shown} result{'s' if n_shown != 1 else ''}"
        if hidden:
            hint += f" ({hidden} hidden by filters)"
        return ui.div(
            ui.tags.p(
                f"{hint} — open a card to inspect and select",
                style="font-size:0.78rem; color:#6c757d; margin:0.3rem 0 0.2rem;",
            ),
            ui.accordion(*panels, open=None, multiple=True, class_="task-accordion"),
        )

    @render.text
    def seq_source_status():
        """One-line indicator showing which protein source will be used at save time."""
        hit = _selected_hit()
        if hit is not None:
            pdb_note = " + AFDB structure" if _selected_pdb() else ""
            return f"Using UniProt {hit['accession']} ({hit['gene'] or hit['organism']}){pdb_note}"
        if input.input_file():
            return "Using uploaded file"
        if input.protein_seq_paste().strip():
            return "Using pasted sequence"
        return "No sequence selected yet."

    # ── preset refresh (called after saving a new preset) ────────────────────

    def _populate_presets():
        """Refresh the Load preset dropdown after a new preset is saved."""
        names = [p.stem for p in sorted(_PARAMS_DIR.glob("*.json"))]
        choices = {"": "— select —", **{n: n for n in names}}
        ui.update_select("load_preset", choices=choices, selected="")

    # ── button enable / disable ───────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.task_label)
    def _toggle_add():
        """Enable Add only when a task label has been typed."""
        ui.update_action_button("add_task", disabled=not bool(input.task_label()))

    @reactive.effect
    def _toggle_load():
        """Enable Load only when a preset is selected."""
        ui.update_action_button("load_preset_btn",
                                disabled=not bool(input.load_preset()))

    @reactive.effect
    def _toggle_save_preset():
        """Enable Save preset only when a name is typed and tasks exist."""
        ui.update_action_button(
            "save_preset_btn",
            disabled=not (bool(input.preset_name()) and len(tasks()) > 0),
        )

    @reactive.calc
    def _ready():
        """True when all required fields for Save Analysis are present."""
        has_sequence = (
            bool(input.input_file())
            or _selected_hit() is not None
            or bool(input.protein_seq_paste().strip())
        )
        return (
            bool(input.email())
            and has_sequence
            and bool(input.working_dir())
            and len(tasks()) > 0
        )

    @reactive.effect
    @reactive.event(_ready)
    def _toggle_save():
        ui.update_action_button("save_analysis", disabled=not _ready())

    # ── auto-fill working directory ───────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.run_name)
    def _autofill_dir():
        """Derive working_dir from run_name when the user types a name."""
        if not input.run_name():
            return
        base = GLOBAL_DEFAULTS.get("working_dir") or "data"
        ui.update_text("working_dir", value=f"{base}/{input.run_name()}/")

    # ── save-status hint ──────────────────────────────────────────────────────

    @render.text
    def save_status():
        """One-line hint next to the Save button."""
        if not _ready():
            return "Requires email, a sequence file or UniProt hit, an analysis name, and at least one task."
        if input.save_analysis() == 0:
            return "Ready — click Save Analysis to write the run configuration."
        return f"Saved → {input.working_dir()}{input.run_name()}.run.json"

    @render.text
    def task_type_desc():
        """One-line description of the currently selected analysis type."""
        return TASK_DESCRIPTIONS.get(input.task_type(), "")

    # ── task state helpers ────────────────────────────────────────────────────

    def _snapshot():
        """Read all current param inputs into task_snap (copy-and-set for reactivity)."""
        snap = {}
        for t in tasks():
            hidden = task_hidden(t["name"])
            snap[t["id"]] = {}
            for p in t["params"]:
                if p in hidden:
                    continue
                wid = f"{t['id']}_{p}"
                if wid in input:
                    snap[t["id"]][p] = input[wid]()
        task_snap.set(snap)

    def _collect_args(task):
        """Read input values for all params in a task, falling back to registry defaults."""
        args = {}
        for p, pcfg in TASK_PARAMETERS[task["name"]]["params"].items():
            wid = f"{task['id']}_{p}"
            # is_set() is False for hidden params (no widget) — fall back to default
            if input[wid].is_set():
                args[p] = input[wid]()
            else:
                args[p] = pcfg["default"]
        return args

    # ── task CRUD ─────────────────────────────────────────────────────────────

    def _register_remove(task_id):
        """Register a one-off reactive effect that removes a task when its button fires."""
        @reactive.effect
        @reactive.event(lambda: input[f"{task_id}_remove"]())
        def _handler():
            _snapshot()
            tasks.set([t for t in tasks() if t["id"] != task_id])

    def _register_table_upload(task_id):
        """Copy an uploaded TSV to /tables/, update snap, and re-render the card."""
        @reactive.effect
        @reactive.event(lambda: input[f"{task_id}_add_table"]())
        def _handler():
            files = input[f"{task_id}_add_table"]()
            if not files:
                return
            f = files[0]
            dest = _TABLES_DIR / f["name"]
            shutil.copy(f["datapath"], dest)
            snap = dict(task_snap())
            snap.setdefault(task_id, {})["scores_file"] = str(dest)
            task_snap.set(snap)
            tasks.set(list(tasks()))
            ui.notification_show(f"Table '{f['name']}' added.", type="message", duration=3)

    @reactive.effect
    @reactive.event(input.add_task)
    def _add_task():
        """Append a new task card for the selected analysis type."""
        ttype = input.task_type()
        label = input.task_label().strip()
        if not ttype or not label:
            return
        _snapshot()
        hidden = task_hidden(ttype)
        params  = {p: v for p, v in task_defaults(ttype).items() if p not in hidden}
        tips    = task_tooltips(ttype)
        choices = task_choices(ttype)
        task = make_task(ttype, label, params, tips, choices)
        task["start_open"] = True   # newly added tasks open by default
        tasks.set(tasks() + [task])
        _register_remove(task["id"])
        if ttype == "scores":
            _register_table_upload(task["id"])
        ui.update_text("task_label", value="")

    # ── preset load / save ────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.load_preset_btn)
    def _load_preset():
        """Load a saved task set from params/ and append to the current list."""
        name = input.load_preset()
        if not name:
            return
        path = _PARAMS_DIR / f"{name}.json"
        try:
            data = json.load(open(path))
        except FileNotFoundError:
            ui.notification_show(f"Preset '{name}' not found.", type="error")
            return
        _snapshot()
        new_tasks = []
        for tid, entry in data.items():
            ttype  = entry.get("type") or entry.get("analysis", "")  # accept old "analysis" key
            hidden = task_hidden(ttype)
            # use saved args as params, supplement tooltips + choices from registry
            params  = {p: v for p, v in entry["args"].items() if p not in hidden and p != "output"}
            tips    = task_tooltips(ttype)
            choices = task_choices(ttype)
            task = {"id": tid, "name": ttype, "params": params, "tooltips": tips,
                    "choices": choices, "start_open": False}   # preloaded tasks start collapsed
            new_tasks.append(task)
            _register_remove(tid)
            if ttype == "scores":
                _register_table_upload(tid)
        tasks.set(tasks() + new_tasks)

    @reactive.effect
    @reactive.event(input.save_preset_btn)
    def _save_preset():
        """Write the current task set (no global block) to a new params/ file."""
        _snapshot()
        name = input.preset_name().strip()
        if not name:
            return
        result = {}
        wd, rn = input.working_dir(), input.run_name()
        for task in tasks():
            args = _collect_args(task)
            result.update(build_defaults_entry(task["id"], task["name"], args,
                                               task_output_suffix(task["name"]), wd, rn))
        path = _PARAMS_DIR / f"{name}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(result, f, indent=4)
        ui.notification_show(f"Saved preset '{name}'.", type="message", duration=3)
        _populate_presets()

    # ── task card rendering ───────────────────────────────────────────────────

    @render.ui
    def task_cards():
        """Render tasks as a collapsible accordion; new tasks open, preloaded closed."""
        current_tasks = tasks()
        snap = task_snap()

        if not current_tasks:
            return ui.p(
                "No tasks yet — choose a type and label above, then click Add.",
                style="color:#6c757d; margin-top:0.4rem; font-size:0.8rem;",
            )

        open_ids = [t["id"] for t in current_tasks if t.get("start_open", False)]

        panels = []
        for task in current_tasks:
            label = task["id"].rsplit("_", 1)[0]
            widgets = _build_param_inputs(task, snap)
            add_table = (
                ui.div(
                    compact_file_input(f"{task['id']}_add_table", "＋ Add table (.tsv)",
                        accept=[".tsv"]),
                    style="margin-top:0.2rem;",
                )
                if task["name"] == "scores" else None
            )
            panels.append(
                ui.accordion_panel(
                    ui.span(
                        label,
                        ui.tags.span(task["name"], class_="task-type-badge"),
                    ),
                    # params grid
                    ui.layout_column_wrap(*widgets, width=1/2) if widgets else
                    ui.p("No configurable parameters.", style="color:#6c757d; font-size:0.8rem;"),
                    # scores-only: upload a new property table
                    *([add_table] if add_table else []),
                    # remove button at bottom of body
                    ui.div(
                        ui.input_action_button(
                            f"{task['id']}_remove", "✕ Remove task",
                            class_="btn-sm btn-outline-danger",
                        ),
                        style="margin-top:0.6rem;",
                    ),
                    value=task["id"],
                )
            )

        return ui.accordion(
            *panels,
            open=open_ids,
            multiple=True,
            class_="task-accordion",
        )

    # ── load run JSON (direct upload in setup tab) ────────────────────────────

    @reactive.effect
    @reactive.event(input.upload_run_json)
    def _handle_json_upload():
        """Store an uploaded run JSON path in shared state (same as other tabs)."""
        info = input.upload_run_json()
        if info:
            shared_json.set(info[0]["datapath"])

    # ── populate setup from any loaded JSON ───────────────────────────────────

    def _tasks_from_run_json(data):
        """Parse a run-JSON dict into a list of setup task dicts, skipping reagents."""
        new_tasks = []
        for tid, entry in data.get("tasks", {}).items():
            ttype = entry.get("type", "")
            if not ttype or ttype == "reagents":
                continue
            hidden = task_hidden(ttype)
            raw_args = entry.get("args", {})
            # strip global keys and internal-only fields
            params = {
                p: v for p, v in raw_args.items()
                if p not in hidden and p not in GLOBAL_KEYS and p != "output"
            }
            tips    = task_tooltips(ttype)
            choices = task_choices(ttype)
            task = {"id": tid, "name": ttype, "params": params,
                    "tooltips": tips, "choices": choices, "start_open": False}
            new_tasks.append(task)
        return new_tasks

    @reactive.effect
    @reactive.event(shared_json)
    def _load_from_shared_json():
        """Populate setup task cards (and global fields) whenever any tab loads a JSON."""
        if _skip_next_load[0]:
            _skip_next_load[0] = False
            return
        path = shared_json.get()
        if not path:
            return
        try:
            with open(path) as f:
                data = json.load(f)
        except Exception:
            return
        if "tasks" not in data:
            return   # preset format, not a run JSON — ignore

        gb = data.get("global", {})

        # populate global text fields when present in the JSON
        if gb.get("email"):
            ui.update_text("email", value=gb["email"])
        if gb.get("run_name"):
            ui.update_text("run_name", value=gb["run_name"])
        if gb.get("working_dir"):
            ui.update_text("working_dir", value=gb["working_dir"])

        new_tasks = _tasks_from_run_json(data)
        _snapshot()
        for t in new_tasks:
            _register_remove(t["id"])
            if t["name"] == "scores":
                _register_table_upload(t["id"])
        tasks.set(new_tasks)
        if new_tasks:
            ui.notification_show(
                f"Loaded {len(new_tasks)} task{'s' if len(new_tasks) != 1 else ''} from JSON.",
                type="message", duration=3,
            )

    # ── save analysis ─────────────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.save_analysis)
    async def _save_analysis():
        """Write the run JSON and publish its path via shared_json.

        Protein-source precedence: UniProt hit > uploaded file > pasted sequence.
        Genomic-source precedence: uploaded file > pasted sequence.
        All branches produce (input_file, pdb_path, pdb_is_valid_plddt).
        """
        wd    = input.working_dir()
        rn    = input.run_name()
        email = input.email()

        Path(wd).mkdir(parents=True, exist_ok=True)

        # ── protein source resolution ─────────────────────────────────────────
        hit = _selected_hit()
        if hit is not None:
            # UniProt branch: write sequence from search result; write AFDB PDB if fetched
            fasta_dest = wd + rn + ".fa"
            save_fasta(hit["accession"], hit["sequence"], fasta_dest)
            input_file = fasta_dest

            pdb_path = ""
            pdb_text = _selected_pdb()
            if pdb_text:
                pdb_dest = wd + rn + ".pdb"
                with open(pdb_dest, "w") as f:
                    f.write(pdb_text)
                pdb_path = pdb_dest
            # AFDB structures always carry valid pLDDT B-factors
            pdb_is_valid_plddt = bool(pdb_path)
        elif input.input_file():
            # upload branch: copy protein file to working dir (existing logic)
            tmp_protein = Path(input.input_file()[0]["datapath"])
            ext = ".pdb" if tmp_protein.suffix.lower() == ".pdb" else ".fa"
            protein_dest = wd + rn + ext
            shutil.copy(tmp_protein, protein_dest)

            pdb_path = ""
            input_file = protein_dest
            if ext == ".pdb":
                pdb_path = protein_dest
                fasta_dest = protein_dest.replace(".pdb", ".fa")
                save_fasta(Path(pdb_path).stem, get_sequence(pdb_path), fasta_dest)
                input_file = fasta_dest
            pdb_is_valid_plddt = _pdb_plddt_valid()
        else:
            # pasted sequence branch
            raw = input.protein_seq_paste().strip()
            fasta_dest = wd + rn + ".fa"
            seq = raw if raw.startswith(">") else f">{rn}\n{raw}"
            with open(fasta_dest, "w") as f:
                f.write(seq + "\n")
            input_file = fasta_dest
            pdb_path = ""
            pdb_is_valid_plddt = False

        # ── copy optional genomic FASTA ───────────────────────────────────────
        genomic_path = ""
        if input.input_genomic():
            tmp_gen = Path(input.input_genomic()[0]["datapath"])
            genomic_dest = wd + rn + ".genomic" + tmp_gen.suffix
            shutil.copy(tmp_gen, genomic_dest)
            genomic_path = genomic_dest
        elif input.genomic_seq_paste().strip():
            raw_gen = input.genomic_seq_paste().strip()
            genomic_dest = wd + rn + ".genomic.fa"
            seq_gen = raw_gen if raw_gen.startswith(">") else f">{rn}_genomic\n{raw_gen}"
            with open(genomic_dest, "w") as f:
                f.write(seq_gen + "\n")
            genomic_path = genomic_dest

        # ── build the global block ────────────────────────────────────────────
        global_block = build_global_block(
            email        = email,
            run_name     = rn,
            working_dir  = wd,
            input_file   = input_file,
            pdb          = pdb_path,
            genomic_file = genomic_path,
        )

        # ── build per-task entries ────────────────────────────────────────────
        task_entries = []
        for task in tasks():
            args = _collect_args(task)
            if task["name"] == "plddt":
                # use the supplied PDB only when it contains valid pLDDT B-factors
                use_supplied_pdb = bool(pdb_path) and pdb_is_valid_plddt
                args["existing_AF2"] = 0 if use_supplied_pdb else 1
                # write pdb into task args so parse_run merge doesn't lose it
                args["pdb"] = pdb_path if use_supplied_pdb else ""
                # always inject the resolved taxid so AFDB lookup can filter by organism
                args["taxid"] = organism_taxid()
            task_entries.append(
                build_task_entry(task["id"], task["name"], args,
                                 task_output_suffix(task["name"]), wd, rn)
            )

        # ── build reagents entry or warn ──────────────────────────────────────
        reagents_entry = None
        if genomic_path:
            reagents_entry = build_reagents_entry(
                defaults     = task_defaults("reagents"),
                genomic_path = genomic_path,
                out_suffix   = task_output_suffix("reagents"),
                working_dir  = wd,
                run_name     = rn,
            )
        else:
            ui.notification_show(
                "No genomic FASTA uploaded — CRISPR reagent design will not run. "
                "Upload a genomic region FASTA and re-save to enable it.",
                type="warning",
                duration=8,
            )

        run_json = build_run_json(
            global_block   = global_block,
            task_entries   = task_entries,
            reagents_entry = reagents_entry,
        )

        out_path = f"{wd}{rn}.run.json"
        with open(out_path, "w") as f:
            json.dump(run_json, f, indent=4)

        _skip_next_load[0] = True
        shared_json.set(out_path)

        # give the save a moment to register before jumping to the Progress tab
        await asyncio.sleep(0.5)
        ui.update_navs("main_tabs", selected="progress", session=session.root_scope())
