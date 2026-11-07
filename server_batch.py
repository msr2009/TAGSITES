import os

from shiny import reactive

from modules.batch_nav_server import batch_nav_server
from modules.batch_progress_server import batch_progress_server
from modules.batch_setup_server import batch_setup_server
from modules.results_server import results_server


def app_server(input, output, session):
    """Wire the batch app: Setup -> Prepare -> Review (the single-protein Results module)."""
    # the open batch: plain dict (worker threads update it in place) plus two signals —
    # rev re-renders views of it, start (re)runs Prepare
    bstate = {"path": None, "data": None, "rev": reactive.Value(0), "start": reactive.Value(0)}

    # the Results module's inputs, exactly as in the single-protein app
    shared_json = reactive.Value("")
    shared_sites = reactive.Value([])
    shared_results_trigger = reactive.Value(0)
    # {"path": run json, "ranks": {position: rank}} for the loaded protein's suggested
    # regions, so batch_nav can record each accepted site's rank, and only its own
    # (see utils.scoring.pick_suggested_regions)
    shared_regions = reactive.Value({})
    # set by results_server once a protein's run is loaded and its sites restored
    loaded_path = reactive.Value("")
    # the open batch's Setup site filters (defaults of the Rescoring switches)
    site_filters = reactive.Value({})

    @reactive.effect
    @reactive.event(bstate["start"])
    def _batch_filters():
        """Take the site filters from the batch Setup just handed over."""
        data = bstate["data"] or {}
        site_filters.set(dict(data.get("settings", {}).get("site_filters") or {}))

    def _kept_sites(run_json):
        """Accepted sites of the batch protein at run_json; None (recompute) for any other."""
        real = os.path.realpath(run_json)
        for p in (bstate["data"] or {}).get("proteins", []):
            if p.get("run_json") and os.path.realpath(p["run_json"]) == real:
                return list(p.get("sites") or []) if p.get("status") == "accepted" else None
        return None

    batch_setup_server("setup", bstate=bstate)
    batch_progress_server("progress", bstate=bstate)
    results_server("results", shared_json=shared_json, shared_sites=shared_sites,
                   shared_results_trigger=shared_results_trigger, auto_suggest=True,
                   loaded_path=loaded_path, site_filters=site_filters,
                   kept_sites=_kept_sites, shared_regions=shared_regions)
    batch_nav_server("nav", bstate=bstate, shared_json=shared_json, shared_sites=shared_sites,
                     loaded_path=loaded_path, shared_regions=shared_regions)
