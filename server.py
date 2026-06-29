from shiny import reactive

from modules.setup_server import setup_server
from modules.progress_server import progress_server
from modules.results_server import results_server
from modules.reagents_server import reagents_server


def app_server(input, output, session):
    # per-session shared state: the path to the active run JSON
    # created inside app_server so each browser session gets its own value
    shared_values = reactive.Value("")

    # chosen tag sites — list of residue positions (ints) from the Results page
    shared_sites = reactive.Value([])

    # incremented by progress_server each time a run completes so results_server
    # reloads even when the JSON path hasn't changed
    shared_results_trigger = reactive.Value(0)

    setup_server("setup", shared_json=shared_values)
    progress_server("progress", shared_json=shared_values,
                    shared_results_trigger=shared_results_trigger)
    results_server("results", shared_json=shared_values, shared_sites=shared_sites,
                   shared_results_trigger=shared_results_trigger)
    reagents_server("reagents", shared_json=shared_values, shared_sites=shared_sites)
