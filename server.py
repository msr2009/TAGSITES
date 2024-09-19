from modules.setup_server import setup_server
from modules.progress_server import progress_server
from modules.results_server import results_server
from modules.reagents_server import reagents_server

from shiny import reactive

shared_values = reactive.Value("")

def app_server(input, output, session):
	#call each server
	setup_server(input, output, session, shared_values)
	progress_server(input, output, session, shared_values)
	results_server(input, output, session, shared_values)
	reagents_server(input, output, session, shared_values)

