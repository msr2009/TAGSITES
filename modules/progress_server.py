from shiny import reactive

def progress_server(input, output, session, shared_values):
	
	#enable run_analysis button

	#load json into dict

	#update load json file to use shared_value

	@reactive.Effect
	@reactive.event(input.input1)
	def _():
		# Example server logic for Page 1
		output.output1.set(f"You entered: {input.input1()}")

