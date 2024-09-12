from shiny import reactive

def reagents_server(input, output, session, shared_values):
	@reactive.Effect
	@reactive.event(input.input1)
	def _():
		# Example server logic for Page 1
		output.output1.set(f"You entered: {input.input1()}")
