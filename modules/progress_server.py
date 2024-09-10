from shiny import reactive

def progress_server(input, output, session):
	@reactive.Effect
	@reactive.event(input.input1)
	def _():
		# Example server logic for Page 1
		output.output1.set(f"You entered: {input.input1()}")

