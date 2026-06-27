from pathlib import Path
from shiny import App
from ui import app_ui
from server import app_server

app = App(app_ui, app_server, static_assets=Path(__file__).parent / "www")

if __name__ == "__main__":
    app.run()
