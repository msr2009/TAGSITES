from shiny import ui
from modules.setup_ui import setup_ui
from modules.progress_ui import progress_ui
from modules.results_ui import results_ui
from modules.reagents_ui import reagents_ui

app_ui = ui.page_fluid(
    ui.tags.style(
        "#download_app{"
        "color:#b0bac4;border-color:#b0bac4;"
        "position:relative;top:-3px;"
        "}"
        "#download_app:hover{color:#fff;background:#b0bac4;border-color:#b0bac4;}"
    ),
    ui.navset_card_tab(
        ui.nav_panel("Analysis Setup", setup_ui("setup")),
        ui.nav_panel("Progress", progress_ui("progress")),
        ui.nav_panel("Results", results_ui("results")),
        ui.nav_panel("Design Reagents", reagents_ui("reagents")),
        ui.nav_control(
            ui.download_button("download_app", "⬇ Download App",
                class_="btn-sm btn-outline-secondary",
                style="margin: 4px 0;"),
        ),
    )
)
