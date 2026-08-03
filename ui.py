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
        "ul#main_tabs.nav-tabs{background:#fff;border-bottom-color:#dee2e6;}"
        # inactive tabs: $blue-200; active tab: $blue-600
        "ul#main_tabs.nav-tabs > li.nav-item > a.nav-link{"
        "color:#212529;font-weight:500;background:#9ec5fe;"
        "border-color:#0a58ca;"
        "}"
        "ul#main_tabs.nav-tabs > li.nav-item > a.nav-link:hover{"
        "background:#0a58ca;border-color:#0a58ca;"
        "}"
        "ul#main_tabs.nav-tabs > li.nav-item > a.nav-link.active{"
        "color:#fff;font-weight:700;background:#0a58ca;"
        "border-color:#0a58ca;"
        "}"
    ),
    ui.navset_card_tab(
        ui.nav_panel("(1) Setup Analyses", setup_ui("setup"), value="setup"),
        ui.nav_panel("(2) Run Analyses", progress_ui("progress"), value="progress"),
        ui.nav_panel("(3) View Results", results_ui("results"), value="results"),
        ui.nav_panel("(4) Create Reagents", reagents_ui("reagents"), value="reagents"),
        ui.nav_control(
            ui.download_button("download_app", "⬇ Download App",
                class_="btn-sm btn-outline-secondary",
                style="margin: 4px 0;"),
        ),
        id="main_tabs",
    )
)
