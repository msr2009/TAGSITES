from shiny import ui

from modules.batch_nav_ui import batch_guides_ui, batch_nav_ui
from modules.batch_progress_ui import batch_progress_ui
from modules.batch_setup_ui import batch_setup_ui
from modules.results_ui import results_ui

# same tab colours as the single-protein app (ui.py)
_TAB_CSS = (
    "ul#main_tabs.nav-tabs{background:#fff;border-bottom-color:#dee2e6;}"
    "ul#main_tabs.nav-tabs > li.nav-item > a.nav-link{"
    "color:#212529;font-weight:500;background:#9ec5fe;border-color:#0a58ca;}"
    "ul#main_tabs.nav-tabs > li.nav-item > a.nav-link:hover{"
    "background:#0a58ca;border-color:#0a58ca;}"
    "ul#main_tabs.nav-tabs > li.nav-item > a.nav-link.active{"
    "color:#fff;font-weight:700;background:#0a58ca;border-color:#0a58ca;}"
)

app_ui = ui.page_fluid(
    ui.tags.style(_TAB_CSS),
    ui.tags.title("TAGSITES batch"),
    # tab id main_tabs: the shared modules switch tabs by this id (ui.update_navs / update_navset)
    ui.navset_card_tab(
        ui.nav_panel("(1) Setup Batch", batch_setup_ui("setup"), value="setup"),
        ui.nav_panel("(2) Prepare", batch_progress_ui("progress"), value="progress"),
        # the guide picker belongs to the nav module ("nav" namespace) but sits in Results
        ui.nav_panel("(3) Review Sites", batch_nav_ui("nav"),
                     results_ui("results", between=batch_guides_ui("nav")), value="results"),
        id="main_tabs",
    ),
)
