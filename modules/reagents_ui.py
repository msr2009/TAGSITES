from shiny import ui, module


@module.ui
def reagents_ui():
    return ui.page_fluid(
        ui.h2("Design Reagents"),
        ui.output_ui("show_or_upload_json_reagents"),
        ui.hr(),
    )
