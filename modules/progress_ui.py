from shiny import ui, module


@module.ui
def progress_ui():
    return ui.page_fluid(
        ui.h3("Load Analysis"),
        ui.output_ui("show_or_upload_json"),
        ui.hr(),
        ui.row(
            ui.column(4,
                ui.input_action_button("run_analysis", "Run Analysis", disabled=True,
                                       class_="btn-primary"),
            ),
            ui.column(4,
                ui.download_button("download_results", "Download Results",
                                   class_="btn-secondary"),
            ),
        ),
        ui.hr(),
        ui.h4("Task Status"),
        ui.output_ui("task_cards"),
    )
