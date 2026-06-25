from shiny import ui, module

from config import INPUT_JSON, SELECTABLE_TASKS


@module.ui
def setup_ui():
    return ui.page_fluid(

        ui.tags.style("""
            .task-card { margin-bottom: 1rem; }
            .task-card .card-header {
                display: flex;
                align-items: center;
                justify-content: space-between;
                font-weight: 600;
            }
            .tip-icon {
                cursor: help;
                color: #6c757d;
                font-size: 0.85em;
                margin-left: 4px;
                vertical-align: middle;
            }
            .sidebar-section { margin-bottom: 1rem; }
        """),

        ui.layout_sidebar(

            # ── Sidebar: global config ───────────────────────────────────────
            ui.sidebar(
                ui.h5("Run Configuration"),
                ui.output_ui("sidebar_global_inputs"),
                ui.hr(),
                ui.h5("Defaults"),
                ui.input_selectize("load_default_tasks", "Load saved defaults",
                                   choices=[], multiple=False, selected=None,
                                   options={"create": False}),
                ui.input_action_button("load_defaults_button", "Load",
                                       disabled=True, class_="btn-sm btn-outline-secondary w-100"),
                ui.hr(),
                ui.input_text("new_default_name", "Save current tasks as…", ""),
                ui.input_action_button("save_default_tasks_button", "Save as Default",
                                       disabled=True, class_="btn-sm btn-outline-secondary w-100"),
                width=320,
            ),

            # ── Main panel: task builder + task cards ────────────────────────
            ui.h4("Add Analysis Tasks"),

            ui.card(
                ui.card_body(
                    ui.layout_column_wrap(
                        ui.input_select("task_selector", "Analysis type", SELECTABLE_TASKS),
                        ui.input_text("task_desc_name", "Task label (required)"),
                        width=1/2,
                    ),
                    ui.input_action_button("add_task", "＋ Add Task",
                                           disabled=True, class_="btn-primary"),
                ),
                class_="mb-3",
            ),

            ui.output_ui("task_params_container"),

            ui.hr(),

            ui.layout_column_wrap(
                ui.input_action_button("save_analysis", "Save Analysis",
                                       disabled=True, class_="btn-success"),
                ui.output_text_verbatim("tip_save_analysis"),
                width=1/2,
            ),
        )
    )
