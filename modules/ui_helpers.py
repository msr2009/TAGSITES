"""ui_helpers.py — shared UI helper functions and CSS snippets."""

from shiny import ui


# CSS for compact_file_input — include in any tab that uses it
COMPACT_FILE_CSS = """
    /* compact-file-input row */
    .compact-file-input {
        display: flex;
        align-items: center;
        gap: 0.5rem;
        flex-wrap: nowrap;
        margin-bottom: 0.25rem;
    }
    .compact-file-input .cfi-label {
        white-space: nowrap;
        flex-shrink: 0;
        font-size: 0.8rem;
        font-weight: 500;
        color: #343a40;
        margin: 0;
    }
    .compact-file-input .cfi-filename {
        font-size: 0.75rem;
        color: #6c757d;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
        max-width: 220px;
    }
    /* Shiny input_file hidden off-screen — functional but invisible */
    .cfi-hidden {
        position: fixed;
        top: -9999px;
        left: -9999px;
    }
    .cfi-hidden .shiny-input-container { margin: 0 !important; }
"""


def compact_file_input(input_id, label, **kwargs):
    """Custom file input: btn-sm btn-outline-secondary button + filename span.

    The actual Shiny input_file is hidden off-screen (stays reactive); a styled
    <button> triggers it via inp.click() and updates the filename display.
    """
    wrapper_id = f"cfi-{input_id}"

    # label may be a plain string or a Shiny Tag (e.g. label_with_tip returns ui.span)
    label_el = ui.tags.span(label, class_="form-label cfi-label") if label is not None else None

    return ui.div(
        *([label_el] if label_el is not None else []),
        ui.tags.button(
            "Choose file",
            type="button",
            class_="btn btn-sm btn-outline-secondary cfi-trigger",
        ),
        ui.tags.span("No file selected", class_="cfi-filename"),
        # Shiny's actual input — off-screen, invisible, but fully reactive
        ui.div(ui.input_file(input_id, None, **kwargs), class_="cfi-hidden"),
        # wire the button → hidden input; update filename on change
        ui.tags.script(
            f"""(function(){{
  function init(){{
    var w=document.getElementById('{wrapper_id}');
    if(!w){{setTimeout(init,60);return;}}
    var btn=w.querySelector('.cfi-trigger');
    var nm=w.querySelector('.cfi-filename');
    var inp=w.querySelector('input[type=file]');
    if(!btn||!inp){{setTimeout(init,60);return;}}
    btn.onclick=function(){{inp.click();}};
    inp.onchange=function(){{
      nm.textContent=inp.files&&inp.files.length?inp.files[0].name:'No file selected';
    }};
  }}
  init();
}})();"""
        ),
        id=wrapper_id,
        class_="compact-file-input",
    )
