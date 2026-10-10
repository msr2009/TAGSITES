#!/usr/bin/env bash
# deploy.sh — push TAGSITES to shinyapps.io (https://mattrich.shinyapps.io/tagsites/)
#
# Run from the repo root:   ./deploy.sh
# Dry run (writes manifest.json, uploads nothing):   ./deploy.sh --manifest
#
# WHY THIS SCRIPT EXISTS, AND WHY .rscignore DOES NOT WORK
#
# .rscignore is an R rsconnect feature. rsconnect-python does NOT read it — the
# string appears nowhere in the package (verified against 1.31.1), and a test
# bundle built with a .rscignore listing two paths included both of them, plus
# the .rscignore file itself. Every exclusion must be passed as -x/--exclude.
#
# rsconnect-python excludes only these on its own: .git/ .svn/ __pycache__/
# node_modules/ renv/ packrat/ rsconnect/ rsconnect-python/ .Rproj.user/
# (bundle.py:69-80). Everything else below has to be named explicitly.
#
# Without these flags the bundle is ~10 GB, nearly all of it data/. With them it
# is about 6 MB.

set -euo pipefail

ENTRYPOINT="app_modular:app"

EXCLUDES=(
  # ── Large directories that must never ship ──────────────────────────────
  -x 'data/'                       # 9.7 GB of run outputs
  -x 'USERDATA/'
  -x 'webservice-clients/'         # 188 MB vendored EBI client repo, 176 MB of
                                   # it Java .jars. Reference only: the app uses
                                   # scripts/ncbiblast.py etc. and clients_folder
                                   # points at scripts/, never here.
  -x 'notebooks/'                  # 7.7 MB dev notebooks
  -x '.ipynb_checkpoints/'         # 5.7 MB Jupyter autosaves
  -x 'scripts/uniprot_species/'    # 7.5 MB species flat file + notebook
  -x 'uniprot_species.flat.txt'
  -x 'yeast_tag_rescue_analysis/'  # 15 MB yeast tagging comparison, analysis only
  -x 'WORMPRO/'                    # local exploration tables

  # NOTE: ucsc.local.json is gitignored but MUST ship — shinyapps.io has no
  # secrets mechanism, so the UCSC BLAT API key travels inside this private
  # bundle (scripts/offtarget_blat.py). rsconnect does not read .gitignore, so
  # it is included automatically; do not add an -x for it.

  # ── Dev-only ────────────────────────────────────────────────────────────
  -x 'tests/'
  -x 'CLAUDE.md'
  -x 'environment.yml'             # requirements.txt is the deploy manifest
  -x '.claude/'
  -x 'docs/'
  -x 'deploy.sh'
  -x 'PLAN.md'

  # ── Caches and editor leftovers ─────────────────────────────────────────
  -x '.pytest_cache/'
  -x '.ruff_cache/'
  # rsconnect's own __pycache__ rule only matches at the top level, so nested
  # ones (modules/, scripts/, utils/) leak: 133 .pyc files, 73 of them stale
  # cpython-312 from before the Python 3.10 downgrade. Only the '**/*.pyc' form
  # catches them — '**/__pycache__/' and '*/__pycache__/' both fail to.
  -x '**/*.pyc'
  -x '*~'                          # vim backups AND .*.un~ undo files
  -x '.rscignore'
  -x '.smbdelete*'
  -x '.DS_Store'

  # ── Scripts the running app never imports ───────────────────────────────
  # Kept deliberately small. scripts/ is only ~2 MB in total, and
  # scripts/providers.py:89 resolves backends via importlib.import_module, so
  # *_remote.py / *_local.py / *_bulk.py are invisible to static analysis and
  # must NOT be excluded on the strength of an import scan. Only entrypoints
  # confirmed unreachable are listed.
  -x 'scripts/batch_run_tag_sites.py'
  -x 'scripts/batch_plot_tagsites.py'
  -x 'scripts/batch_rescore_tagsites.py'
  -x 'scripts/proteome_run.py'
  -x 'scripts/build_score_testing_set.py'
  -x 'scripts/build_heatmap_images.py'
  -x 'scripts/local_AF2_prediction.py'
  -x 'scripts/design_guides_across_region.py'
  -x 'scripts/get_species_taxonomy.py'
  -x 'scripts/run_tag_sites.py'
  -x 'scripts/run_tag_sites_from_json.py'
  -x 'scripts/dpy-5.guides'
  -x 'scripts/tmp*'
  -x 'scripts/distributions/'
)

if [[ "${1:-}" == "--manifest" ]]; then
  # dry run: list what would be bundled, then clean up. A stray manifest.json
  # changes how later deploys resolve dependencies, so it must not be left behind.
  rsconnect write-manifest shiny . --entrypoint "$ENTRYPOINT" --overwrite "${EXCLUDES[@]}"
  python - <<'PY'
import json, os
m = json.load(open("manifest.json"))
files = sorted(m["files"])
total = sum(os.path.getsize(f) for f in files if os.path.exists(f))
print("\n{} files, {:.1f} MB".format(len(files), total / 1e6))
for f in files:
    if os.path.exists(f) and os.path.getsize(f) > 500_000:
        print("  large: {:.1f} MB  {}".format(os.path.getsize(f) / 1e6, f))
PY
  rm -f manifest.json
  exit 0
fi

# app_id comes from rsconnect-python/TAGSITES.json, so this updates the existing
# app in place rather than creating a second one
rsconnect deploy shiny . --entrypoint "$ENTRYPOINT" "${EXCLUDES[@]}"
