#!/usr/bin/env bash
# deploy-blank.sh — template for pushing TAGSITES to shinyapps.io
#
# This is the tracked, account-free copy of the deploy script. The working copy, deploy.sh,
# is gitignored because it carries your account name and app id. To set one up:
#
#   1. Copy it:            cp deploy-blank.sh deploy.sh && chmod +x deploy.sh
#
#   2. Install the client in the Python 3.10 environment you deploy from (shinyapps.io
#      provisions the interpreter version it sees, and rs3 caps Python at 3.10 — see
#      CLAUDE.md "Version pinning"):
#          pip install rsconnect-python
#
#   3. Register your shinyapps.io account once per machine. On shinyapps.io go to
#      Account -> Tokens -> Show, then:
#          rsconnect add --account <ACCOUNT> --name <ACCOUNT> --token <TOKEN> --secret <SECRET>
#      The token and secret stay in rsconnect's own config, never in this script.
#
#   4. Fill in the two values below (or export them as environment variables):
#        RSCONNECT_NAME    the server nickname given to --name in step 3 (usually the
#                          account name, the <ACCOUNT> in https://<ACCOUNT>.shinyapps.io/)
#        RSCONNECT_APP_ID  the numeric id of an existing app, to update it in place. Leave
#                          it empty for the first deploy: a new app is created and its id
#                          is printed (also visible in the app's URL on the shinyapps.io
#                          dashboard). Put that id here so later deploys update the same
#                          app instead of creating another one.
#
#   5. Optional: a UCSC BLAT API key in ucsc.local.json (gitignored) enables the fast
#      off-target region screen; without it the app falls back to EBI. See the NOTE in
#      EXCLUDES below: that file must ship inside the bundle.
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
# is a few MB.

set -euo pipefail

# ── fill these in (see step 4 above), or export them before running ──────────
RSCONNECT_NAME="${RSCONNECT_NAME:-}"        # e.g. myaccount
RSCONNECT_APP_ID="${RSCONNECT_APP_ID:-}"    # e.g. 12345678; empty = create a new app

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
  -x 'deploy-blank.sh'
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

if [[ -z "$RSCONNECT_NAME" ]]; then
  echo "Set RSCONNECT_NAME (the server nickname from 'rsconnect add'); see the header." >&2
  exit 1
fi

# with an app id the existing app is updated in place; without one a new app is created
APP_ARGS=()
if [[ -n "$RSCONNECT_APP_ID" ]]; then
  APP_ARGS=(-a "$RSCONNECT_APP_ID")
else
  echo "RSCONNECT_APP_ID is empty: deploying as a NEW app. Record the id it prints." >&2
fi
rsconnect deploy shiny . -n "$RSCONNECT_NAME" ${APP_ARGS[@]+"${APP_ARGS[@]}"} --entrypoint "$ENTRYPOINT" \
  "${EXCLUDES[@]}"
