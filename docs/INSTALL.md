# Installing TAGSITES

Use this once per machine. Every other guide assumes the `tagsites` environment exists.

## 1. Prerequisites

- macOS or Linux (Windows works for the web-services install only).
- [mamba](https://mamba.readthedocs.io/) or conda. Mamba is faster and is used in these docs;
  replace `mamba` with `conda` if that is what is installed.
- Internet access. Local workflows need it only for one-time reference downloads.
- Disk: a web-services install needs well under 1 GB. Local backends add about 9 GB for one
  organism ([LOCAL_SETUP.md](LOCAL_SETUP.md)); a database build adds tens of GB
  ([BUILD_LOCAL_DATABASE.md](BUILD_LOCAL_DATABASE.md)). Put reference data on **local disk**,
  because SQLite indexes are unreliable on network volumes.

## 2. Get the code

```bash
git clone https://github.com/msr2009/TAGSITES.git
cd TAGSITES
```

Without git: click **⬇ Download App** in the running app, then
`unzip tagsites-app.zip && cd tagsites-app`. The zip includes `environment.yml`, `docs/` and
the command-line entry points but not `tests/`.

## 3. Create the environment

```bash
mamba env create -f environment.yml
mamba activate tagsites
```

`environment.yml` installs the Python packages plus the local-backend tools: `diamond`, `mafft`
and BLAST+ (`blastn`, `makeblastdb`, `blastdbcmd`). `requirements.txt` is the shinyapps.io
manifest and lacks those tools, so use it only for a web-services install:

```bash
python3.10 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### Version pins

Python 3.10, `numpy=1.26.*`, `scikit-learn=1.0.2` and `lightgbm=3.3.5` are deliberate. All four
exist only because `rs3` (Rule Set 3 guide scoring) caps them. RS3 fails soft: without it the
pipeline writes blank `rs3_score` columns and everything else is unchanged. If a new dependency
needs a newer Python or numpy 2, drop `rs3` instead of changing the pins (see "Version pinning"
in [CLAUDE.md](../CLAUDE.md)).

## 4. Check the install

```bash
python -c "import shiny, Bio, plotly, pandas; print('python packages ok')"
python -c "import rs3" && echo "RS3 available"        # optional; guide scores only
which diamond mafft blastn makeblastdb blastdbcmd      # local backends only
python -c "import pyhmmer; print(pyhmmer.__version__)" # local domain scanning only
pytest                                                  # offline test suite
```

Missing `diamond`, `mafft`, `pyhmmer` or `blastn` does not affect the web-services app.

## 5. Accounts and keys

| Item | Needed for | Where |
|---|---|---|
| E-mail address | Every EBI job (BLAST, InterPro, Genewise, Clustal). The Setup tab requires one even where no EBI job is submitted. | Setup tab, `Email` field |
| UCSC BLAT API key | Fast off-target screening through UCSC BLAT when no local BLAST+ database exists | `ucsc.local.json` (`{"api_key": "..."}`, gitignored). Without it the screen falls back to EBI blastn. |
| DeepTMHMM checkout | Topology prediction for database builds only (separate academic license) | `deeptmhmm` block of the config ([BUILD_LOCAL_DATABASE.md](BUILD_LOCAL_DATABASE.md)) |

## 6. Start the app

```bash
mamba activate tagsites
python app_modular.py          # http://127.0.0.1:8000
```

Follow [LOCAL_SETUP.md](LOCAL_SETUP.md), [LOCAL_DATABASE.md](LOCAL_DATABASE.md) or
[BATCH_APP.md](BATCH_APP.md) for local modes.

## 7. Update

```bash
git pull
mamba env update -f environment.yml --prune
```

Restart any running app afterwards. The batch app changes often, so a downloaded copy that is
more than a few weeks old may not open a current database.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| Solver conflict creating the environment | Another package wants numpy 2 or Python 3.11+. Keep the pins; install the package with `pip` only if it does not need to change them. |
| `RS3 unavailable` in logs | `rs3` or `lightgbm` is missing or the wrong version. Reinstall from `environment.yml`, or accept blank RS3 scores. |
| Port 8000 in use | Another app instance is running. Stop it, or start on another port with `python -m shiny run --port 8010 app_modular:app`. |
