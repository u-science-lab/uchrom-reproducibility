#!/bin/bash
# Linked contact maps (Supplementary Fig. 1c-f): linking
# against copying the scHiCAR per-cell .scool, and the web browser's pseudo-bulk queries, on the scHiCAR store with
# its maps linked only, opened backed as the browser opens it (e_contacts_query.py, e_contacts.py).  Replaces the
# 2026-09-26 laptop runs, whose queries read the legacy .h5cd file and whose link / copy stores were format 2.1.
# Data: $SCRATCH/uchrom/data/schicar_mop (synced from the laptop's data folder: the store, the per-cell .scool,
# the 22 per-type .cool, the RNA .h5ad), staged on node-local disk; results to $W/runs/fig2/linked_<job>/.
#   sbatch extended_data_fig1/slurm/linked.sbatch
#SBATCH --job-name=fig2-linked
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=%x-%j.out
export REPO=${REPO:-$SCRATCH/uchrom/repo-linked}
source "$REPRO/extended_data_fig1/slurm/_env.sh"
# cooler reads the published per-type maps (the identity check); pip in the job, not on a login node
$PY -c "import cooler" 2>/dev/null || $PY -m pip install -q --prefer-binary "cooler==0.10.4" \
    "numpy==$($PY -c 'import numpy; print(numpy.__version__)')" \
    "pandas==$($PY -c 'import pandas; print(pandas.__version__)')"
LOCAL=${TMPDIR:-/tmp}/data
mkdir -p "$LOCAL"
cp -r "$UCHROM_DATA/schicar_mop" "$LOCAL/"              # node-local disk, as the laptop's SSD
export UCHROM_DATA=$LOCAL
RUN=$OUT/linked_${SLURM_JOB_ID:-local}
export UCHROM_REPRO=$RUN/repro                          # results_dir(): $UCHROM_REPRO/extended_data_fig1/results
mkdir -p "$UCHROM_REPRO/extended_data_fig1/results"
cd "$REPRO"
$PY -u extended_data_fig1/e_contacts_query.py --repeats 3
$PY -u extended_data_fig1/e_contacts.py --repeats 3
ls -la "$UCHROM_REPRO/extended_data_fig1/results"
