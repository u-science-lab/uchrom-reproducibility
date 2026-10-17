#!/bin/bash
# The unit tests of chromdata, u-chrom and the atlas on the repo copy (venv-sc; the laptop only runs small tests).
#SBATCH --job-name=sv-tests
#SBATCH --partition=xiaojie,normal
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=%x-%j.out
source "${REPRO:-$SCRATCH/uchrom-nd/repro}/extended_data_fig2/sherlock/_env.sh"
cd "$REPO/packages/chromdata" && $PY -m pytest tests -q -p no:cacheprovider 2>&1 | tail -5
cd "$REPO/packages/uchrom" && $PY -m pytest tests -q -p no:cacheprovider 2>&1 | tail -15
cd "$REPO/apps/atlas" && $PY -m pytest tests -q -p no:cacheprovider 2>&1 | tail -5
