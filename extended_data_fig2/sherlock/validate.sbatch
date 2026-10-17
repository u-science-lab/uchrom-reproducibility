#!/bin/bash
# uchrom.spatial against SPARK (R) and esda on the spot RNA of Guo et al. 2026 mouse brain r1 (validate_stats.py);
# the numbers go to extended_data_fig2/results/validate_stats.json of the repo copy.
#   sbatch extended_data_fig2/sherlock/validate.sbatch       (after setup.sbatch)
#SBATCH --job-name=sv-validate
#SBATCH --partition=xiaojie,normal
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=%x-%j.out
source "${REPRO:-$SCRATCH/uchrom-nd/repro}/extended_data_fig2/sherlock/_env.sh"
cd "$REPRO"
$PY -u extended_data_fig2/validate_stats.py --out "$OUT/validate"
