#!/bin/bash
# Stages of spatial_variable.py (cache: $OUT/cache; figure data: $OUT/figdata, copied to the paper repository
# from the laptop).    sbatch analysis.sbatch regions,compartments [brain | melanoma]
# (stages as the argument: sbatch --export splits a value at its commas)
#SBATCH --job-name=sv-analysis
#SBATCH --partition=xiaojie,normal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=%x-%j.out
source "${REPRO:-$SCRATCH/uchrom-nd/repro}/extended_data_fig2/sherlock/_env.sh"
cd "$REPRO"
$PY -u extended_data_fig2/spatial_variable.py --tissue "${2:-brain}" --stats "$STATS" --out "$OUT/cache" \
    --figdata "$OUT/figdata" \
    --stages "${1:?stages}" --workers "${SLURM_CPUS_PER_TASK:-8}"
