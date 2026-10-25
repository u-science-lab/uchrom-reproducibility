#!/bin/bash
# Precision / recall / F1 / APA of every call set of snaphic_run against the references, and the agreement
# with the authors' published calls (supp_table3/sc_structures/snaphic_mesc.py summary --apa).
#
#   sbatch --dependency=afterok:<run>:<bonev> --export=ALL,REPO=...,ENGINE=snaphic .../snaphic_summary.sbatch
#SBATCH --job-name=snaphic-summary
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --output=%x-%j.out
source "${REPRO:-$SCRATCH/uchrom/repro}/supp_table3/sc_structures/slurm/_env.sh"
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK
"$PY" "$REPRO/supp_table3/sc_structures/snaphic_mesc.py" summary --out "$OUT/snaphic/mesc" --apa
