#!/bin/bash
# Single-cell domain boundaries (single_cell_boundaries) of the mES cells against the Bonev bulk insulation
# (supp_table3/sc_structures/impute_domains.py), then the GM12878 check (SnapHiC + boundaries on chr21 vs Rao
# 2014; supp_table3/sc_structures/snaphic_gm12878.py).
#
#   sbatch --dependency=afterok:<bonev> --export=ALL,REPO=...,ENGINE=snaphic .../snaphic_domains.sbatch
#SBATCH --job-name=snaphic-domains
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=%x-%j.out
source "${REPRO:-$SCRATCH/uchrom/repro}/supp_table3/sc_structures/slurm/_env.sh"
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK
lscpu | grep "Model name" || true
"$PY" "$REPRO/supp_table3/sc_structures/impute_domains.py" --out "$OUT/snaphic/domains" --mesc "$OUT/snaphic/mesc" \
  --n 742 200 50 --threads "$SLURM_CPUS_PER_TASK"
"$PY" "$REPRO/supp_table3/sc_structures/snaphic_gm12878.py" --out "$OUT/snaphic/gm12878"
