#!/bin/bash
# Run time / memory / agreement of the RWR imputation and SnapHiC's sums, numpy vs native, by threads
# (supp_table3/sc_structures/impute_speed.py).  One whole node's worth of cores for the thread scaling.
#
#   sbatch --export=ALL,REPO=$SCRATCH/uchrom-nd/repo-sc-snaphic,ENGINE=snaphic .../impute_speed.sbatch
#SBATCH --job-name=impute-speed
#SBATCH --partition=xiaojie,normal
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=%x-%j.out
source "${REPRO:-$SCRATCH/uchrom-nd/repro}/supp_table3/sc_structures/sherlock/_env.sh"
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK
lscpu | grep -E "Model name|^CPU\(s\)|Thread|Socket" || true
"$PY" "$REPRO/supp_table3/sc_structures/impute_speed.py" all --out "$OUT/snaphic/speed" --threads 1 4 16 32
