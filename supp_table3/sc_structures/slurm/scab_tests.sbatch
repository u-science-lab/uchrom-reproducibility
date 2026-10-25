#!/bin/bash
# The uchrom unit tests on the scab repo copy, in venv-sc and in the env that produced the scab numbers
# (SCAB_PY, e.g. venv-consensus); every failure listed (no -x).
#
#   sbatch --export=ALL,REPO=$SCRATCH/uchrom/repo-sc-scab,SCAB_PY=... supp_table3/sc_structures/slurm/scab_tests.sbatch
#SBATCH --job-name=scab-tests
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=%x-%j.out
source "${REPRO:-$SCRATCH/uchrom/repro}/supp_table3/sc_structures/slurm/_env.sh"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$REPO/packages/uchrom"
for P in "$PY" ${SCAB_PY:+"$SCAB_PY"}; do
  echo "=== $P"
  "$P" -c 'import numpy, pandas, scipy, cooler; print("numpy", numpy.__version__, "pandas", pandas.__version__, "scipy", scipy.__version__, "cooler", cooler.__version__)'
  "$P" -m pytest tests -q -p no:cacheprovider -rf 2>&1 | grep -E "^(FAILED|ERROR)|passed|failed" || true
done
