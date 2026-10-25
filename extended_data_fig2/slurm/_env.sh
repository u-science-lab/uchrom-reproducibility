# Shared settings of the spatial_variable jobs on a compute cluster (sourced by every *.sbatch here).
#
#   $HOME/uchrom/venv-sc         Python env of supp_table3/sc_structures/slurm (+ esda / libpysal,
#                                               the reference of Moran's I / Geary's C: setup.sbatch)
#   $HOME/uchrom/ref/Rlib        R library with SPARK (the reference of SPARK-X; setup.sbatch)
#   $SCRATCH/uchrom/repo-sv                  repo copy synced from the laptop (tracked files + REVISION)
#   $UCHROM_ARCHIVE/spatial_3dgenome/<sample>/      per-spot contact statistics (contact_stats.sbatch of the recipe)
#   $SCRATCH/uchrom/runs/sv                  outputs
set -euo pipefail
ml purge
ml python/3.12.1 R/4.4.2
BASE=$HOME/uchrom
W=$SCRATCH/uchrom
REPO=${REPO:-$W/repo-sv}
REPRO=${REPRO:-$W/repro}       # a copy of uchrom-reproducibility: the benchmarks (REPO: U-Chrom, the packages)
export UCHROM_REPO=$REPO           # paper_paths.UCHROM_REPO: the benchmarks import its packages
PY=$BASE/venv-sc/bin/python
export R_LIBS_USER=$BASE/ref/Rlib
export UCHROM_DATA=$W/data
export PIP_CACHE_DIR=$SCRATCH/.cache/pip XDG_CACHE_HOME=$SCRATCH/.cache MPLBACKEND=Agg
export PYTHONPATH=$REPO/packages/uchrom:$REPO/packages/chromdata:$W/engine-sc/${ENGINE:-final}${PYTHONPATH:+:$PYTHONPATH}
STATS=$UCHROM_ARCHIVE/spatial_3dgenome
# the code revision as the job starts (shared/_record.py records it; a later sync does not change it)
export UCHROM_REVISION=$(cat "$REPO/REVISION" 2>/dev/null || echo unknown)
OUT=$W/runs/sv
mkdir -p "$OUT"
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} repo=$REPO revision=$(cat "$REPO/REVISION" 2>/dev/null || echo ?)"
