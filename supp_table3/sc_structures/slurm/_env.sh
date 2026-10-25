# Shared settings of the single-cell structure-calling jobs on a compute cluster (sourced by every *.sbatch here).
#
#   $HOME/uchrom/venv-sc          Python env: chromdata + u-chrom (editable, repo-sc) and their
#                                                dependencies, pytest (setup.sbatch)
#   $SCRATCH/uchrom/repo-sc[-<part>]           repo copies, synced from the laptop (tracked + new files only:
#                                                `git describe --always --dirty=+dirty --abbrev=12 > REVISION;
#                                                (git ls-files -z --cached --others --exclude-standard; printf 'REVISION\0')
#                                                | rsync -a --from0 --files-from=- ./ <cluster>:...`; REVISION is
#                                                what shared/_record.py records as the code revision);
#                                                REPO= picks one (default repo-sc)
#   $SCRATCH/uchrom/engine-sc/<tag>           uchrom_maps + uchrom_recon built from a repo copy (setup.sbatch)
#   $SCRATCH/uchrom/data                      UCHROM_DATA (shared with the other uchrom jobs)
#   $SCRATCH/uchrom/runs/sc                   outputs
set -euo pipefail
ml purge
ml python/3.12.1
BASE=$HOME/uchrom
W=$SCRATCH/uchrom
REPO=${REPO:-$W/repo-sc}
REPRO=${REPRO:-$W/repro}       # a copy of uchrom-reproducibility: the benchmarks (REPO: U-Chrom, the packages)
export UCHROM_REPO=$REPO
PY=$BASE/venv-sc/bin/python
ENGINE=${ENGINE:-sc}
export UCHROM_DATA=$W/data
export PIP_CACHE_DIR=$SCRATCH/.cache/pip XDG_CACHE_HOME=$SCRATCH/.cache MPLBACKEND=Agg
# the repo copy's packages first, then its native kernels
export PYTHONPATH=$REPO/packages/uchrom:$REPO/packages/chromdata:$W/engine-sc/$ENGINE${PYTHONPATH:+:$PYTHONPATH}
# the code revision as the job starts (shared/_record.py records it; a later sync does not change it)
export UCHROM_REVISION=$(cat "$REPO/REVISION" 2>/dev/null || echo unknown)
OUT=$W/runs/sc
mkdir -p "$OUT"
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} repo=$REPO engine=$ENGINE $($PY -c 'import sys; print(sys.version.split()[0])')"
