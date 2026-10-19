#!/bin/bash
# Version 3, family schic, question R7 (conditional): (1) the 50-cell chr19 input (data_v3/nagano2017_dip_chr19,
# schic_build_data.py --only snaphic); (2) U-Chrom's SnapHiC (schic.py --only R7, venv-uchrom); (3) the authors'
# SnapHiC @ 782dc2b in an environment an agent would build on Python 3.12 (venv-snaphic: the newest numpy, scipy,
# pandas, scikit-image, statsmodels, networkx, h5py, cooler; install time recorded), 8 processes, defaults.
#
#   sbatch supp_fig7/v3/sherlock/schic_snaphic.sbatch
#SBATCH --job-name=ab3-schic-snaphic
#SBATCH --partition=xiaojie,normal
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -uo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
SIF=$GH/py312full.sif
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is)"
mkdir -p "$L_SCRATCH/build_tmp" "$AB/bench/v3/results/schic_work"
PIPENV="--env TMPDIR=/tmp --env PIP_CACHE_DIR=/tmp/pip --env HOME=/tmp"
W=/bench/v3/results/schic_work
run() {  # venv, rw data_v3 (0/1), python args...
  local v=$1 rw=$2; shift 2
  local d3="$AB/data_v3:/data_v3:ro"; [ "$rw" = 1 ] && d3="$AB/data_v3:/data_v3"
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=8 \
    --bind "$SCRATCH/uchrom-nd/data:/shared:ro" --bind "$AB/data:/data:ro" --bind "$d3" \
    --bind "$AB/bench:/bench" --bind "$GH/$v:/opt/venv:ro" --bind "$L_SCRATCH:/tmp" "$SIF" /opt/venv/bin/python -u "$@"
}
# (1) input
[ -f "$AB/data_v3/nagano2017_dip_chr19/cells.txt" ] || \
  run venv-uchrom 1 /bench/v3/truth/schic_build_data.py --shared /shared --data /data --out /data_v3 --only snaphic
# (2) U-Chrom
[ -n "${SKIP_UCHROM:-}" ] || run venv-uchrom 0 /bench/v3/truth/schic.py --data /data --data3 /data_v3 --shared /shared --out $W/schic_uchrom_r7.json --only R7
# (3) the authors' SnapHiC.  First run (job 47732313): the newest packages; SnapHiC fails in its interaction step
#     under pandas 3 (interaction_caller.py:253 `result.loc[:,'i'] = result['i'].astype(int)` -> TypeError: Invalid
#     value ... for dtype 'int32').  PINS (":"-separated pip specifiers, e.g. PINS='pandas<3') rebuilds the environment.
PINS=$(echo "${PINS:-}" | tr ":" " ")
if [ -n "$PINS" ] && [ "$(cat "$GH/venv-snaphic/PINS" 2>/dev/null)" != "$PINS" ]; then rm -f "$GH/venv-snaphic/OK"; fi
if [ ! -f "$GH/venv-snaphic/OK" ]; then
  rm -rf "$GH/venv-snaphic"; mkdir -p "$GH/venv-snaphic"
  t0=$(date +%s)
  apptainer exec --containall --bind "$GH/venv-snaphic:/opt/venv" --bind "$L_SCRATCH/build_tmp:/tmp" $PIPENV "$SIF" \
    bash -c "python -m venv /opt/venv && /opt/venv/bin/pip install -q numpy scipy pandas scikit-image statsmodels networkx \
             h5py cooler $(printf "'%s' " $PINS) && cd /opt/venv && git clone -q https://github.com/HuMingLab/SnapHiC src/SnapHiC && \
             cd src/SnapHiC && git checkout -q 782dc2bc7e5a9e3c376bac4a61666a1a02d061b8 && git log -1 --format='SnapHiC %H %cd'" \
    && touch "$GH/venv-snaphic/OK" && echo "$PINS" > "$GH/venv-snaphic/PINS"
  echo "venv-snaphic built in $(( $(date +%s) - t0 )) s"
fi
apptainer exec --containall --bind "$GH/venv-snaphic:/opt/venv:ro" "$SIF" /opt/venv/bin/pip freeze > "$GH/venv-snaphic.freeze.txt"
grep -i -E "^(numpy|scipy|pandas|scikit-image|statsmodels|networkx|h5py|cooler)==" "$GH/venv-snaphic.freeze.txt"
apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --bind "$SCRATCH/uchrom-nd/data:/shared:ro" --bind "$AB/data_v3:/data_v3:ro" --bind "$AB/bench:/bench" \
  --bind "$GH/venv-snaphic:/opt/venv:ro" --bind "$L_SCRATCH:/tmp" "$SIF" /opt/venv/bin/python -u \
  /bench/v3/truth/schic_ref_snaphic.py --data3 /data_v3 --shared /shared --work $W --out $W/schic_ref_snaphic.json --nproc 8
echo "end=$(date -Is)"
