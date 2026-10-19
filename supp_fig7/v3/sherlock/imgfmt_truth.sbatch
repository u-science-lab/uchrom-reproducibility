#!/bin/bash
# Version 3, family imgfmt: the reference runs of T6 and F1-F4 (v3/truth/imgfmt.py describes them), into
# $SCRATCH/uchrom-nd/agentbench/imgfmt_work/truth/, then results/imgfmt.json.
#
# The authors' SnapFISH-IMPUTE (sfimpute 0.0.2 = github @7f7f1a7, 2023-11-06) fails under pandas 3 (its to_darr()
# returns a read-only array that boxcox_by1d writes into), so it runs in venv-sfimpute2023: the package with the
# releases of its time (numpy 1.26, pandas 2.1, scipy 1.11, scikit-image 0.22, scikit-learn 1.3), built here once.
#
#   sbatch supp_fig7/v3/sherlock/imgfmt_truth.sbatch            (after syncing v3/ to bench/v3/)
#   CELL=<cell id of the F2/F3 store> sbatch --export=ALL ...
#SBATCH --job-name=ab3-imgfmt-truth
#SBATCH --partition=xiaojie,normal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -uo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
U=${VENV_U:-venv-uchrom-v3}      # the agents' U-Chrom environment (v3 snapshot; the first references: venv-uchrom = fb43ac4)
SIF=$GH/py312full.sif
W=$AB/imgfmt_work/${TAG:-truth_v3}           # truth = the runs with venv-uchrom (fb43ac4), job 47751628
mkdir -p "$W" "$L_SCRATCH/build_tmp"
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is)"

# ---- the authors' environment of T6 (once)
SFI=venv-sfimpute2023
if ! apptainer exec --containall --cleanenv --bind "$GH/$SFI:/opt/venv:ro" "$SIF" \
     /opt/venv/bin/python -c "import sfimpute.impute, mpi4py, pandas; assert pandas.__version__.startswith('2.1')" 2>/dev/null; then
  rm -rf "$GH/$SFI"; mkdir -p "$GH/$SFI"
  apptainer exec --containall --bind "$GH/$SFI:/opt/venv" --bind "$L_SCRATCH/build_tmp:/tmp" \
    --env TMPDIR=/tmp --env PIP_CACHE_DIR=/tmp/pip --env HOME=/tmp "$SIF" bash -c \
    "python -m venv /opt/venv && /opt/venv/bin/pip install -q --upgrade pip && /opt/venv/bin/pip install -q \
     sfimpute==0.0.2 numpy==1.26.4 pandas==2.1.4 scipy==1.11.4 scikit-image==0.22.0 scikit-learn==1.3.2 \
     matplotlib==3.8.2 seaborn==0.13.0 mpi4py mpich"
  apptainer exec --containall --bind "$GH/$SFI:/opt/venv" "$SIF" /opt/venv/bin/pip freeze > "$GH/$SFI.freeze.txt"
fi
grep -i -E "^(sfimpute|numpy|pandas|scipy|scikit-image|scikit-learn|mpi4py|mpich)==" "$GH/$SFI.freeze.txt"

run() {  # environment, python arguments...
  local v=$1; shift
  local t0=$(date +%s)
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=1 --env HOME=/tmp \
    --bind "$AB/data_v3:/data_v3:ro" --bind "$W:/work" --bind "$GH/src_imgfmt:/src:ro" \
    --bind "$AB/bench:/bench:ro" --bind "$GH/$v:/opt/venv:ro" --bind "$L_SCRATCH:/tmp" \
    "$SIF" /opt/venv/bin/python -u "$@"
  echo "== exit $? after $(( $(date +%s) - t0 )) s: $v $*"
}
T6=/bench/v3/truth/imgfmt_t6.py
FM=/bench/v3/truth/imgfmt_fmt.py
URL_F1=https://uchrom-atlas-r2.u-science.org/gageseq_mouse_cortex.chromdata.zarr
URL_F4=https://uchrom-atlas-r2.u-science.org/hires_embryo.chromdata.zarr
STORE=/data_v3/imgfmt_takei2021_fov4/takei2021_fov4.chromdata.zarr

# ---- F1-F4 (each twice: chromdata / U-Chrom, and the raw Zarr + Parquet readers in the plain environment)
run venv-plain $FM f1 --impl pick --url $URL_F1 --types "L4/5 IT CTX" --out /work/f1_pick.json
C1=$(python3 -c "import json; print(json.load(open('$W/f1_pick.json'))['values']['cell'])")
echo "F1 cell: $C1"
run $U $FM f1 --impl chromdata --url $URL_F1 --cell "$C1" --chrom chr7 --bin1 30 --bin2 31 --out /work/f1_chromdata.json
run venv-plain $FM f1 --impl raw --url $URL_F1 --cell "$C1" --chrom chr7 --bin1 30 --bin2 31 --out /work/f1_raw.json
C2=${CELL:-$(apptainer exec --containall --cleanenv --bind "$AB/data_v3:/data_v3:ro" --bind "$GH/$U:/opt/venv:ro" \
  "$SIF" /opt/venv/bin/python -c "
from chromdata import ChromData
cd = ChromData.read('$STORE', backed=True)
ids = sorted(map(str, cd.cells.index)); print(ids[len(ids) // 2])" 2>/dev/null)}
echo "F2/F3 cell: $C2"
run $U $FM f2 --impl chromdata --store $STORE --key nucleus --cell "$C2" --out /work/f2_chromdata.json
run venv-plain $FM f2 --impl raw --store $STORE --key nucleus --cell "$C2" --out /work/f2_raw.json
run $U $FM f3 --impl uchrom --store $STORE --key separation --cell "$C2" --chrom chr7 --out /work/f3_uchrom.json
run venv-plain $FM f3 --impl raw --store $STORE --key separation --cell "$C2" --chrom chr7 --out /work/f3_raw.json
run $U $FM f4 --impl chromdata --url $URL_F4 --chrom chr19 --types "early neurons|blood" --out /work/f4_chromdata.json
run venv-plain $FM f4 --impl raw --url $URL_F4 --chrom chr19 --types "early neurons|blood" --out /work/f4_raw.json

# ---- T6: the 129 allele, 1,416 haploids, 10 % of the observed loci held out
D6=/data_v3/imgfmt_sox2
run venv-imgfmt $T6 run --impl linear --data $D6 --out /work/t6_linear.json
run $U $T6 run --impl uchrom --data $D6 --jobs 16 --out /work/t6_uchrom.json
run $SFI $T6 run --impl sfi_func --data $D6 --jobs 16 --out /work/t6_sfi_func.json
run $SFI $T6 run --impl sfi_mpi --data $D6 --jobs 16 --run-impute /src/SnapFISH-IMPUTE/run_impute.py \
  --work /tmp/t6_mpi --out /work/t6_sfi_mpi.json
ls -la "$W"
echo "end=$(date -Is)"
