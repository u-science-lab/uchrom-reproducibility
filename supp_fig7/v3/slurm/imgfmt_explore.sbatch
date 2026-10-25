#!/bin/bash
# Version 3, family imgfmt: first runs of every implementation on the built data (design of T6 / T7, choice of the
# cells of F1-F4).  Outputs to $SCRATCH/uchrom/agentbench/imgfmt_work/explore/.
#
#   sbatch --dependency=afterok:<imgfmt_setup> supp_fig7/v3/slurm/imgfmt_explore.sbatch
#SBATCH --job-name=ab3-imgfmt-explore
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=$SCRATCH/uchrom/agentbench/logs/%x-%j.out
set -uo pipefail
AB=$SCRATCH/uchrom/agentbench
GH=$HOME/uchrom/agentbench
SIF=$GH/py312full.sif
W=$AB/imgfmt_work/explore
mkdir -p "$W"
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is)"
run() {  # environment, python arguments...
  local v=$1; shift
  local t0=$(date +%s)
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=1 --env HOME=/tmp \
    --bind "$AB/data_v3:/data_v3:ro" --bind "$W:/work" --bind "$GH/src_imgfmt:/src:ro" \
    --bind "$AB/bench:/bench:ro" --bind "$GH/$v:/opt/venv:ro" --bind "${TMPDIR:-/tmp}:/tmp" \
    "$SIF" /opt/venv/bin/python -u "$@"
  echo "== exit $? after $(( $(date +%s) - t0 )) s: $v $*"
}
T6=/bench/v3/truth/imgfmt_t6.py
T7=/bench/v3/truth/imgfmt_t7.py
FM=/bench/v3/truth/imgfmt_fmt.py
URL_F1=https://uchrom-atlas-r2.u-science.org/gageseq_mouse_cortex.chromdata.zarr
URL_F4=https://uchrom-atlas-r2.u-science.org/hires_embryo.chromdata.zarr
STORE=/data_v3/imgfmt_takei2021_fov4/takei2021_fov4.chromdata.zarr

# ---- format questions first (fast)
run venv-plain $FM f1 --impl pick --url $URL_F1 --types "L4/5 IT CTX" --out /work/f1_pick.json
CELL=$(python3 -c "import json; print(json.load(open('$W/f1_pick.json'))['values']['cell'])" 2>/dev/null || echo none)
echo "F1 cell: $CELL"
run venv-uchrom $FM f1 --impl chromdata --url $URL_F1 --cell "$CELL" --chrom chr7 --bin1 30 --bin2 31 --out /work/f1_chromdata.json
run venv-plain $FM f1 --impl raw --url $URL_F1 --cell "$CELL" --chrom chr7 --bin1 30 --bin2 31 --out /work/f1_raw.json
run venv-uchrom $FM f2 --impl chromdata --store $STORE --key nucleus --cell 4_7 --out /work/f2_chromdata.json
run venv-plain $FM f2 --impl raw --store $STORE --key nucleus --cell 4_7 --out /work/f2_raw.json
run venv-uchrom $FM f3 --impl uchrom --store $STORE --key separation --cell 4_7 --chrom chr7 --out /work/f3_uchrom.json
run venv-plain $FM f3 --impl raw --store $STORE --key separation --cell 4_7 --chrom chr7 --out /work/f3_raw.json
run venv-uchrom $FM f4 --impl chromdata --url $URL_F4 --chrom chr19 --types "early neurons|blood" --out /work/f4_chromdata.json
run venv-plain $FM f4 --impl raw --url $URL_F4 --chrom chr19 --types "early neurons|blood" --out /work/f4_raw.json

# ---- T7: assignments of chr3, chr7, chr11 (201 cells) by every method; chr19 (the chromosome U-Chrom's tutorial
#      tuned its settings on) for the record
D7=/data_v3/imgfmt_seqfish
for CH in chr3,chr7,chr11 chr19; do
  tag=${CH//,/_}
  run venv-uchrom $T7 run --impl uchrom_default --data $D7 --chroms $CH --out /work/t7_${tag}_uchrom_default.json
  run venv-uchrom $T7 run --impl uchrom_tuned --data $D7 --chroms $CH --out /work/t7_${tag}_uchrom_tuned.json
  run venv-jie $T7 run --impl jie --data $D7 --chroms $CH --jobs 16 --out /work/t7_${tag}_jie_default.json
  run venv-jie $T7 run --impl jie --data $D7 --chroms $CH --jobs 16 --param num_skip=15 \
    --param total_num_skip_frac=0.8 --param init_skip_frac=0.2 --out /work/t7_${tag}_jie_nb.json
  run venv-imgfmt $T7 run --impl kmeans2 --data $D7 --chroms $CH --out /work/t7_${tag}_kmeans2.json
  run venv-imgfmt $T7 run --impl dbscan --data $D7 --chroms $CH --param eps_nm=1000 --out /work/t7_${tag}_dbscan1000.json
  run venv-imgfmt $T7 run --impl nnchain --data $D7 --chroms $CH --out /work/t7_${tag}_nnchain.json
  run venv-imgfmt $T7 run --impl all_one --data $D7 --chroms $CH --out /work/t7_${tag}_all_one.json
done

# ---- T6: the 129 allele, 1,416 haploids, 10 % of the observed loci held out
D6=/data_v3/imgfmt_sox2
run venv-imgfmt $T6 run --impl linear --data $D6 --out /work/t6_linear.json
run venv-uchrom $T6 run --impl uchrom --data $D6 --jobs 8 --out /work/t6_uchrom.json
run venv-sfimpute $T6 run --impl sfi_func --data $D6 --jobs 8 --out /work/t6_sfi_func.json
run venv-sfimpute $T6 run --impl sfi_mpi --data $D6 --jobs 8 --run-impute /src/SnapFISH-IMPUTE/run_impute.py \
  --work /tmp/t6_mpi --out /work/t6_sfi_mpi.json
ls -la "$W"
echo "end=$(date -Is)"
