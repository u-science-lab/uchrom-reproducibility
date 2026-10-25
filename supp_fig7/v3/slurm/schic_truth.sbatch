#!/bin/bash
# Version 3, family schic: the reference answers with U-Chrom (v3/truth/schic.py, the agents' environment venv-uchrom,
# commit fb43ac4) -> bench/v3/results/schic_work/schic_uchrom.json; also the scHiCAR cell table the independent scripts
# read (schic_cells.tsv).
#
#   sbatch --export=ALL,ONLY=R1:R2:R3,CELLS=<cellA>:<cellB> supp_fig7/v3/slurm/schic_truth.sbatch
# (lists separated by ":" -- sbatch --export splits on commas)
#SBATCH --job-name=ab3-schic-truth
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=$SCRATCH/uchrom/agentbench/logs/%x-%j.out
set -euo pipefail
AB=$SCRATCH/uchrom/agentbench
GH=$HOME/uchrom/agentbench
U=${VENV_U:-venv-uchrom-v3}      # the agents' U-Chrom environment (v3 snapshot; the first references: venv-uchrom = fb43ac4)
SIF=$GH/py312full.sif
ONLY=$(echo "${ONLY:-R1:R2:R3:R4:R5:R6}" | tr ":" ",")
CELLS=$(echo "${CELLS:-}" | tr ":" ",")
mkdir -p "$AB/bench/v3/results/schic_work"
S=$AB/bench/v3/results/schic_work      # keep the references of the fb43ac4 environment once
[ -f "$S/schic_uchrom_fb43ac4.json" ] || cp "$S/schic_uchrom.json" "$S/schic_uchrom_fb43ac4.json"
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} snapshot=$(cat "$AB/snapshot_v3/COMMIT") only=$ONLY cells=${CELLS:-} start=$(date -Is)"
run() {
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=16 --env OPENBLAS_NUM_THREADS=16 \
    --bind "$AB/data:/data:ro" --bind "$AB/data_v3:/data_v3:ro" --bind "$AB/bench:/bench" \
    --bind "$GH/$U:/opt/venv:ro" --bind "${TMPDIR:-/tmp}:/tmp" "$SIF" /opt/venv/bin/python -u "$@"
}
W=/bench/v3/results/schic_work
if [ -n "${BUILD:-}" ]; then   # BUILD=types (or other parts of schic_build_data.py, ":"-separated): data_v3 writable
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=16 --bind "$SCRATCH/uchrom/data:/shared:ro" \
    --bind "$AB/data:/data:ro" --bind "$AB/data_v3:/data_v3" --bind "$AB/bench:/bench:ro" \
    --bind "$GH/$U:/opt/venv:ro" --bind "${TMPDIR:-/tmp}:/tmp" "$SIF" /opt/venv/bin/python -u \
    /bench/v3/truth/schic_build_data.py --shared /shared --data /data --out /data_v3 --only "$(echo "$BUILD" | tr ":" ",")"
fi
run /bench/v3/truth/schic.py --data /data --data3 /data_v3 --out $W/schic_uchrom.json --export-cells $W/schic_cells.tsv
run /bench/v3/truth/schic.py --data /data --data3 /data_v3 --out $W/schic_uchrom.json --only "$ONLY" ${CELLS:+--cells "$CELLS"}
echo "end=$(date -Is)"
