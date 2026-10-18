#!/bin/bash
# Version 3, family schic: the independent references of R3-R6 on 8 CPUs (what an agent's session has):
#   the authors' scHiCluster (venv-schicluster; v3/truth/schic_ref_schicluster.py) for R3, R4, R5 (compartment score),
#   cooltools / numpy (venv-verify; v3/truth/schic_ref_numpy.py) for R5 (E1, correlation) and R6.
# No U-Chrom in either environment.  Needs schic_cells.tsv from schic_truth.sbatch.
#
#   sbatch --export=ALL,CELLS=<cellA>:<cellB>[,ONLY_S=R3:R4:R5][,ONLY_N=R5:R6] benchmarks/fig6/agent_bench/v3/sherlock/schic_ref.sbatch
# (lists separated by ":" -- sbatch --export splits on commas; ONLY_S=none / ONLY_N=none skip a part)
#SBATCH --job-name=ab3-schic-ref
#SBATCH --partition=xiaojie,normal
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -euo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
SIF=$GH/py312full.sif
ONLY_S=$(echo "${ONLY_S:-R3:R4:R5}" | tr ":" ","); [ "$ONLY_S" = none ] && ONLY_S=""
ONLY_N=$(echo "${ONLY_N:-R5:R6}" | tr ":" ","); [ "$ONLY_N" = none ] && ONLY_N=""
CELLS=$(echo "${CELLS:-}" | tr ":" ",")
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} cells=${CELLS:-} start=$(date -Is)"
run() {  # venv, threads, python args...
  local v=$1 t=$2; shift 2
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=$t --env OPENBLAS_NUM_THREADS=$t --env MKL_NUM_THREADS=$t \
    --bind "$AB/data:/data:ro" --bind "$AB/data_v3:/data_v3:ro" --bind "$SCRATCH/uchrom-nd/data:/shared:ro" \
    --bind "$AB/bench:/bench" --bind "$GH/$v:/opt/venv:ro" --bind "$L_SCRATCH:/tmp" "$SIF" /opt/venv/bin/python -u "$@"
}
W=/bench/v3/results/schic_work
if [ -n "$ONLY_S" ]; then
  run venv-schicluster 1 /bench/v3/truth/schic_ref_schicluster.py --data /data --data3 /data_v3 \
    --cells-tsv $W/schic_cells.tsv --out $W/schic_ref_schicluster.json --only "$ONLY_S" --r3-cells "${CELLS:-GTGAAATCTACGAGATAC,TAGACCGGTATGTACCCG}" --workers 8
fi
if [ -n "$ONLY_N" ]; then
  run venv-verify 8 /bench/v3/truth/schic_ref_numpy.py --data /data --data3 /data_v3 --shared /shared --work $W \
    --out $W/schic_ref_numpy.json --only "$ONLY_N"
fi
echo "end=$(date -Is)"
