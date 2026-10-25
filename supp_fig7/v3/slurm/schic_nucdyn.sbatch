#!/bin/bash
# Version 3, family schic: the authors' NucDynamics (github.com/tjs23/nuc_dynamics @ 00eaf0b, venv-nucdyn, no U-Chrom)
# on Stevens 2017 Cell 1 (R1) and Cell 2 (R2), 1-Mb final particles, seeds 0 and 1, on 8 CPUs (an agent's session);
# runtimes recorded (v3/truth/schic_ref_nucdyn.py).
#
#   sbatch [--export=ALL,SEEDS=0:1,CELLSEL=R1:R2,EXTRA="-temps 100",TAG=_t100] supp_fig7/v3/slurm/schic_nucdyn.sbatch
#SBATCH --job-name=ab3-schic-nucdyn
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=$SCRATCH/uchrom/agentbench/logs/%x-%j.out
set -euo pipefail
AB=$SCRATCH/uchrom/agentbench
GH=$HOME/uchrom/agentbench
SIF=$GH/py312full.sif
SEEDS=$(echo "${SEEDS:-0:1}" | tr ":" ",")
EXTRA=${EXTRA:-}            # e.g. EXTRA="-temps 100" (a shorter anneal); results under <label>${TAG}
TAG=${TAG:-}
CELLSEL=${CELLSEL:-R1:R2}
mkdir -p "$AB/bench/v3/results/schic_work"
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is) seeds=$SEEDS"
run() {
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=1 \
    --bind "$AB/data:/data:ro" --bind "$AB/data_v3:/data_v3:ro" --bind "$AB/bench:/bench" \
    --bind "$GH/venv-nucdyn:/opt/venv:ro" --bind "${TMPDIR:-/tmp}:/tmp" "$SIF" /opt/venv/bin/python -u "$@"
}
W=/bench/v3/results/schic_work
case ":$CELLSEL:" in *:R1:*) run /bench/v3/truth/schic_ref_nucdyn.py --cell-dir /data/stevens2017_cell1 --label "R1$TAG" \
  --out $W/schic_ref_nucdyn$TAG.json --seeds "$SEEDS" --extra "$EXTRA";; esac
case ":$CELLSEL:" in *:R2:*) run /bench/v3/truth/schic_ref_nucdyn.py --cell-dir /data_v3/stevens2017_cell2 --label "R2$TAG" \
  --out $W/schic_ref_nucdyn$TAG.json --seeds "$SEEDS" --extra "$EXTRA";; esac
echo "end=$(date -Is)"
