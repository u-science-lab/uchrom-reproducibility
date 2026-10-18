#!/bin/bash
# Version 3, family bulk: reference answers of H1-H4, XM1, P1, P2 with U-Chrom (venv-uchrom, the agents' version),
# then the independent implementations (venv-hicref: mustache-hic, cooltools, scikit-learn; P1 / P2 in venv-uchrom
# with chromdata + numpy only, checked not to import uchrom), then the graded values -> bench/v3/results/bulk.json.
#
#   sbatch [--dependency=afterok:<bulk_setup job>] benchmarks/fig6/agent_bench/v3/sherlock/bulk_truth.sbatch [ONLY]
#   (ONLY: comma-separated question ids, default all)
#SBATCH --job-name=ab3-bulk-truth
#SBATCH --partition=xiaojie,normal
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -uo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
U=${VENV_U:-venv-uchrom-v3}      # the agents' U-Chrom environment (v3 snapshot; version 2's venv-uchrom = fb43ac4)
ONLY=${1:-H1,H2,H3,H4,XM1,P1,P2}
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} snapshot=$(cat "$AB/snapshot_v3/COMMIT") only=$ONLY start=$(date -Is)"
mkdir -p "$AB/data_v3/bulk_work" "$AB/bench/v3/results"
OUT=/bench/v3/results/bulk.json
[ -f "$AB/bench/v3/results/bulk_fb43ac4.json" ] || cp "$AB/bench/v3/results/bulk.json" "$AB/bench/v3/results/bulk_fb43ac4.json"
run() {  # environment, python arguments...
  local v=$1; shift
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=16 --env OPENBLAS_NUM_THREADS=16 --env TMPDIR=/tmp \
    --bind "$AB/data:/data:ro" --bind "$AB/data_v3:/data_v3" --bind "$AB/dl:/dl:ro" --bind "$AB/bench:/bench" \
    --bind "$GH/$v:/opt/venv:ro" --bind "$L_SCRATCH:/tmp" "$GH/py312full.sif" /opt/venv/bin/python -u "$@"
}
pick() {  # keep the ids of $ONLY that are in the list $1
  python3 -c "import sys; a=sys.argv[1].split(','); b=sys.argv[2].split(','); print(','.join(x for x in b if x in a))" "$1" "$ONLY"
}
run $U /bench/v3/truth/bulk.py --impl uchrom --only "$ONLY" --out $OUT || echo "uchrom path failed"
Q=$(pick H1,H2,H3,H4,XM1)
[ -n "$Q" ] && { run venv-hicref /bench/v3/truth/bulk.py --impl independent --only "$Q" --out $OUT || echo "independent (hicref) failed"; }
Q=$(pick P1,P2)
[ -n "$Q" ] && { run $U /bench/v3/truth/bulk.py --impl independent --only "$Q" --out $OUT || echo "independent (numpy) failed"; }
run $U /bench/v3/truth/bulk.py --impl finalize --out $OUT
echo "end=$(date -Is)"
