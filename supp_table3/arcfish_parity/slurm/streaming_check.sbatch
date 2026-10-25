#!/bin/bash
# The ArcFISH part of the streaming validation (Extended Data Fig. 1e, Supplementary Note spec) re-run with U-Chrom
# origin/main b6d8873 (venv-uchrom-v3: wheels of that snapshot, after the ArcFISH parity re-implementation de8384e).
# Copy streaming_check.py into $W first (W below); results/streaming.json is the output of job 47904579.
#SBATCH --job-name=arcfish-stream
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=$SCRATCH/uchrom/arcfish_parity/logs/%x-%j.out
set -uo pipefail
AP=$SCRATCH/uchrom/arcfish_parity
W=$AP/streaming
GH=$HOME/uchrom/agentbench
D=$SCRATCH/uchrom/data
REV=$(cat $SCRATCH/uchrom/agentbench/snapshot_v3/COMMIT)
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} rev=$REV start=$(date -Is)"
mkdir -p "${TMPDIR:-/tmp}/tmp"
apptainer exec --containall --cleanenv --bind "$W:$W" --bind "$D:$D:ro" --bind "${TMPDIR:-/tmp}/tmp:/tmp" \
  --bind "$GH/venv-uchrom-v3:/opt/venv:ro" --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --env UCHROM_REVISION="$REV" --env HOME=/tmp \
  "$GH/py312full.sif" /opt/venv/bin/python -u "$W/streaming_check.py" "$D/4DNFIHF3JCBY.csv" \
  "$W/streaming_${SLURM_JOB_ID}.json"
echo "exit=$? end=$(date -Is)"
md5sum "$D/4DNFIHF3JCBY.csv"
