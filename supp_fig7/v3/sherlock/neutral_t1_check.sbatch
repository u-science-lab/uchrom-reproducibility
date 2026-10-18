#!/bin/bash
# Version 3: T1 from the FOF-CT export (v3/data/neutral_t1_check.py): the input built from the core table with pandas,
# then the authors' ArcFISH (venv-arcfish, a21d889) and U-Chrom (venv-uchrom-v3) on it, compared with the reference.
#
#   sbatch benchmarks/fig6/agent_bench/v3/sherlock/neutral_t1_check.sbatch      (after neutral_build.sbatch)
#SBATCH --job-name=ab3-neutral-t1
#SBATCH --partition=xiaojie,normal
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -uo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
AP=$SCRATCH/uchrom-nd/arcfish_parity
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
SIF=$GH/py312full.sif
REPO=$SCRATCH/uchrom-nd/repo-main
REV=$(cat "$REPO/REVISION")
W=$AP/neutral_t1_${SLURM_JOB_ID}
mkdir -p "$W" "$L_SCRATCH/tmp"
echo "host=$(hostname) job=$SLURM_JOB_ID revision=$REV start=$(date -Is)"
COMMON=(--containall --cleanenv --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4
        --env UCHROM_REVISION="$REV" --env SLURM_JOB_ID="$SLURM_JOB_ID"
        --bind "$L_SCRATCH/tmp:/tmp" --bind "$AP:$AP" --bind "$REPO:$REPO:ro" --bind "$AB:$AB:ro")
ours() { apptainer exec "${COMMON[@]}" --bind "$GH/venv-uchrom-v3:/opt/venv:ro" --env PYTHONPATH="$REPO/benchmarks" \
  "$SIF" /opt/venv/bin/python -u "$@"; }
theirs() { apptainer exec "${COMMON[@]}" --bind "$GH/venv-arcfish:/opt/venv:ro" --env PYTHONPATH=/opt/venv/src \
  "$SIF" /opt/venv/bin/python -u "$@"; }
B=$REPO/benchmarks/arcfish_parity
C=$AB/bench/v3/data/neutral_t1_check.py
ours "$C" prepare --core "$AB/data_v3/takei2025_fov0_fofct/takei2025_fov0_fofct_core.csv" \
  --store "$AB/data/takei2025_fov0/takei2025_fov0.chromdata.zarr" --manifest "$AP/inputs_$REV/manifest.json" --out "$W"
theirs "$B/run_upstream.py" --manifest "$W/manifest.json" --name takei2025_chr15_1mb --out "$W" \
  || echo "upstream FAILED"
ours "$B/run_ours.py" --manifest "$W/manifest.json" --name takei2025_chr15_1mb --out "$W" || echo "ours FAILED"
ours "$C" compare --runs "$W" --reference "$AB/bench/v3/results/arcfish.json"
echo "exit $? end=$(date -Is)"
