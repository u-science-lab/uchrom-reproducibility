#!/bin/bash
# The ArcFISH parity check: inputs (U-Chrom), every dataset through the authors' package and through U-Chrom, then
# the comparison.  U-Chrom of the repo copy $AP/repo (its packages on PYTHONPATH ahead of the agents' environment); the
# benchmark from a copy of uchrom-reproducibility ($REPRO).
#
#   sbatch supp_table3/arcfish_parity/slurm/run.sbatch        (after slurm/setup.sbatch)
#SBATCH --job-name=arcfish-parity
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=$SCRATCH/uchrom/arcfish_parity/logs/%x-%j.out
set -euo pipefail
AP=$SCRATCH/uchrom/arcfish_parity
REPRO=${REPRO:-$SCRATCH/uchrom/repro}
GH=$HOME/uchrom/agentbench
SIF=$GH/py312full.sif
D=$SCRATCH/uchrom/data
IN=$AP/inputs; OUT=$AP/runs/${SLURM_JOB_ID}
mkdir -p "$IN" "$OUT" "${TMPDIR:-/tmp}/tmp"
export UCHROM_REVISION=$(cat "$AP/repo/REVISION")
COMMON=(--containall --cleanenv --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 --env MKL_NUM_THREADS=8
        --bind "${TMPDIR:-/tmp}/tmp:/tmp" --bind "$D:$D" --bind "$AP:$AP" --bind "$REPRO:$REPRO" --bind "$SCRATCH/uchrom/agentbench/data:/abdata:ro")
ours() { apptainer exec "${COMMON[@]}" --bind "$GH/venv-uchrom:/opt/venv:ro" --env UCHROM_REVISION="$UCHROM_REVISION" \
  --env PYTHONPATH="$AP/repo/packages/uchrom:$AP/repo/packages/chromdata:$REPRO/shared" "$SIF" /opt/venv/bin/python -u "$@"; }
theirs() { apptainer exec "${COMMON[@]}" --bind "$GH/venv-arcfish:/opt/venv:ro" --env PYTHONPATH=/opt/venv/src \
  "$SIF" /opt/venv/bin/python -u "$@"; }
B=$REPRO/supp_table3/arcfish_parity
[ -f "$IN/manifest.json" ] || ours "$B/inputs.py" --data "$D" --takei2025 /abdata/takei2025_fov0/takei2025_fov0.chromdata.zarr --out "$IN"
for name in ${NAMES:-takei2021_chr6 su2020_chr2p takei2025_chr15_1mb bintu2018 su2020_chr21}; do
  theirs "$B/run_upstream.py" --manifest "$IN/manifest.json" --name "$name" --out "$OUT" || echo "upstream $name FAILED"
  ours "$B/run_ours.py" --manifest "$IN/manifest.json" --name "$name" --out "$OUT" || echo "ours $name FAILED"
done
ours "$B/compare.py" --runs "$OUT" --out "$OUT/parity.json"
