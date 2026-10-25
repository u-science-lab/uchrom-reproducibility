#!/bin/bash
# Version 3, family arcfish: T1-T5 again with the agents' environment of the snapshot (venv-uchrom-v3: U-Chrom from
# site-packages, main b6d8873) and the authors' package (venv-arcfish, a21d889), on the inputs of main's
# benchmarks/arcfish_parity/inputs.py (Bintu 2018 insulation at prominence 0.01, as T5 asks), all five datasets; then
# the comparison (the parity file of benchmarks/arcfish_parity) and v3/truth/arcfish.py.
#
#   (local) git archive <sha> benchmarks/arcfish_parity benchmarks/_record.py -> $SCRATCH/uchrom/repo-main/ (+ REVISION)
#   sbatch [--dependency=afterok:<setup_v3>] supp_fig7/v3/slurm/arcfish_truth.sbatch
#SBATCH --job-name=ab3-arcfish-truth
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=$SCRATCH/uchrom/agentbench/logs/%x-%j.out
set -uo pipefail
AB=$SCRATCH/uchrom/agentbench
AP=$SCRATCH/uchrom/arcfish_parity
GH=$HOME/uchrom/agentbench
SIF=$GH/py312full.sif
D=$SCRATCH/uchrom/data
REPO=$SCRATCH/uchrom/repo-main
REV=$(cat "$REPO/REVISION")
IN=$AP/inputs_$REV; OUT=$AP/runs/${SLURM_JOB_ID}
mkdir -p "$IN" "$OUT" "${TMPDIR:-/tmp}/tmp"
echo "host=$(hostname) job=$SLURM_JOB_ID revision=$REV snapshot=$(cat $AB/snapshot_v3/COMMIT) start=$(date -Is)"
COMMON=(--containall --cleanenv --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 --env MKL_NUM_THREADS=8
        --env UCHROM_REVISION="$REV" --env SLURM_JOB_ID="$SLURM_JOB_ID" --env UCHROM_DATA="$D"
        --bind "${TMPDIR:-/tmp}/tmp:/tmp" --bind "$D:$D" --bind "$AP:$AP" --bind "$REPO:$REPO:ro" --bind "$AB:$AB")
ours() { apptainer exec "${COMMON[@]}" --bind "$GH/venv-uchrom-v3:/opt/venv:ro" --env PYTHONPATH="$REPO/benchmarks" \
  "$SIF" /opt/venv/bin/python -u "$@"; }
theirs() { apptainer exec "${COMMON[@]}" --bind "$GH/venv-arcfish:/opt/venv:ro" --env PYTHONPATH=/opt/venv/src \
  "$SIF" /opt/venv/bin/python -u "$@"; }
B=$REPO/benchmarks/arcfish_parity
ours -c "import uchrom, sys; print('uchrom', uchrom.__file__); assert '/opt/venv/' in uchrom.__file__"
[ -f "$IN/manifest.json" ] || ours "$B/inputs.py" --data "$D" --takei2025 "$AB/data/takei2025_fov0/takei2025_fov0.chromdata.zarr" --out "$IN"
for name in ${NAMES:-takei2021_chr6 su2020_chr2p takei2025_chr15_1mb bintu2018 su2020_chr21}; do
  t0=$SECONDS; theirs "$B/run_upstream.py" --manifest "$IN/manifest.json" --name "$name" --out "$OUT" || echo "upstream $name FAILED"
  echo "upstream $name $((SECONDS - t0)) s"
  t0=$SECONDS; ours "$B/run_ours.py" --manifest "$IN/manifest.json" --name "$name" --out "$OUT" || echo "ours $name FAILED"
  echo "ours $name $((SECONDS - t0)) s"
done
ours "$B/compare.py" --runs "$OUT" --out "$OUT/parity.json"
ours "$AB/bench/v3/truth/arcfish.py" --runs "$OUT" --out "$OUT/arcfish_v3.json"
echo "end=$(date -Is)"
