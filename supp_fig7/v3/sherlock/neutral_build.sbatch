#!/bin/bash
# Version 3: the data of T1, P1 and P2 in community formats next to the shared stores (v3/data/neutral_build.py,
# v3/data/neutral.md): Takei 2025 field of view 0 as FOF-CT core + cell tables, the melanoma section r1 as a per-spot
# .scool with trans contacts + a spot table; each checked against its store and against the references of P1 / P2.
#
#   [WHAT=takei|melanoma] sbatch --export=ALL supp_fig7/v3/sherlock/neutral_build.sbatch
#   (after syncing v3/ to bench/v3/; both by default)
#SBATCH --job-name=ab3-neutral
#SBATCH --partition=xiaojie,normal
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -uo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
SIF=$GH/py312full.sif
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is)"
run() {  # python arguments...
  local t0=$(date +%s)
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=4 --env HOME=/tmp \
    --bind "$AB/data:/data:ro" --bind "$AB/data_v3:/data_v3" --bind "$AB/bench:/bench:ro" \
    --bind "$GH/venv-uchrom-v3:/opt/venv:ro" --bind "$L_SCRATCH:/tmp" \
    "$SIF" /opt/venv/bin/python -u "$@"
  local rc=$?
  echo "== exit $rc after $(( $(date +%s) - t0 )) s: $*"
  return $rc
}
B=/bench/v3/data/neutral_build.py
R=/bench/v3/results/bulk.json
WHAT=${WHAT:-takei melanoma}
[[ $WHAT == *takei* ]] && run $B takei --store /data/takei2025_fov0/takei2025_fov0.chromdata.zarr \
  --rna /data_v3/takei2025_fov0_rna/takei2025_fov0_rna_counts.tsv --results $R --out /data_v3/takei2025_fov0_fofct \
  --report /data_v3/_build_info/neutral_build_takei.json
[[ $WHAT == *melanoma* ]] && run $B melanoma --store /data/hicrna_human_melanoma/hicrna_human_melanoma.chromdata.zarr \
  --section human_melanoma_r1 --results $R --out /data_v3/hicrna_melanoma_r1 \
  --report /data_v3/_build_info/neutral_build_melanoma.json
du -sh "$AB"/data_v3/takei2025_fov0_fofct "$AB"/data_v3/hicrna_melanoma_r1
ls -la "$AB"/data_v3/takei2025_fov0_fofct "$AB"/data_v3/hicrna_melanoma_r1
echo "end=$(date -Is)"
