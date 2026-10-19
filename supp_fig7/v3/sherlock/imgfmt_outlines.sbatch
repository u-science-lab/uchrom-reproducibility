#!/bin/bash
# Version 3, family imgfmt: the store of F2 / F3 (Takei et al. 2021 replicate 1, field of view 4, with nuclear outlines
# and two separation_profile results).  Zenodo serves the IF zip (2.47 GB) at ~0.3 MB/s and times out on range
# requests far into it, so the whole file is downloaded (resumable, md5-checked by the build), then the outlines and
# the store are built (v3/data/imgfmt_build.py outlines / store).  Needs the environments of imgfmt_setup.sbatch.
#
#   sbatch supp_fig7/v3/sherlock/imgfmt_outlines.sbatch
#SBATCH --job-name=ab3-imgfmt-outlines
#SBATCH --partition=xiaojie,normal
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -euo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
SIF=$GH/py312full.sif
DL=$AB/dl/imgfmt
mkdir -p "$DL" "$AB/data_v3"
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is)"
Z=$DL/IF-binned1x1x1-E14-replicate1-Pos4.zip
URL="https://zenodo.org/records/3735329/files/IF-binned1x1x1-E14-replicate1-Pos4.zip?download=1"
t0=$(date +%s)
for attempt in $(seq 1 30); do
  size=$( [ -f "$Z" ] && stat -c %s "$Z" || echo 0 )
  [ "$size" = 2471856846 ] && break
  echo "attempt $attempt at $size bytes ($(date -Is))"
  curl -sSL -C - --retry 5 --retry-delay 30 --speed-time 300 --speed-limit 1000 -o "$Z" "$URL" || sleep 60
done
echo "IF zip: $(stat -c %s "$Z") bytes after $(( $(date +%s) - t0 )) s"
run() {  # environment, python arguments...
  local v=$1; shift
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=4 --env HOME=/tmp \
    --bind "$SCRATCH/uchrom-nd/data:/shared:ro" --bind "$AB/data_v3:/data_v3" --bind "$DL:/dl" \
    --bind "$AB/bench:/bench:ro" --bind "$GH/$v:/opt/venv:ro" --bind "$L_SCRATCH:/tmp" \
    "$SIF" /opt/venv/bin/python -u "$@"
}
run venv-imgfmt /bench/v3/data/imgfmt_build.py outlines --ifzip /dl/IF-binned1x1x1-E14-replicate1-Pos4.zip \
  --shared /shared --fov 4 --out /dl/takei2021_fov4_outlines.json
run venv-uchrom /bench/v3/data/imgfmt_build.py store --shared /shared --outlines /dl/takei2021_fov4_outlines.json \
  --fov 4 --out /data_v3/imgfmt_takei2021_fov4/takei2021_fov4.chromdata.zarr
du -sh "$AB"/data_v3/imgfmt_takei2021_fov4
echo "end=$(date -Is)"
