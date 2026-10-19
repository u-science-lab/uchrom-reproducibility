#!/bin/bash
# Version 3, family bulk, once: (1) the environment of the independent implementations of H1-H3 (the authors' /
# Open2C code, no U-Chrom): mustache-hic 1.3.3 and cooltools 0.7.1, which need numpy < 2, pandas < 3, scipy < 1.15
# (benchmarks/hic_structures/README.md), + scikit-learn for H4; (2) the data folders of the family under data_v3/
# (v3/data/bulk_build.py, from registered sources, without U-Chrom's code).
#
#   sbatch supp_fig7/v3/sherlock/bulk_setup.sbatch      (after syncing v3/ to bench/v3/)
#SBATCH --job-name=ab3-bulk-setup
#SBATCH --partition=xiaojie,normal
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -euo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
SIF=$GH/py312full.sif
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is)"
mkdir -p "$L_SCRATCH/build_tmp" "$AB/data_v3" "$AB/dl"
PIPENV="--env TMPDIR=/tmp --env PIP_CACHE_DIR=/tmp/pip --env HOME=/tmp"

# -- (1) venv-hicref: what an agent installs to run Mustache / cooltools dots + pileup on Python 3.12
if ! apptainer exec --containall --cleanenv --bind "$GH/venv-hicref:/opt/venv:ro" "$SIF" /opt/venv/bin/python -c \
     "import mustache.mustache, cooltools, sklearn" 2>/dev/null; then
  rm -rf "$GH/venv-hicref"; mkdir -p "$GH/venv-hicref"
  t0=$(date +%s)
  apptainer exec --containall --bind "$GH/venv-hicref:/opt/venv" --bind "$L_SCRATCH/build_tmp:/tmp" $PIPENV "$SIF" \
    bash -c "python -m venv /opt/venv && /opt/venv/bin/pip install -q 'numpy<2' 'pandas<3' 'scipy<1.15' \
             cooler cooltools==0.7.1 mustache-hic==1.3.3 scikit-learn"
  echo "venv-hicref built in $(( $(date +%s) - t0 )) s"
fi
apptainer exec --containall --bind "$GH/venv-hicref:/opt/venv" "$SIF" /opt/venv/bin/pip freeze > "$GH/venv-hicref.freeze.txt"
grep -i -E "^(numpy|pandas|scipy|cooler|cooltools|mustache-hic|hic-straw|scikit-learn|statsmodels)==" "$GH/venv-hicref.freeze.txt"

# -- (2) the data folders (cooler / pandas / pyarrow of venv-uchrom; bulk_build.py never imports uchrom; the RNA zip
#    is downloaded here, inside the job)
apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=8 \
  --bind "$SCRATCH/uchrom-nd/data:/shared:ro" --bind "$AB/data:/data:ro" --bind "$AB/data_v3:/data_v3" \
  --bind "$AB/dl:/dl" --bind "$AB/bench:/bench:ro" --bind "$GH/venv-uchrom:/opt/venv:ro" --bind "$L_SCRATCH:/tmp" \
  "$SIF" /opt/venv/bin/python -u /bench/v3/data/bulk_build.py --shared /shared --data /data --out /data_v3 --dl /dl
ls -la "$AB/data_v3"/rao2014_imr90_chr21_10kb "$AB/data_v3"/su2020_chr21_tsv "$AB/data_v3"/takei2025_fov0_rna
echo "end=$(date -Is)"
