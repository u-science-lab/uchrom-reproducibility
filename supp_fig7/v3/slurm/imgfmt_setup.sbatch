#!/bin/bash
# Version 3, family imgfmt, once: (1) the environments of the independent implementations (no U-Chrom):
#   venv-sfimpute  the authors' SnapFISH-IMPUTE (PyPI sfimpute 0.0.2 = github hyuyu104/SnapFISH-IMPUTE @7f7f1a7) with
#                  mpi4py + the MPICH wheel, so that its run_impute.py runs under mpiexec as the README says (T6);
#   venv-jie       the authors' jie (github b2jia/jie @358f365, python-igraph; one line patched for pandas >= 2:
#                  Series.between(inclusive=False) -> inclusive="neither", the same interval) (T7);
#   venv-imgfmt    numpy / pandas / pyarrow / zarr / fsspec / shapely / scikit-learn: outline polygons (F2 build),
#                  the baselines of T7;
# (2) the data folders of T6 and T7 under data_v3/ (v3/data/imgfmt_build.py, from registered sources); the store of
# F2 / F3 is built by imgfmt_outlines.sbatch.
#
#   sbatch supp_fig7/v3/slurm/imgfmt_setup.sbatch      (after syncing v3/ to bench/v3/)
#SBATCH --job-name=ab3-imgfmt-setup
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=$SCRATCH/uchrom/agentbench/logs/%x-%j.out
set -euo pipefail
AB=$SCRATCH/uchrom/agentbench
GH=$HOME/uchrom/agentbench
SIF=$GH/py312full.sif
SRC=$GH/src_imgfmt
DL=$AB/dl/imgfmt
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is)"
mkdir -p "${TMPDIR:-/tmp}/build_tmp" "$AB/data_v3" "$DL" "$SRC"
PIPENV="--env TMPDIR=/tmp --env PIP_CACHE_DIR=/tmp/pip --env HOME=/tmp"

build() {  # name, bash command run inside the image with /opt/venv the new environment, /src the sources
  local name=$1; shift
  rm -rf "$GH/$name"; mkdir -p "$GH/$name"
  local t0=$(date +%s)
  apptainer exec --containall --bind "$GH/$name:/opt/venv" --bind "$SRC:/src" --bind "${TMPDIR:-/tmp}/build_tmp:/tmp" \
    $PIPENV "$SIF" bash -c "python -m venv /opt/venv && /opt/venv/bin/pip install -q --upgrade pip && $*"
  echo "$name built in $(( $(date +%s) - t0 )) s"
  apptainer exec --containall --bind "$GH/$name:/opt/venv" "$SIF" /opt/venv/bin/pip freeze > "$GH/$name.freeze.txt"
}
check() {  # name, python statement
  apptainer exec --containall --cleanenv --bind "$GH/$1:/opt/venv:ro" "$SIF" /opt/venv/bin/python -c "$2" 2>/dev/null
}

# -- (1) environments ----------------------------------------------------------------------------------------
apptainer exec --containall --bind "$SRC:/src" --env HOME=/tmp "$SIF" bash -c '
set -e
cd /src
[ -d SnapFISH-IMPUTE ] || git clone -q https://github.com/hyuyu104/SnapFISH-IMPUTE.git
(cd SnapFISH-IMPUTE && git checkout -q 7f7f1a703aaf6db98147183588a9da6ee1d1c7be && git log -1 --format="SnapFISH-IMPUTE %H %cs")
[ -d jie ] || git clone -q https://github.com/b2jia/jie.git
(cd jie && git checkout -q 358f365 && git log -1 --format="jie %H %cs")
rm -rf jie_patched && cp -r jie jie_patched && rm -rf jie_patched/.git
sed -i "s/inclusive = False/inclusive = \"neither\"/" jie_patched/jie/aligner.py
diff jie/jie/aligner.py jie_patched/jie/aligner.py || true
'
if ! check venv-sfimpute "import sfimpute.impute, mpi4py, skimage"; then
  build venv-sfimpute "/opt/venv/bin/pip install -q sfimpute==0.0.2 mpi4py mpich"
fi
grep -i -E "^(sfimpute|mpi4py|mpich|numpy|pandas|scipy|scikit-image|scikit-learn)==" "$GH/venv-sfimpute.freeze.txt" || true
# the PyPI package is the GitHub commit's sfimpute/ (file by file)
apptainer exec --containall --cleanenv --bind "$GH/venv-sfimpute:/opt/venv:ro" --bind "$SRC:/src:ro" "$SIF" bash -c \
  'for f in __init__ impute parallel preprocess; do cmp /opt/venv/lib/python3.12/site-packages/sfimpute/$f.py /src/SnapFISH-IMPUTE/sfimpute/$f.py && echo "same $f.py"; done'
if ! check venv-jie "import jie.aligner, igraph, sklearn"; then
  build venv-jie "/opt/venv/bin/pip install -q numpy pandas scipy scikit-learn matplotlib seaborn tqdm python-igraph /src/jie_patched"
fi
grep -i -E "^(jie|python-igraph|igraph|numpy|pandas|scipy|scikit-learn)( @|==)" "$GH/venv-jie.freeze.txt" || true
if ! check venv-imgfmt "import shapely, geopandas, zarr, pyarrow, fsspec, aiohttp, sklearn"; then
  build venv-imgfmt "/opt/venv/bin/pip install -q numpy pandas scipy pyarrow zarr fsspec aiohttp requests shapely geopandas scikit-learn"
fi
grep -i -E "^(numpy|pandas|pyarrow|zarr|fsspec|aiohttp|shapely|geopandas|scikit-learn)==" "$GH/venv-imgfmt.freeze.txt" || true

# -- (2) data ------------------------------------------------------------------------------------------------
run() {  # environment, python arguments...
  local v=$1; shift
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=8 --env HOME=/tmp \
    --bind "$SCRATCH/uchrom/data:/shared:ro" --bind "$AB/data_v3:/data_v3" --bind "$DL:/dl" \
    --bind "$AB/bench:/bench:ro" --bind "$GH/$v:/opt/venv:ro" --bind "${TMPDIR:-/tmp}:/tmp" \
    "$SIF" /opt/venv/bin/python -u "$@"
}
run venv-imgfmt /bench/v3/data/imgfmt_build.py sox2 --dl /dl --out /data_v3/imgfmt_sox2
run venv-imgfmt /bench/v3/data/imgfmt_build.py seqfish --shared /shared --out /data_v3/imgfmt_seqfish \
  --chroms chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX
# the store of F2 / F3 (outlines from the 2.47-GB IF zip of Zenodo): imgfmt_outlines.sbatch, a job of its own
du -sh "$AB"/data_v3/imgfmt_*
ls -la "$AB"/data_v3/imgfmt_*
echo "end=$(date -Is)"
