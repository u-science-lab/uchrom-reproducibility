#!/bin/bash
# Version 3, once: the U-Chrom environment of the agents at the snapshot's commit (venv-uchrom-v3; version 2's
# venv-uchrom, fb43ac4, stays for the record).  The native engines uchrom-recon and uchrom-maps of the same commit are
# built here with Sherlock's Rust (abi3 wheels, glibc 2.17) and added to the snapshot's wheels, as `u-chrom[recon,maps]`
# would install them; the rest as version 2 (torch CPU + the same libraries, then chromdata + u-chrom wheels).
#
#   (local) python benchmarks/fig6/agent_bench/build_snapshot.py --commit <sha> --out snapshot_v3
#           git archive -o uchrom-native-<sha>.tar <sha> packages/uchrom-recon packages/uchrom-maps
#           rsync both to sherlock-dtn: $AB/snapshot_v3/, $AB/src/
#   sbatch benchmarks/fig6/agent_bench/v3/sherlock/setup_v3.sbatch
#SBATCH --job-name=ab3-setup
#SBATCH --partition=xiaojie,normal
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -euo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
SNAP=$AB/snapshot_v3
COMMIT=$(cat "$SNAP/COMMIT")
SIF=$GH/py312full.sif
VENV=venv-uchrom-v3
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} commit=$COMMIT start=$(date -Is)"
mkdir -p "$L_SCRATCH/build_tmp"
PIPENV="--env TMPDIR=/tmp --env PIP_CACHE_DIR=/tmp/pip --env HOME=/tmp"

# -- the native engines of the commit (host build: Sherlock's Rust; abi3 wheels for any Python >= 3.11)
ml python/3.12.1 rust/1.90.0
B=$L_SCRATCH/native_build; rm -rf "$B"; mkdir -p "$B/src" "$B/wheels"
tar -xf "$AB/src/uchrom-native-$COMMIT.tar" -C "$B/src"
python3 -m venv "$B/venv" && "$B/venv/bin/pip" install -q maturin
export CARGO_TARGET_DIR=$B/target CARGO_HOME=$B/cargo
for crate in uchrom-recon uchrom-maps; do
  t0=$SECONDS
  (cd "$B/src/packages/$crate" && "$B/venv/bin/maturin" build --release -o "$B/wheels")
  echo "$crate built in $((SECONDS - t0)) s"
done
rm -f "$SNAP"/wheels/uchrom_recon-*.whl "$SNAP"/wheels/uchrom_maps-*.whl
cp "$B"/wheels/*.whl "$SNAP/wheels/"
ls -la "$SNAP/wheels"

# -- the environment, built inside the image at the path the runs bind it to (/opt/venv)
LIBS="numpy pandas scipy zarr pyarrow fsspec aiohttp anndata cooler h5py scikit-learn matplotlib statsmodels"
rm -rf "$GH/$VENV"; mkdir -p "$GH/$VENV"
t0=$SECONDS
apptainer exec --containall --bind "$GH/$VENV:/opt/venv" --bind "$SNAP:/opt/snapshot:ro" \
  --bind "$L_SCRATCH/build_tmp:/tmp" $PIPENV "$SIF" bash -c "
set -e
python -m venv /opt/venv && /opt/venv/bin/pip install -q --upgrade pip
/opt/venv/bin/pip install -q torch --index-url https://download.pytorch.org/whl/cpu
/opt/venv/bin/pip install -q $LIBS \
  '/opt/snapshot/wheels/chromdata-0.1.0-py3-none-any.whl[contacts,anndata,remote]' \
  /opt/snapshot/wheels/u_chrom-0.3.0.dev0-py3-none-any.whl /opt/snapshot/wheels/uchrom_recon-*.whl \
  /opt/snapshot/wheels/uchrom_maps-*.whl"
echo "$VENV built in $((SECONDS - t0)) s"
apptainer exec --containall --bind "$GH/$VENV:/opt/venv" "$SIF" /opt/venv/bin/pip freeze > "$GH/$VENV.freeze.txt"

# -- check: the packages of the commit, from site-packages; the plain environment still without U-Chrom
apptainer exec --containall --cleanenv --bind "$GH/$VENV:/opt/venv:ro" "$SIF" /opt/venv/bin/python -c "
import inspect, uchrom, chromdata, uchrom_recon, uchrom_maps
from uchrom.emb.hic import schicluster_impute
print('uchrom', uchrom.__version__, uchrom.__file__)
print('chromdata', chromdata.__file__)
print('uchrom_recon', uchrom_recon.__version__, '| uchrom_maps', uchrom_maps.__file__)
print('scHiCluster fix present:', 'n_threads' in inspect.signature(schicluster_impute).parameters)"
apptainer exec --containall --cleanenv --bind "$GH/venv-plain:/opt/venv:ro" "$SIF" /opt/venv/bin/python -c \
  "import importlib.util as u; print('venv-plain', {m: bool(u.find_spec(m)) for m in ('uchrom', 'chromdata', 'zarr', 'cooler')})"
grep -i -E "^(numpy|pandas|scipy|scikit-learn|torch|cooler|zarr|u-chrom|chromdata|uchrom-recon|uchrom-maps)==" "$GH/$VENV.freeze.txt"
echo "end=$(date -Is)"
