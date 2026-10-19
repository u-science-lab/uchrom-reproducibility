#!/bin/bash
# Version 3, family schic, once: (1) the environments of the authors' code, built as an agent would build them in the
# benchmark's image (python:3.12-bookworm, pip, git, gcc; install times recorded):
#   venv-schicluster  scHiCluster (Zhou et al. 2019): `pip install schicluster==1.3.2` (PyPI's newest), falling back to
#                     --no-deps + the modules its 2019 code imports when the full install fails; the GitHub tag v1.3.4
#                     (a248677, newer than PyPI: the row-normalised compartment score) checked out in src/;
#   venv-nucdyn       NucDynamics (Stevens et al. 2017): github.com/tjs23/nuc_dynamics @ 00eaf0b in src/, its Cython
#                     module built in place (numpy, cython, setuptools);
# (2) the data folders of the family under data_v3/ (v3/truth/schic_build_data.py, from registered sources);
# (3) an inspection of the inputs (v3/truth/schic_explore.py).
#
#   sbatch supp_fig7/v3/sherlock/schic_setup.sbatch      (after syncing v3/ to bench/v3/)
#SBATCH --job-name=ab3-schic-setup
#SBATCH --partition=xiaojie,normal
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -uo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
SIF=$GH/py312full.sif
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is)"
mkdir -p "$L_SCRATCH/build_tmp" "$AB/data_v3"
PIPENV="--env TMPDIR=/tmp --env PIP_CACHE_DIR=/tmp/pip --env HOME=/tmp"
inenv() {  # venv dir, command...
  local v=$1; shift
  apptainer exec --containall --bind "$GH/$v:/opt/venv" --bind "$L_SCRATCH/build_tmp:/tmp" $PIPENV "$SIF" bash -c "$*"
}

# -- (1a) scHiCluster
if [ ! -f "$GH/venv-schicluster/OK" ]; then
  rm -rf "$GH/venv-schicluster"; mkdir -p "$GH/venv-schicluster"
  t0=$(date +%s)
  inenv venv-schicluster "python -m venv /opt/venv && /opt/venv/bin/pip install -q schicluster==1.3.2" \
    > "$GH/venv-schicluster.full_install.log" 2>&1
  rc=$?
  echo "schicluster==1.3.2 full install: exit $rc after $(( $(date +%s) - t0 )) s"; tail -5 "$GH/venv-schicluster.full_install.log"
  if [ $rc -ne 0 ]; then
    rm -rf "$GH/venv-schicluster"; mkdir -p "$GH/venv-schicluster"
    t0=$(date +%s)
    inenv venv-schicluster "python -m venv /opt/venv && /opt/venv/bin/pip install -q --no-deps schicluster==1.3.2 && \
      /opt/venv/bin/pip install -q numpy scipy scikit-learn pandas cooler h5py xarray && \
      /opt/venv/bin/pip install -q torch --index-url https://download.pytorch.org/whl/cpu"
    echo "schicluster==1.3.2 --no-deps + numpy scipy scikit-learn pandas cooler h5py xarray torch(cpu): exit $? after $(( $(date +%s) - t0 )) s"
  fi
  inenv venv-schicluster "cd /opt/venv && git clone -q https://github.com/zhoujt1994/scHiCluster src/scHiCluster && \
    cd src/scHiCluster && git checkout -q a248677867fa4b250bc5b2c5f858871f001e3bb1 && git log -1 --format='scHiCluster %H %cd'"
  inenv venv-schicluster "/opt/venv/bin/python -c 'import schicluster, schicluster.schicluster as s; print(schicluster.__version__, s.__file__)'" \
    && touch "$GH/venv-schicluster/OK"
fi
inenv venv-schicluster "/opt/venv/bin/pip freeze" > "$GH/venv-schicluster.freeze.txt"
grep -i -E "^(schicluster|numpy|scipy|scikit-learn|pandas|cooler|torch)==" "$GH/venv-schicluster.freeze.txt"

# -- (1b) NucDynamics: first as the README says (pip numpy cython, newest), which fails to compile dyn_util.pyx with
#    NumPy 2 / Cython 3 (numpy.math cimport, int division); then with numpy<2 and cython<3 (wheels for Python 3.12)
if [ ! -f "$GH/venv-nucdyn/OK" ]; then
  for pins in "numpy cython setuptools" "numpy<2 cython<3 setuptools"; do
    rm -rf "$GH/venv-nucdyn"; mkdir -p "$GH/venv-nucdyn"
    t0=$(date +%s)
    inenv venv-nucdyn "python -m venv /opt/venv && /opt/venv/bin/pip install -q $(printf "'%s' " $pins) && \
      cd /opt/venv && git clone -q https://github.com/tjs23/nuc_dynamics src/nuc_dynamics && cd src/nuc_dynamics && \
      git checkout -q 00eaf0b57385365c645eee0fda9dd5a2055550a0 && git log -1 --format='nuc_dynamics %H %cd' && \
      /opt/venv/bin/python setup_cython.py build_ext --inplace > /opt/venv/build.log 2>&1; echo build exit \$?; \
      grep -m 5 -E 'Error|error:|cimport' /opt/venv/build.log"
    echo "venv-nucdyn [$pins] (pip, clone, build_ext) in $(( $(date +%s) - t0 )) s"
    if inenv venv-nucdyn "cd /opt/venv/src/nuc_dynamics && /opt/venv/bin/python -c 'import dyn_util; print(dyn_util.__file__)'"; then
      touch "$GH/venv-nucdyn/OK"; break
    fi
  done
fi
inenv venv-nucdyn "/opt/venv/bin/pip freeze" > "$GH/venv-nucdyn.freeze.txt"
cat "$GH/venv-nucdyn.freeze.txt"
[ -n "${ENV_ONLY:-}" ] && exit 0

# -- (2) data folders, (3) inspection (U-Chrom's environment for cooler / pandas)
run() {
  apptainer exec --containall --cleanenv --env OMP_NUM_THREADS=8 \
    --bind "$SCRATCH/uchrom-nd/data:/shared:ro" --bind "$AB/data:/data:ro" --bind "$AB/data_v3:/data_v3" \
    --bind "$AB/bench:/bench:ro" --bind "$GH/venv-uchrom:/opt/venv:ro" --bind "$L_SCRATCH:/tmp" \
    "$SIF" /opt/venv/bin/python -u "$@"
}
run /bench/v3/truth/schic_build_data.py --shared /shared --data /data --out /data_v3
run /bench/v3/truth/schic_explore.py /data
ls -la "$AB/data_v3/stevens2017_cell2" "$AB/data_v3/mm10_cpg" "$AB/data_v3/nagano2017_hap_2i" | head -20
du -sh "$AB/data_v3/nagano2017_hap_2i"
echo "end=$(date -Is)"
