#!/bin/bash
# The unit tests of the changed packages on the repo copy $AP/repo (uchrom from source, the agents' environment for
# everything else, pytest into the job's scratch).
#
#   TESTS="packages/uchrom/tests/strc packages/uchrom/tests/fea" sbatch supp_table3/arcfish_parity/slurm/tests.sbatch
#SBATCH --job-name=arcfish-tests
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=$SCRATCH/uchrom/arcfish_parity/logs/%x-%j.out
set -uo pipefail
AP=$SCRATCH/uchrom/arcfish_parity
GH=$HOME/uchrom/agentbench
mkdir -p "${TMPDIR:-/tmp}/tmp" "${TMPDIR:-/tmp}/pylib"
apptainer exec --containall --cleanenv --bind "${TMPDIR:-/tmp}/tmp:/tmp" --bind "${TMPDIR:-/tmp}/pylib:/pylib" --bind "$AP:$AP" \
  --bind "$GH/venv-uchrom:/opt/venv:ro" --env TMPDIR=/tmp --env PIP_CACHE_DIR=/tmp/pip --env HOME=/tmp \
  --env OMP_NUM_THREADS=8 "$GH/py312full.sif" bash -c "
/opt/venv/bin/pip install -q --target /pylib pytest pytest-timeout >/dev/null 2>&1
cd $AP/repo/packages/uchrom && PYTHONPATH=/pylib:$AP/repo/packages/uchrom:$AP/repo/packages/chromdata \
  /opt/venv/bin/python -m pytest -q -x --timeout=600 -p no:cacheprovider ${TESTS:-tests/strc tests/fea} 2>&1 | tail -40"
