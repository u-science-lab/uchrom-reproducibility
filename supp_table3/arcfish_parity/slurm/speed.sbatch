#!/bin/bash
# The ArcFISH TAD caller before (the agents' environment, fb43ac4) and after the parity fix (the repo copy), on the same
# node: Su 2020 chr21, all 7,591 traces, in memory and streaming (backed store), CPU.
#SBATCH --job-name=arcfish-speed
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=$SCRATCH/uchrom/arcfish_parity/logs/%x-%j.out
AP=$SCRATCH/uchrom/arcfish_parity
GH=$HOME/uchrom/agentbench
D=$SCRATCH/uchrom/data
cat > "${TMPDIR:-/tmp}/speed.py" <<'PY'
import sys, time, warnings
warnings.simplefilter("ignore")
import uchrom
from chromdata import ChromData
from uchrom.strc.tad import TADCallerParams, call_tads_by_pval
path = sys.argv[1]
w = 1e5 if "parity" in uchrom.__file__ else 2e5          # +/-100 kb in both conventions
for backed in (False, True):
    cd = ChromData.read(path, backed=backed)
    t0 = time.time()
    t = call_tads_by_pval(cd, chrom="chr21", params=TADCallerParams(window_bp=w), device="cpu", key_added=None)
    print(("new" if "parity" in uchrom.__file__ else "old"), "backed" if backed else "memory", f"{time.time() - t0:.1f} s", len(t), flush=True)
PY
for which in old new; do
  if [ $which = old ]; then PP=""; else PP="$AP/repo/packages/uchrom:$AP/repo/packages/chromdata"; fi
  apptainer exec --containall --cleanenv --bind "${TMPDIR:-/tmp}:${TMPDIR:-/tmp}" --bind "$AP:$AP" --bind "$D:$D" \
    --bind "$GH/venv-uchrom:/opt/venv:ro" --env OMP_NUM_THREADS=8 --env PYTHONPATH="$PP" \
    "$GH/py312full.sif" /opt/venv/bin/python -u "${TMPDIR:-/tmp}/speed.py" "$D/su2020_chr21/su2020_chr21.chromdata.zarr"
done
