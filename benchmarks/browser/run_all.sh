#!/usr/bin/env bash
# Fig 6c benchmarks.  Run on an idle machine: each script waits (up to 30 min) for load < 2.0, then refuses;
# close browsers, games, VMs.
#   bash benchmarks/browser/run_all.sh <data directory>
# Results go to benchmarks/browser/results (see benchmarks/_paper.py).
set -euo pipefail
cd "$(dirname "$0")/../.."
DATA=${1:-${UCHROM_DATA:-$HOME/.cache/uchrom}}
PY=${PY:-.venv/bin/python}
EXTRA=${EXTRA:-}   # e.g. EXTRA=--force to skip the idle-machine check (the load is still recorded)
D=benchmarks/browser/_data
[ -e $D/takei2025_fov0_x30.chromdata.zarr ] || \
    $PY benchmarks/browser/make_scaled.py --src "$DATA/takei2025_fov0/takei2025_fov0.chromdata.zarr" --factors 1 3 10 30
$PY -u benchmarks/browser/bench_frontend.py $D/takei2025_fov0_x{1,3,10,30}.chromdata.zarr --reps 5 $EXTRA
cd benchmarks/browser
../../$PY -u bench_server.py $EXTRA --imaging "../../$DATA/takei2025_fov0/takei2025_fov0.chromdata.zarr" \
    --scaled _data/takei2025_fov0_x10.chromdata.zarr _data/takei2025_fov0_x30.chromdata.zarr \
    --hic "../../$DATA/schicar_mop/schicar_mop.chromdata.zarr"
