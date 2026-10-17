#!/usr/bin/env bash
# The browser benchmarks (Supplementary Fig. 6a).  Run on an idle machine: each script waits (up to 30 min) for load < 2.0, then refuses;
# close browsers, games, VMs.
#   bash supp_fig6/browser/run_all.sh <data directory>
# Results go to supp_fig6/browser/results (see shared/paper_paths.py).
set -euo pipefail
cd "$(dirname "$0")/../.."
DATA=${1:-${UCHROM_DATA:-$HOME/.cache/uchrom}}
PY=${PY:-.venv/bin/python}
EXTRA=${EXTRA:-}   # e.g. EXTRA=--force to skip the idle-machine check (the load is still recorded)
D=supp_fig6/browser/_data
[ -e $D/takei2025_fov0_x30.chromdata.zarr ] || \
    $PY supp_fig6/browser/make_scaled.py --src "$DATA/takei2025_fov0/takei2025_fov0.chromdata.zarr" --factors 1 3 10 30
$PY -u supp_fig6/browser/bench_frontend.py $D/takei2025_fov0_x{1,3,10,30}.chromdata.zarr --reps 5 $EXTRA
cd supp_fig6/browser
../../$PY -u bench_server.py $EXTRA --imaging "../../$DATA/takei2025_fov0/takei2025_fov0.chromdata.zarr" \
    --scaled _data/takei2025_fov0_x10.chromdata.zarr _data/takei2025_fov0_x30.chromdata.zarr \
    --hic "../../$DATA/schicar_mop/schicar_mop.chromdata.zarr"
