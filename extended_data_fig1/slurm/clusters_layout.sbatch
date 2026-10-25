#!/bin/bash
# The locus index's chunk size and codec on real SPRITE clusters (g_clusters_layout.py; U-Chrom PR #153, format
# 2.4): Quinodoz et al. 2018 mouse ES cells, GEO GSE114242 (fetch quinodoz2018_sprite).  Stores on node-local disk,
# results to $W/runs/clusters/<job>/ (JSON + log archived in long-term storage).
#   sbatch --export=ALL,REPO=$SCRATCH/uchrom/repo-mw,REPRO=$SCRATCH/uchrom/repro extended_data_fig1/slurm/clusters_layout.sbatch
#SBATCH --job-name=clusters-layout
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --output=%x-%j.out
source "$REPRO/extended_data_fig1/slurm/_env.sh"
ARCHIVE=$UCHROM_DATA/quinodoz2018_sprite/GSE114242_mouse_combined_clusters.tar.gz
RUN=$W/runs/clusters/${SLURM_JOB_ID:-local}
WORK=${TMPDIR:-/tmp}/layout
mkdir -p "$RUN" "$WORK"
cd "$REPRO"
$PY -u extended_data_fig1/g_clusters_layout.py --archive "$ARCHIVE" --work "$WORK" --out "$RUN/clusters_layout.json"
ARCH=$UCHROM_ARCHIVE/results/clusters_sprite/$(date +%Y-%m-%d)_${UCHROM_REVISION}_layout
mkdir -p "$ARCH" && cp "$RUN"/*.json "$ARCH"/ && cp "$SLURM_SUBMIT_DIR/clusters-layout-${SLURM_JOB_ID}.out" "$ARCH"/ 2>/dev/null || true
echo "archived to $ARCH"
