#!/bin/bash
# Multi-way contacts in the store (U-Chrom PR #153, format 2.4) on real SPRITE clusters (g_clusters_sprite.py):
# Quinodoz et al. 2018 mouse ES cells, GEO GSE114242 (python -m uchrom.datasets fetch quinodoz2018_sprite, or the
# archive in $UCHROM_DATA/quinodoz2018_sprite).  A short run on the first 2 M clusters, then the whole member;
# store on node-local disk, results to $W/runs/clusters/<job>/ (JSON + log archived in long-term storage).
#   sbatch --export=ALL,REPO=$SCRATCH/uchrom/repo-mw,REPRO=$SCRATCH/uchrom/repro extended_data_fig1/slurm/clusters.sbatch
#SBATCH --job-name=clusters-sprite
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --output=%x-%j.out
source "$REPRO/extended_data_fig1/slurm/_env.sh"
ARCHIVE=$UCHROM_DATA/quinodoz2018_sprite/GSE114242_mouse_combined_clusters.tar.gz
RUN=$W/runs/clusters/${SLURM_JOB_ID:-local}
WORK=${TMPDIR:-/tmp}/clusters
mkdir -p "$RUN"
cd "$REPRO"
$PY -u extended_data_fig1/g_clusters_sprite.py --archive "$ARCHIVE" --limit 2000000 --work "$WORK/limit" \
    --out "$RUN/clusters_sprite_2M.json"
rm -rf "$WORK/limit"
/usr/bin/time -v $PY -u extended_data_fig1/g_clusters_sprite.py --archive "$ARCHIVE" --work "$WORK/full" \
    --out "$RUN/clusters_sprite.json"
ARCH=$UCHROM_ARCHIVE/results/clusters_sprite/$(date +%Y-%m-%d)_${UCHROM_REVISION}
mkdir -p "$ARCH" && cp "$RUN"/*.json "$ARCH"/ && cp "$SLURM_SUBMIT_DIR/clusters-sprite-${SLURM_JOB_ID}.out" "$ARCH"/ 2>/dev/null || true
echo "archived to $ARCH"
