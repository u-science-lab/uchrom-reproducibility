#!/bin/bash
# Archive the large outputs of a benchmark run from $SCRATCH (purged after 90 days) to the archive, one tar per part
# (the archive counts files: bundle them), with a manifest.  The headline numbers are committed in the benchmark's results/
# folder in uchrom-reproducibility; this keeps what they were computed from (call tables, maps, models, logs).
#
#   sbatch --export=ALL,RUN=sc,DEST=sc_structures/2026-10-08_ad141a7 supp_table3/sc_structures/slurm/archive.sbatch
#
# EXCLUDE: one more tar --exclude pattern (files the archive keeps another copy of, e.g. an unbalanced map next
# to its balanced copy).
#SBATCH --job-name=archive
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out
set -euo pipefail
SRC=$SCRATCH/uchrom/runs/${RUN:?}
OUT=$UCHROM_ARCHIVE/results/${DEST:?}
mkdir -p "$OUT" "$UCHROM_ARCHIVE/results/hic_structures"
cd "$SRC"
for part in */; do
  part=${part%/}
  # SnapHiC's own run (the reference implementation): its final loops and log, not its 19 GB of per-cell files
  tar -cf "$OUT/$part.tar" --exclude="reference/*/input" --exclude="reference/*/snaphic_out/binned" \
      --exclude="reference/*/snaphic_out/hic" --exclude="reference/*/snaphic_out/interactions" \
      --exclude="reference/*/snaphic_out/rwr" ${EXCLUDE:+--exclude="$EXCLUDE"} "$part"
done
{
  echo "# $DEST — archived $(date -u +%FT%TZ) from \$SCRATCH/uchrom/runs/$RUN (job $SLURM_JOB_ID)"
  echo "# headline numbers: the results/ of the $(dirname "$DEST") benchmark in uchrom-reproducibility"
  for t in "$OUT"/*.tar; do echo "$(basename "$t")  $(du -h "$t" | cut -f1)  $(tar -tf "$t" | wc -l) files"; done
} > "$OUT/MANIFEST.txt"
cat "$OUT/MANIFEST.txt"
