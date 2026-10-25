#!/bin/bash
# Archive agent benchmark version 3 from $SCRATCH (purged after 90 days) to the archive: the code and references as run
# (bench/v3 + the merged questions / answers / grader), the main run and pilots (session transcripts, answers, grades)
# and the inputs built for v3 (data_v3), one tar per part, with a manifest.
#
#   sbatch bench/v3/slurm/archive_v3.sbatch
#SBATCH --job-name=ab3-archive
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=$SCRATCH/uchrom/agentbench/logs/%x-%j.out
set -euo pipefail
AB=$SCRATCH/uchrom/agentbench
OUT=$UCHROM_ARCHIVE/results/agent_bench/2026-10-14_v3
mkdir -p "$OUT"
cd "$AB"
tar -cf "$OUT/bench_v3.tar" bench/v3 bench/questions_v3.yaml bench/answers_v3.json bench/merge_v3.py bench/grade_v3.sbatch bench/build_snapshot.py
tar -cf "$OUT/runs_main_v3.tar" runs/main_v3
tar -cf "$OUT/runs_pilots_v3.tar" runs/pilot_v3b_noclass runs/pilot_v3c runs/pilot_v3c_http_before_fix runs/pilot_v3d_neutral runs/pilot_v3_nopath
tar -cf "$OUT/data_v3.tar" data_v3
{
  echo "# agent benchmark v3 -- archived $(date -u +%FT%TZ) from \$SCRATCH/uchrom/agentbench (job $SLURM_JOB_ID)"
  echo "# code and headline results: uchrom-reproducibility supp_fig7/v3/, supp_fig7/results/agent_bench_v3.json"
  for t in "$OUT"/*.tar; do echo "$(basename "$t")  $(du -h "$t" | cut -f1)  $(tar -tf "$t" | wc -l) files"; done
} > "$OUT/MANIFEST.txt"
cat "$OUT/MANIFEST.txt"
