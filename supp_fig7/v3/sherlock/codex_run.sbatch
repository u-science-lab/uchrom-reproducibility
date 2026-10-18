#!/bin/bash
# One bounded batch of Codex sessions of benchmark version 3 (run_codex.py): the questions x conditions x models x
# repeats given, every session to completion (or its 1,800-s limit), then the job ends.  Sessions already done are
# skipped, so the same command started again later fills in the ones a usage limit stopped -- started by a person,
# never by this script (Sherlock's rules on unattended agent sessions, /etc/agents/policy.md).
#
#   sbatch --export=ALL,RUN=pilot_v3,MODELS=gpt-6-luna,REPS=1,CONDITIONS=none:skills[,IDS=T1:R3,JOBS=8,THREADS=4,MEM_GB=24] \
#       benchmarks/fig6/agent_bench/v3/sherlock/codex_run.sbatch
# (lists separated by ":" -- sbatch --export splits on commas)
#SBATCH --job-name=ab3-codex
#SBATCH --partition=xiaojie
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
#SBATCH --no-requeue
# --no-requeue: on the owners partition a preempted batch ends instead of starting its agent sessions again by
# itself; a person starts the same command again (done sessions are skipped).
set -uo pipefail
AB=$SCRATCH/uchrom-nd/agentbench
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
ml python/3.12.1
PY=$GROUP_HOME/$USER/uchrom-nd/venv-consensus/bin/python             # the runner (pyyaml), on the host
CODEX=/share/software/user/open/codex/0.156.0/lib/node_modules/@openai/codex/node_modules/@openai/codex-linux-x64/vendor/x86_64-unknown-linux-musl
L() { echo "${1:-}" | tr ":" ","; }
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} run=$RUN models=$(L "$MODELS") reps=${REPS:-3} ids=$(L "${IDS:-all}") \
conditions=$(L "${CONDITIONS:-none:bare:docs:skills}") snapshot=$(cat "$AB/snapshot_v3/COMMIT") start=$(date -Is)"
"$PY" -u "$AB/bench/run_codex.py" --snapshot "$AB/snapshot_v3" --data-root "$AB/data" --data-root-v3 "$AB/data_v3" \
  --image "$GH/py312full.sif" --venv-uchrom "$GH/venv-uchrom-v3" --venv-plain "$GH/venv-plain" --codex-dir "$CODEX" \
  --questions "$AB/bench/questions_v3.yaml" --out "$AB/runs" --run "$RUN" --models "$(L "$MODELS")" \
  --reps "${REPS:-3}" --conditions "$(L "${CONDITIONS:-none:bare:docs:skills}")" ${IDS:+--ids "$(L "$IDS")"} \
  --jobs "${JOBS:-8}" --threads "${THREADS:-4}" --mem-gb "${MEM_GB:-24}" --timeout 1800
echo "end=$(date -Is)"
