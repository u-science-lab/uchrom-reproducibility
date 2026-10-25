#!/bin/bash
# Grade one run of benchmark version 3 (v3/grade.py with answers_v3.json / questions_v3.yaml; the judge is a
# `claude -p` session without tools, memory, user settings or MCP servers; numbers tested in code), then print the
# mean score per model x condition x kind and per question x condition.  Grades are written as they come, so a
# preempted job loses nothing: the same command grades the rest.
#
#   sbatch --export=ALL,RUN=pilot_v3c[,JUDGE=claude-opus-5-5] supp_fig7/v3/slurm/grade_v3.sbatch
#SBATCH --job-name=ab3-grade
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=$SCRATCH/uchrom/agentbench/logs/%x-%j.out
set -uo pipefail
AB=$SCRATCH/uchrom/agentbench
ml python/3.12.1 claude-code
PY=$HOME/uchrom/venv-consensus/bin/python
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} run=$RUN judge=${JUDGE:-claude-opus-5-5} claude=$(claude --version 2>/dev/null)"
cd "$AB/bench"
"$PY" -u v3/grade.py --runs "$AB/runs/$RUN" --answers answers_v3.json --questions questions_v3.yaml \
  --judge-model "${JUDGE:-claude-opus-5-5}" --jobs "${JOBS:-8}"
"$PY" - "$AB/runs/$RUN/grades.json" <<'PY'
import json, sys, collections
g = json.load(open(sys.argv[1]))
t, q = collections.defaultdict(list), collections.defaultdict(dict)
for v in g.values():
    t[(v["model_requested"], v["condition"], v.get("kind"))].append(v["score"])
    t[(v["model_requested"], v["condition"], "all")].append(v["score"])
    q[v["id"]][v["condition"]] = round(v["score"], 2)
for k in sorted(t, key=str):
    print(f"{k[0]:12s} {k[1]:7s} {str(k[2]):7s} n={len(t[k]):4d} mean={sum(t[k]) / len(t[k]):.3f}")
for i in sorted(q):
    print(f"{i:4s} " + "  ".join(f"{c}={s}" for c, s in sorted(q[i].items())))
PY
