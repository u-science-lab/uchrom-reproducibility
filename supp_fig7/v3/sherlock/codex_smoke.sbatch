#!/bin/bash
# Codex CLI inside the benchmark's container, once: which skill folders it discovers (CODEX_HOME/skills,
# <work>/.codex/skills, <work>/.agents/skills), its JSON events, that a shell command runs, no memory / web search.
#SBATCH --job-name=ab3-codex-smoke
#SBATCH --partition=xiaojie,normal
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -uo pipefail
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
SIF=$GH/py312full.sif
V=/share/software/user/open/codex/0.156.0/lib/node_modules/@openai/codex/node_modules/@openai/codex-linux-x64/vendor/x86_64-unknown-linux-musl
T=$L_SCRATCH/codex_smoke; rm -rf "$T"; mkdir -p "$T/home/.codex" "$T/work" "$T/tmp"
cp "$HOME/.codex/auth.json" "$T/home/.codex/auth.json"; chmod 600 "$T/home/.codex/auth.json"
cat > "$T/home/.codex/config.toml" <<'TOML'
approval_policy = "never"
sandbox_mode = "danger-full-access"
web_search = "disabled"
check_for_update_on_startup = false
[features]
memories = false
multi_agent = false
TOML
mk() { mkdir -p "$1/$2"; printf -- "---\nname: %s\ndescription: A probe skill named %s (benchmark smoke test). Use it when asked about probes.\n---\nSay PROBE-%s.\n" "$2" "$2" "$2" > "$1/$2/SKILL.md"; }
mk "$T/home/.codex/skills" probe-home
mk "$T/work/.codex/skills" probe-project
mk "$T/work/.agents/skills" probe-agents
apptainer exec --containall --cleanenv --home "$T/home:/home/agent" --pwd /work \
  --bind "$T/work:/work" --bind "$T/tmp:/tmp" --bind "$V:/opt/codex:ro" \
  --bind /etc/codex:/etc/codex:ro --bind /etc/agents:/etc/agents:ro \
  --env HOME=/home/agent --env CODEX_HOME=/home/agent/.codex --env TMPDIR=/tmp \
  --env PATH=/usr/local/bin:/usr/bin:/bin "$SIF" \
  /opt/codex/bin/codex exec --json --skip-git-repo-check -C /work -m gpt-6-luna \
  "List the names of every skill available to you (exact names). Then run the shell command: python -c 'print(6*7)' and report its output. Finally say whether you can search the web." \
  > "$T/events.jsonl" 2> "$T/stderr.txt"
echo "rc=$?"
python3 - "$T/events.jsonl" <<'PY'
import json, sys, collections
c = collections.Counter(); last = None
for line in open(sys.argv[1]):
    try: j = json.loads(line)
    except ValueError: continue
    t = j.get("type"); c[t] += 1
    it = j.get("item") or {}
    if it: c["item:" + str(it.get("type"))] += 1
    if t in ("turn.completed", "turn.failed", "error", "thread.started"): print(json.dumps(j)[:600])
    if it.get("type") == "agent_message": last = it.get("text")
    if it.get("type") == "command_execution" and t == "item.completed": print("CMD", json.dumps(it)[:400])
print(dict(c)); print("LAST:", last)
PY
tail -5 "$T/stderr.txt"; ls -la "$T/home/.codex"
