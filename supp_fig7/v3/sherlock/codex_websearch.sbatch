#!/bin/bash
# Does the Codex config turn web search off?  Two candidate settings, one session each (gpt-6-luna): the events show
# whether a web_search item is made.
#SBATCH --job-name=ab3-codex-web
#SBATCH --partition=xiaojie,normal
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -uo pipefail
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
SIF=$GH/py312full.sif
V=/share/software/user/open/codex/0.156.0/lib/node_modules/@openai/codex/node_modules/@openai/codex-linux-x64/vendor/x86_64-unknown-linux-musl
run() {  # name, config.toml body
  local T=$L_SCRATCH/web_$1; rm -rf "$T"; mkdir -p "$T/home/.codex" "$T/work" "$T/tmp"
  cp "$HOME/.codex/auth.json" "$T/home/.codex/auth.json"; chmod 600 "$T/home/.codex/auth.json"
  printf '%s\n' "$2" > "$T/home/.codex/config.toml"
  apptainer exec --containall --cleanenv --home "$T/home:/home/agent" --pwd /work --bind "$T/work:/work" \
    --bind "$T/tmp:/tmp" --bind "$V:/opt/codex:ro" --bind /etc/codex:/etc/codex:ro --bind /etc/agents:/etc/agents:ro \
    --env HOME=/home/agent --env CODEX_HOME=/home/agent/.codex --env TMPDIR=/tmp --env PATH=/usr/local/bin:/usr/bin:/bin \
    "$SIF" /opt/codex/bin/codex exec --json --skip-git-repo-check -C /work -m gpt-6-luna \
    "Use your web search tool (not the shell, no curl) to find the latest cooltools release on PyPI and report it. If you have no web search tool, say NO_WEB_TOOL and stop." \
    < /dev/null > "$T/events.jsonl" 2> "$T/stderr.txt"
  echo "== $1 rc=$?"
  python3 - "$T/events.jsonl" <<'PY'
import json, sys, collections
c = collections.Counter(); last = None
for line in open(sys.argv[1]):
    try: j = json.loads(line)
    except ValueError: continue
    it = j.get("item") or {}
    if it: c[str(it.get("type"))] += 1
    if it.get("type") == "agent_message": last = it.get("text")
    if it.get("type") == "mcp_tool_call" and j.get("type") == "item.completed":
        print("MCP", {k: it.get(k) for k in ("server", "tool", "status")}, str(it.get("arguments"))[:120])
print(dict(c)); print("LAST:", (last or "")[:300])
PY
}
BASE='approval_policy = "never"
sandbox_mode = "danger-full-access"
check_for_update_on_startup = false'
run all_off "$BASE
web_search = \"disabled\"
[features]
memories = false
multi_agent = false
apps = false
browser_use = false
browser_use_external = false
in_app_browser = false
computer_use = false
image_generation = false
plugins = false
remote_plugin = false"
