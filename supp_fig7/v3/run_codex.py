"""Version 3 of the agent benchmark: each question x condition x model x repeat is one `codex exec` session (Codex
CLI, `module load codex`) in its own Apptainer container and its own Slurm job step, inside a Slurm allocation.

Conditions (what the agent gets besides the data and a Python with the scientific stack):

  none    no U-Chrom at all: numpy, pandas, scipy, zarr, pyarrow, anndata, cooler ... but not chromdata / uchrom
  bare    the `uchrom` package installed (site-packages), no documentation, no skills
  docs    + ./docs: the documentation of the commit (README, every published page, the tutorials as Markdown)
  skills  + ./docs + ./.agents/skills: the five U-Chrom skills (the agent-skills folder Codex reads)

    python supp_fig7/v3/run_codex.py --snapshot SNAP --data-root DATA --data-root-v3 DATA3 \\
        --image py312full.sif --venv-uchrom VENV_U --venv-plain VENV_P --codex-dir CODEX_VENDOR \\
        --questions questions_v3.yaml --ids T1,R3 --conditions none,bare,docs,skills \\
        --models gpt-6-luna --reps 1 --jobs 8 --run pilot_v3 --out OUT

Isolation: the container sees only its working directory (/work), the question's data (read-only, under
/work/data), the Python environment (/opt/venv, read-only), the Codex executable (/opt/codex, read-only), a fresh
home (/home/agent) whose Codex folder holds a copy of the user's credentials (auth.json) and this benchmark's
config.toml (no memories, no sub-agents, no web search, no apps / browser / computer use / plugins, no approvals:
the container is the sandbox), and Sherlock's managed Codex settings and agent rules (/etc/codex, /etc/agents,
read-only).  No checkout of U-Chrom, no questions, no answers.  Every session runs as `srun --exact` with its own
CPUs and memory (``mem_gb`` of the question, default --mem-gb): a session cannot take another's memory.  Packages
an agent installs with `pip install` go to its home and come after the environment's on the path (they never
replace an installed numpy); it may also make its own virtual environment.  The shell has internet access in every
condition.

Each run writes OUT/<run>/<model>/<id>__<condition>__r<k>/{events.jsonl, answer.md, stderr.txt, meta.json};
events.jsonl is Codex's JSON event stream with the time each line arrived (``_t``, s from the start).  meta.json:
wall time, the time spent in commands and the rest (the model), turns, tool calls, tokens, the skills read, whether
the docs were opened, U-Chrom imported, packages installed, and whether a usage limit stopped it (``rate_limited``:
not graded; run the same command again to fill them in).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent

#: data key -> (root, folder under the root (bound read-only at /work/data/<folder>), file the prompt names: "" = the
#: folder).  root "base" = --data-root (the shared stores), "v3" = --data-root-v3 (v3/data/*.md)
DATA = {
    # the shared stores
    "takei": ("base", "takei2025_fov0", "takei2025_fov0.chromdata.zarr"),
    "schicar": ("base", "schicar_mop", "schicar_mop.chromdata.zarr"),
    "sim": ("base", "sim_cell1", "sim_cell1.cdz"),
    "brain": ("base", "hicrna_mouse_brain", "hicrna_mouse_brain.chromdata.zarr"),
    "melanoma": ("base", "hicrna_human_melanoma", "hicrna_human_melanoma.chromdata.zarr"),
    "liu": ("base", "liu2025_mop", "liu2025_mop.chromdata.zarr"),
    "gc": ("base", "mm10_gc", "mm10_gc_500kb.tsv"),
    "stevens": ("base", "stevens2017_cell1", ""),
    # built for these questions (v3/data/*.md)
    "takei2021_fofct": ("v3", "takei2021_fofct", ""),
    "su2020_chr2p": ("v3", "su2020_chr2p", "su2020_chr2p_first100.csv"),
    "bintu2018": ("v3", "bintu2018", "bintu2018_IMR90_chr21.csv"),
    "rao_chr21_10kb": ("v3", "rao2014_imr90_chr21_10kb", ""),
    "su2020_tsv": ("v3", "su2020_chr21_tsv", ""),
    "takei_rna": ("v3", "takei2025_fov0_rna", "takei2025_fov0_rna_counts.tsv"),
    "stevens_cell2": ("v3", "stevens2017_cell2", ""),
    "mm10_cpg_1mb": ("v3", "mm10_cpg", "mm10_cpg_1Mb.tsv"),
    "schicar_types": ("v3", "schicar_cell_types", "schicar_cell_types.tsv"),
    "nagano_hap_2i": ("v3", "nagano2017_hap_2i", ""),
    "sox2_holdout": ("v3", "imgfmt_sox2", ""),
    "takei_fov4": ("v3", "imgfmt_takei2021_fov4", "takei2021_fov4.chromdata.zarr"),
    # the stores of T1 / P2 and P1 in community formats (v3/data/neutral.md)
    "takei_fofct": ("v3", "takei2025_fov0_fofct", ""),
    "melanoma_r1": ("v3", "hicrna_melanoma_r1", ""),
}
URLS = {"atlas_hires": "https://uchrom-atlas-r2.u-science.org/hires_brain.chromdata.zarr",
        "atlas_gageseq": "https://uchrom-atlas-r2.u-science.org/gageseq_mouse_cortex.chromdata.zarr",
        "atlas_hires_embryo": "https://uchrom-atlas-r2.u-science.org/hires_embryo.chromdata.zarr"}
LIBS = "numpy, pandas, scipy, zarr, pyarrow, fsspec, anndata, cooler, h5py, scikit-learn and matplotlib"
INTRO = {
    "none": f"You are analysing 3-D genome data. Python (run it with `python`) has {LIBS} installed.",
    "bare": f"You are analysing 3-D genome data with U-Chrom: the Python package `uchrom` is installed (run it with "
            f"`python`), together with {LIBS}.",
    "docs": f"You are analysing 3-D genome data with U-Chrom: the Python package `uchrom` is installed (run it with "
            f"`python`), together with {LIBS}, and its documentation is in ./docs.",
}
INTRO["skills"] = INTRO["docs"]
#: every condition may add public packages (cooltools, scanpy, ...) or make its own environment
PIP = (" You may install other public Python packages with `pip install` (they go to your home) or make your own "
       "virtual environment.")
INTRO = {k: v + PIP for k, v in INTRO.items()}
#: every data set of the benchmark is published and public (GEO, Zenodo, 4DN, the public U-Chrom atlas): Low Risk
#: under Stanford's data classification.  Said in the prompt because the agent reads Sherlock's rules (/etc/agents),
#: which tell it to stop when the classification of (human) data is unclear.
PROMPT = """{intro}

Data: {data}  (published, public data: Low Risk under Stanford's data classification)

Question: {question}

Work it out with code, then give a concise final answer: the key numbers with units, how you computed them, \
and anything that limits the conclusion."""

#: Codex configuration of every session (CODEX_HOME/config.toml)
CONFIG = """approval_policy = "never"
sandbox_mode = "danger-full-access"
web_search = "disabled"
check_for_update_on_startup = false

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
remote_plugin = false
"""
#: the agent's login profile: Codex runs every command with `bash -lc`, and a login shell's /etc/profile resets PATH
#: to the image's (whose Python has no packages); this puts the environment and Codex first again
BASH_PROFILE = """export PATH=/opt/venv/bin:/opt/codex/bin:$PATH
export VIRTUAL_ENV=/opt/venv
"""
#: appended to sys.path at start-up (PYTHONPATH=/opt/bench/site): what an agent pip-installs into its home comes
#: after the environment, so it adds packages but never replaces one
SITECUSTOMIZE = """import sys
sys.path.append("/home/agent/.pylib")
"""
LIMIT_MARKS = ("usage limit", "usage_limit", "rate limit", "rate_limit", "too many requests", "429", "quota",
               "exceeded your", "limit reached")


def data_paths(key: str, a):
    """(text for the prompt, [(host folder, container path)]) of a question's data key ('a+b' for two)."""
    texts, binds = [], []
    for k in key.split("+"):
        if k in URLS:
            texts.append(f"{URLS[k]} (served over HTTP)")
            continue
        root, folder, name = DATA[k]
        base = a.data_root if root == "base" else a.data_root_v3
        texts.append(f"./data/{folder}/{name}" if name else f"./data/{folder}/ (a folder)")
        binds.append((base / folder, f"/work/data/{folder}"))
    return "; ".join(texts), binds


def workdir(base: Path, condition: str, snapshot: Path) -> Path:
    w = Path(tempfile.mkdtemp(prefix="analysis_", dir=base))  # neutral name: no hint of a benchmark or condition
    (w / "data").mkdir()
    if condition in ("docs", "skills"):
        shutil.copytree(snapshot / "docs", w / "docs")
    if condition == "skills":
        shutil.copytree(snapshot / "skills", w / ".agents" / "skills")
    return w


def read_stream(stream, sink: list, t0: float):
    for raw in iter(stream.readline, b""):
        sink.append((time.time() - t0, raw.decode("utf-8", "replace").rstrip("\n")))


def summarise(lines, wall: float) -> dict:
    """Turns, tool calls, tokens, time in commands, what the session did (module docstring)."""
    started, cmd_s, n_cmd, n_tools, usage, msgs, errors = {}, 0.0, 0, 0, {}, [], []
    skills, docs, pip, imports, turns = set(), False, False, False, 0
    for t, line in lines:
        try:
            j = json.loads(line)
        except ValueError:
            continue
        typ, it = j.get("type"), j.get("item") or {}
        if typ == "turn.started":
            turns += 1
        if typ in ("turn.completed",):
            for k, v in (j.get("usage") or {}).items():
                usage[k] = usage.get(k, 0) + (v or 0)
        if typ in ("turn.failed", "error"):
            errors.append(json.dumps(j)[:500])
        if typ == "item.started" and it.get("id"):
            started[it["id"]] = t
        if typ == "item.completed":
            kind = it.get("type")
            if kind == "agent_message":
                msgs.append(it.get("text") or "")
            elif kind in ("command_execution", "file_change", "mcp_tool_call", "web_search"):
                n_tools += 1
            if kind == "command_execution":
                n_cmd += 1
                if it.get("id") in started:
                    cmd_s += max(0.0, t - started[it["id"]])
                c = it.get("command") or ""
                skills.update(re.findall(r"\.agents/skills/([\w.-]+)", c))
                # U-Chrom's docs only: ./docs or /work/docs (not a docs/ folder of another package's checkout)
                docs |= bool(re.search(r"(^|[\s'\"=(])(\./)?docs/|/work/docs", c))
                imports |= bool(re.search(r"\b(import|from)\s+(uchrom|chromdata)\b", c))
                pip |= "pip install" in c
    return {"answer": msgs[-1] if msgs else "", "num_turns": turns, "n_tool_calls": n_tools, "n_commands": n_cmd,
            "command_s": round(cmd_s, 1), "model_s": round(max(0.0, wall - cmd_s), 1), "tokens": usage,
            "skills_used": sorted(skills), "read_docs": docs, "pip_install": pip, "imports_uchrom": imports,
            "errors": errors}


def run_one(q, condition, rep, model, a):
    d = a.out / a.run / model / f"{q['id']}__{condition}__r{rep}"
    if (d / "meta.json").exists():
        meta = json.loads((d / "meta.json").read_text())
        if not meta.get("rate_limited"):
            return meta
    d.mkdir(parents=True, exist_ok=True)
    base = Path(a.work_root)
    w = workdir(base, condition, a.snapshot)
    home = Path(tempfile.mkdtemp(prefix="home_", dir=base))
    tmpd = Path(tempfile.mkdtemp(prefix="tmp_", dir=base))   # the container's /tmp (a contained one is tiny)
    ch = home / ".codex"
    ch.mkdir()
    shutil.copy(Path.home() / ".codex" / "auth.json", ch / "auth.json")
    os.chmod(ch / "auth.json", 0o600)
    (ch / "config.toml").write_text(CONFIG)
    for f in (".bash_profile", ".profile"):
        (home / f).write_text(BASH_PROFILE)
    text, binds = data_paths(q["data"], a)
    prompt = PROMPT.format(intro=INTRO[condition], data=text, question=q["prompt"])
    venv = a.venv_plain if condition == "none" else a.venv_uchrom
    bind = [f"{w}:/work", f"{tmpd}:/tmp", f"{venv}:/opt/venv:ro", f"{a.codex_dir}:/opt/codex:ro",
            f"{a.cfg}:/opt/bench:ro", "/etc/codex:/etc/codex:ro", "/etc/agents:/etc/agents:ro"]
    bind += [f"{src}:{target}:ro" for src, target in binds]
    env = {"HOME": "/home/agent", "CODEX_HOME": "/home/agent/.codex", "TMPDIR": "/tmp",
           "PIP_TARGET": "/home/agent/.pylib", "PYTHONPATH": "/opt/bench/site", "PIP_DISABLE_PIP_VERSION_CHECK": "1",
           "PATH": "/opt/venv/bin:/opt/codex/bin:/usr/local/bin:/usr/bin:/bin", "LANG": "C.UTF-8",
           "OMP_NUM_THREADS": str(a.threads), "OPENBLAS_NUM_THREADS": str(a.threads),
           "MKL_NUM_THREADS": str(a.threads)}
    mem = int(q.get("mem_gb") or a.mem_gb)
    # the container hides Slurm's limit (--cleanenv), so U-Chrom's streaming paths would size their batches from
    # the node's RAM: give them half the session's memory, as they would take under Slurm (no effect without uchrom)
    env["UCHROM_MEMORY_BUDGET"] = f"{max(1, mem // 2)}GB"
    # --cleanenv hides the job: the agent reads Sherlock's rules (/etc/agents) and, without SLURM_JOB_ID, takes the
    # node for a login node and refuses to compute.  The session is a job step with these resources: say so.
    if os.environ.get("SLURM_JOB_ID"):
        env.update({"SLURM_JOB_ID": os.environ["SLURM_JOB_ID"], "SLURM_CPUS_PER_TASK": str(a.threads),
                    "SLURM_MEM_PER_NODE": str(mem * 1024), "SLURM_NTASKS": "1"})
    cmd = []
    if a.srun:
        cmd += ["srun", "--exact", "--ntasks=1", f"--cpus-per-task={a.threads}", f"--mem={mem}G", "--quiet",
                f"--job-name={q['id']}-{condition}"]
    cmd += ["apptainer", "exec", "--containall", "--cleanenv", "--home", f"{home}:/home/agent", "--pwd", "/work"]
    for b in bind:
        cmd += ["--bind", b]
    for k, v in env.items():
        cmd += ["--env", f"{k}={v}"]
    cmd += [str(a.image), "/opt/codex/bin/codex", "exec", "--json", "--skip-git-repo-check", "-C", "/work",
            "-m", model, prompt]
    t0 = time.time()
    # its own process group: a time-out ends the whole group (Sherlock's Apptainer gives no PID namespace)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.DEVNULL,
                         start_new_session=True)
    lines: list = []
    errs: list = []
    readers = [threading.Thread(target=read_stream, args=(p.stdout, lines, t0), daemon=True),
               threading.Thread(target=read_stream, args=(p.stderr, errs, t0), daemon=True)]
    for r in readers:
        r.start()
    try:
        p.wait(timeout=a.timeout)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        p.wait()
        rc = "timeout"
    for r in readers:
        r.join(timeout=30)
    err = "\n".join(x for _, x in errs)
    wall = time.time() - t0
    with open(d / "events.jsonl", "w") as fh:
        for t, line in lines:
            try:
                j = json.loads(line)
                j["_t"] = round(t, 3)
                fh.write(json.dumps(j) + "\n")
            except ValueError:
                fh.write(json.dumps({"_raw": line, "_t": round(t, 3)}) + "\n")
    (d / "stderr.txt").write_text(err[-20000:])
    s = summarise(lines, wall)
    (d / "answer.md").write_text(s.pop("answer") + "\n")
    blob = (" ".join(s["errors"]) + " " + err[-4000:]).lower()
    limited = rc != "timeout" and not (d / "answer.md").read_text().strip() and any(m in blob for m in LIMIT_MARKS)
    meta = {"id": q["id"], "category": q.get("category"), "split": q.get("split"), "kind": q.get("kind"),
            "condition": condition, "rep": rep, "model_requested": model, "returncode": rc, "elapsed_s": wall,
            "rate_limited": limited, "mem_gb": mem, "threads": a.threads, "host": os.uname().nodename,
            "slurm_job": os.environ.get("SLURM_JOB_ID"), "snapshot_commit": (a.snapshot / "COMMIT").read_text().strip(),
            "codex": a.codex_version, **s}
    (d / "meta.json").write_text(json.dumps(meta, indent=1))
    for x in (w, home, tmpd):
        shutil.rmtree(x, ignore_errors=True)
    print(f"{model:12s} {q['id']:4s} {condition:6s} r{rep}  {wall:6.0f} s  rc {rc}  turns {s['num_turns']}  "
          f"cmds {s['n_commands']}  skills {s['skills_used']}  {'RATE-LIMITED' if limited else ''}", flush=True)
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--data-root-v3", type=Path, required=True)
    ap.add_argument("--image", type=Path, required=True)
    ap.add_argument("--venv-uchrom", type=Path, required=True)
    ap.add_argument("--venv-plain", type=Path, required=True)
    ap.add_argument("--codex-dir", type=Path, required=True,
                    help="the codex module's vendor folder (bin/codex, a static binary), bound at /opt/codex")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--work-root", default=os.environ.get("L_SCRATCH") or tempfile.gettempdir())
    ap.add_argument("--questions", type=Path, default=HERE / "questions_v3.yaml")
    ap.add_argument("--ids", default="")
    ap.add_argument("--conditions", default="none,bare,docs,skills")
    ap.add_argument("--models", default="gpt-6-astra,gpt-6-sol,gpt-6-luna")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--threads", type=int, default=4, help="CPUs per session")
    ap.add_argument("--mem-gb", type=int, default=24, help="memory per session unless the question says")
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--no-srun", dest="srun", action="store_false", help="run sessions without a job step each")
    ap.add_argument("--run", required=True)
    a = ap.parse_args()
    a.snapshot, a.data_root, a.data_root_v3 = a.snapshot.resolve(), a.data_root.resolve(), a.data_root_v3.resolve()
    a.cfg = Path(tempfile.mkdtemp(prefix="cfg_", dir=a.work_root))
    (a.cfg / "site").mkdir()
    (a.cfg / "site" / "sitecustomize.py").write_text(SITECUSTOMIZE)
    a.codex_version = json.loads((a.codex_dir / "codex-package.json").read_text()).get("version") \
        if (a.codex_dir / "codex-package.json").exists() else None
    qs = yaml.safe_load(a.questions.read_text())["questions"]
    ids = set(filter(None, a.ids.split(",")))
    qs = [q for q in qs if not ids or q["id"] in ids]
    tasks = [(q, c, r, m) for m in a.models.split(",") for r in range(a.reps) for q in qs
             for c in a.conditions.split(",")]
    print(f"{len(tasks)} sessions, {a.jobs} at a time, codex {a.codex_version}", flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        metas = list(ex.map(lambda t: run_one(*t, a), tasks))
    shutil.rmtree(a.cfg, ignore_errors=True)
    n_lim = sum(bool(m.get("rate_limited")) for m in metas)
    print(f"done: {len(metas)} sessions, {n_lim} rate-limited (run again to fill them in)", flush=True)


if __name__ == "__main__":
    main()
