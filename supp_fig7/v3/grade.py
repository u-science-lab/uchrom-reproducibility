"""Grade the agent benchmark's sessions against answers_v3.json with an LLM judge.

    python supp_fig7/v3/grade.py --runs OUT/<run> [--judge-model claude-opus-5-5] [--jobs 8]
    python supp_fig7/v3/grade.py --runs OUT/<run> --judge-model claude-sonnet-5-5 --sample 0.1 --out grades_sonnet.json
        (a second judge on a random tenth)

Each check is binary: a number within tolerance of the reference, a number that reaches a floor, a named text answer,
or a rubric item (RUBRIC below).  For numbers the judge only extracts the value the answer reports, converted to the
reference's unit and scale; the tolerance test is done here in code (`number_pass`), and the judge's own verdict is
kept as `judge_pass` for auditing.  The judge sees the question, the reference values and the answer, not the
condition.  It is a fresh `claude -p` session without tools, memory, user settings or MCP servers.  Runs stopped by a
usage limit are left out (listed in the output); a run that ended without an answer (time-out, error) scores zero.
Writes OUT/<run>/grades.json: one entry per session with its checks, score, model, condition and split (--out:
another file in OUT/<run>; --sample: a random fraction of the sessions, seed 0).  Score of a session = the fraction of
its checks passed.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml


#: rubric items (the questions name them; the prompt gives the judge their definitions)
RUBRIC = {
    "boundaries_match": "the boundary positions the answer reports match the reference's \"boundaries\" (bp): the "
                        "same number of boundaries, and each within one bin of the data of a reference boundary "
                        "(1 Mb for 1-Mb bins, 25 kb / 30 kb for 25-kb / 30-kb loci); a boundary given as the end of a "
                        "domain equal to the next domain's start is the same boundary",
}


def _build_prompt(q, ref, answer):
    g = q.get("grade", {})
    checks = []
    for name in g.get("numbers", []):
        checks.append(f'- number "{name}": reference {ref.get(name)!r}; correct if within ±{g.get("tol", 0) * 100:g}% '
                      f'(or exact if tolerance is 0).  Also give "value": the number the answer reports as its result '
                      f'for this quantity (its headline value, not an alternative it mentions), converted to the same '
                      f'unit and scale as the reference (e.g. 25% -> 0.25 if the reference is a fraction); null if the '
                      f'answer gives none')
    txt = g.get("text")
    for name in ([txt] if isinstance(txt, str) else (txt or [])):
        checks.append(f'- text "{name}": reference {ref.get(name, name)!r}; correct if the answer names it')
    for r in g.get("rubric", []):
        checks.append(f'- rubric "{r}": {RUBRIC.get(r, r)}')
    return f"""You are grading an answer to a data-analysis question. Be strict and literal.

Question: {q['prompt']}

Reference values (computed by the benchmark authors): {json.dumps(ref)}

Checks:
{chr(10).join(checks)}

Answer to grade:
<<<
{answer}
>>>

Reply with JSON only: {{"checks": {{"<check name>": {{"pass": true|false, "given": "<what the answer said, short>", "value": <number or null, number checks only>}}}}, "comment": "<one sentence>"}}"""


def number_pass(value, ref, tol: float) -> bool:
    """Tolerance test for a number check: relative ``tol``; tol 0 means exact (integers) or 1e-9 relative."""
    if value is None or isinstance(value, bool):
        return False
    try:
        v, r = float(value), float(ref)
    except (TypeError, ValueError):
        return False
    if tol == 0:
        return round(v) == r if float(r).is_integer() else abs(v - r) <= 1e-9 * max(abs(r), 1e-12)
    return abs(v - r) <= tol * abs(r)


def build_prompt(q, ref, answer):
    """the prompt, plus the "numbers_min" checks (passed when the value reaches a floor, e.g. a reconstruction
    that agrees with the published structure at least as well as the reference's floor)"""
    p = _build_prompt(q, ref, answer)
    extra = [f'- number "{n}": correct if at least {ref.get(n + "_floor")!r}.  Also give "value": the number the answer '
             f'reports for this quantity; null if none' for n in q.get("grade", {}).get("numbers_min", [])]
    return p.replace("\nAnswer to grade:", "\n" + "\n".join(extra) + "\n\nAnswer to grade:") if extra else p

HERE = Path(__file__).resolve().parent


class JudgeError(RuntimeError):
    pass


def judge(prompt: str, model: str, cfg: Path, tries: int = 3) -> dict:
    cmd = ["claude", "-p", prompt, "--output-format", "json", "--model", model, "--setting-sources", "project",
           "--settings", str(cfg / "settings.json"), "--strict-mcp-config", "--mcp-config", str(cfg / "mcp.json"),
           "--disallowedTools", "Bash", "Read", "Write", "Edit", "Glob", "Grep", "WebFetch", "WebSearch", "Skill"]
    err, text = None, ""
    for _ in range(tries):
        with tempfile.TemporaryDirectory(prefix="judge_") as cwd:
            try:
                out = subprocess.run(cmd, capture_output=True, text=True, timeout=300, cwd=cwd).stdout
                text = json.loads(out).get("result", "")
                return json.loads(re.search(r"\{.*\}", text, re.S).group(0))
            except (subprocess.TimeoutExpired, ValueError, AttributeError) as e:
                err = e
    raise JudgeError(f"judge failed {tries} times: {err!r}; last reply: {text[:300]!r}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--answers", type=Path, default=HERE / "answers_v3.json")
    ap.add_argument("--questions", type=Path, default=HERE / "questions_v3.yaml")
    ap.add_argument("--judge-model", default="claude-opus-5-5")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", default="grades.json", help="file name in the run's folder")
    ap.add_argument("--sample", type=float, help="grade a random fraction of the sessions (seed 0)")
    a = ap.parse_args()
    qs = {q["id"]: q for q in yaml.safe_load(a.questions.read_text())["questions"]}
    answers = json.loads(a.answers.read_text())
    out = a.runs / a.out
    grades = json.loads(out.read_text()) if out.exists() else {}
    cfg = Path(tempfile.mkdtemp(prefix="judgecfg_"))
    (cfg / "settings.json").write_text(json.dumps({"autoMemoryEnabled": False}))
    (cfg / "mcp.json").write_text(json.dumps({"mcpServers": {}}))
    runs = sorted(p for p in a.runs.glob("*/*") if (p / "meta.json").exists())
    if a.sample:
        import random
        runs = sorted(random.Random(0).sample(runs, round(a.sample * len(runs))))
    todo, limited = [], []
    for p in runs:
        key = f"{p.parent.name}/{p.name}"
        meta = json.loads((p / "meta.json").read_text())
        if meta.get("rate_limited"):
            limited.append(key)
        elif key not in grades:
            todo.append((key, p, meta))

    def grade_one(item):
        key, p, meta = item
        q, ref = qs[meta["id"]], answers[meta["id"]]["values"]
        answer = (p / "answer.md").read_text().strip()
        spec = q.get("grade", {})
        if not answer:
            names = list(spec.get("numbers", [])) + list(spec.get("numbers_min", [])) + ([spec["text"]] if isinstance(spec.get("text"), str) else
                                                     list(spec.get("text") or [])) + list(spec.get("rubric", []))
            checks = {n: {"pass": False, "given": "(no answer)"} for n in names}
            return key, {**meta, "score": 0.0, "checks": checks, "comment": "no answer"}
        for _ in range(3):   # a reply whose checks are not objects (e.g. "pass") is judged again
            g = judge(build_prompt(q, ref, answer), a.judge_model, cfg)
            checks = g.get("checks", {}) if isinstance(g, dict) else None
            if isinstance(checks, dict) and all(isinstance(v, dict) for v in checks.values()):
                break
        else:
            raise JudgeError(f"judge returned malformed checks 3 times: {str(g)[:300]!r}")
        for name in spec.get("numbers", []):
            c = checks.setdefault(name, {"pass": False, "given": "(judge returned no check)"})
            c["judge_pass"] = bool(c.get("pass"))
            c["pass"] = number_pass(c.get("value"), ref.get(name), (spec.get("tol_by_number") or {}).get(name,
                                                                                         spec.get("tol", 0)))
        for name in spec.get("numbers_min", []):
            c = checks.setdefault(name, {"pass": False, "given": "(judge returned no check)"})
            c["judge_pass"] = bool(c.get("pass"))
            try:
                c["pass"] = c.get("value") is not None and float(c["value"]) >= float(ref[name + "_floor"])
            except (TypeError, ValueError):
                c["pass"] = False
        score = sum(bool(v.get("pass")) for v in checks.values()) / max(len(checks), 1)
        return key, {**meta, "score": score, "checks": checks, "comment": g.get("comment", ""),
                     "judge_model": a.judge_model}

    def try_one(item):
        try:
            return grade_one(item)
        except JudgeError as e:       # e.g. the judge's own usage limit: left ungraded, a later call fills it in
            return item[0], e

    failed = []
    with ThreadPoolExecutor(a.jobs) as ex:
        for key, g in ex.map(try_one, todo):
            if isinstance(g, JudgeError):
                failed.append(key)
                print(f"{key}: {g}", flush=True)
                continue
            grades[key] = g
            out.write_text(json.dumps(grades, indent=1))
    print(f"graded {len(todo) - len(failed)} new, {len(grades)} in all; {len(failed)} left ungraded (judge failed): "
          f"{failed[:10]}; {len(limited)} rate-limited runs left out: {limited[:10]}")


if __name__ == "__main__":
    main()
