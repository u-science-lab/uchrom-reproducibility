"""Summarise the graded sessions of agent benchmark version 3 into results/agent_bench_v3.json.

    python supp_fig7/v3/summarize_v3.py [--grades supp_fig7/v3/results/main_v3_grades.json]

The question is the unit: a question's score in a condition is the mean over its runs, first within each model, then
over the models.  Conditions are compared over questions with two-sided Wilcoxon signed-rank tests (questions whose
two values are equal drop out), for all 27 questions, Part A (anchor and format questions: reading the data) and Part
B (method questions: a published method, referenced by U-Chrom and by the authors' code).  The same summary is also
given for the models with every session, and per model; timings, tokens and skill use per condition.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "shared"))
from _record import save  # noqa: E402

CONDS = ["none", "bare", "docs", "skills"]
PAIRS = [("none", "bare"), ("none", "docs"), ("none", "skills"), ("bare", "docs"), ("docs", "skills")]
PARTS = {"all": None, "A": ("anchor", "format"), "B": ("method",)}


def question_means(grades: list[dict], models: set[str]) -> dict[str, dict[str, float]]:
    runs = defaultdict(list)
    for g in grades:
        if g["model_requested"] in models:
            runs[(g["id"], g["condition"], g["model_requested"])].append(g["score"])
    per_model = defaultdict(lambda: defaultdict(list))
    for (q, c, _m), v in runs.items():
        per_model[q][c].append(float(np.mean(v)))
    return {q: {c: float(np.mean(v)) for c, v in cs.items()} for q, cs in per_model.items()}


def compare(qm: dict, kind: dict) -> dict:
    out = {}
    for part, kinds in PARTS.items():
        qs = sorted(q for q in qm if kinds is None or kind[q] in kinds)
        res = {"n_questions": len(qs), "mean": {c: float(np.mean([qm[q][c] for q in qs])) for c in CONDS}, "tests": {}}
        for a, b in PAIRS:
            x = np.array([qm[q][a] for q in qs]); y = np.array([qm[q][b] for q in qs])
            differ = int(np.sum(~np.isclose(x, y)))
            p = float(wilcoxon(x, y).pvalue) if differ else 1.0
            res["tests"][f"{a}_vs_{b}"] = {"mean_diff": float(np.mean(y - x)), "questions_differing": differ, "p_wilcoxon": p}
        out[part] = res
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--grades", type=Path, default=HERE / "results" / "main_v3_grades.json")
    ap.add_argument("--out", type=Path, default=HERE.parent / "results" / "agent_bench_v3.json")
    a = ap.parse_args()
    grades = list(json.loads(a.grades.read_text()).values())
    kind = {g["id"]: g["kind"] for g in grades}
    models = sorted({g["model_requested"] for g in grades})
    n = defaultdict(int)
    for g in grades:
        n[(g["model_requested"], g["condition"])] += 1
    complete = {m for m in models if all(n[(m, c)] == max(n.values()) for c in CONDS)}
    per_q = {}
    for q in sorted(kind):
        cell = defaultdict(list)
        for g in grades:
            if g["id"] == q:
                cell[g["condition"]].append(g["score"])
        per_q[q] = {"kind": kind[q], **{c: float(np.mean(cell[c])) for c in CONDS}}
    per_model = {m: {c: {"mean": float(np.mean([g["score"] for g in grades if g["model_requested"] == m and g["condition"] == c])),
                         "n": n[(m, c)]} for c in CONDS} for m in models}
    usage = {}
    for c in CONDS:
        gs = [g for g in grades if g["condition"] == c]
        tok = [sum(v for v in g["tokens"].values() if isinstance(v, (int, float))) if isinstance(g.get("tokens"), dict)
               else (g.get("tokens") or 0) for g in gs]
        usage[c] = {"median_elapsed_s": float(np.median([g.get("elapsed_s") or 0 for g in gs])),
                    "median_tokens": float(np.median(tok)),
                    "share_pip_install": float(np.mean([bool(g.get("pip_install")) for g in gs])),
                    "share_skill_loaded": float(np.mean([bool(g.get("skills_used")) for g in gs]))}
    res = {
        "sessions": len(grades), "sessions_per_model_condition": {f"{m}/{c}": n[(m, c)] for m in models for c in CONDS},
        "models": models, "judge_models": sorted({g.get("judge_model", "") for g in grades}),
        "all_models": compare(question_means(grades, set(models)), kind),
        "complete_models": {"models": sorted(complete), **compare(question_means(grades, complete), kind)},
        "per_model": per_model, "per_question": per_q, "usage": usage,
    }
    save(a.out, res, script=__file__, params=vars(a) | {"grades": str(a.grades), "out": str(a.out)},
         datasets=["takei2025_fov0", "schicar_mop", "atlas stores and benchmark inputs of v3/data/*.md"])
    for part in PARTS:
        r = res["all_models"][part]
        print(f"{part:3s} n={r['n_questions']:2d} " + " ".join(f"{c}={r['mean'][c]:.3f}" for c in CONDS)
              + "  | " + " ".join(f"{k}: P={t['p_wilcoxon']:.3f}" for k, t in r["tests"].items()))
    print(a.out)


if __name__ == "__main__":
    main()
