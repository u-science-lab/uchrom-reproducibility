"""Agreement of two judges of benchmark version 2 on the sessions both graded (grade_v2.py --out).

    python benchmarks/fig6/agent_bench/judge_agreement.py --runs OUT/<run> [--a grades.json] [--b grades_sonnet.json]

Per check (number, named answer, rubric item): the fraction on which the judges agree and Cohen's kappa; numbers are
the tolerance test in code on the value each judge extracted, so a disagreement there is a different extraction.  Per
session: the fraction with the same score and the mean absolute difference; per condition the mean score by either
judge (a judge that favours one condition would show here).
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


def kappa(pairs):
    n = len(pairs)
    po = sum(x == y for x, y in pairs) / n
    pa, pb = sum(x for x, _ in pairs) / n, sum(y for _, y in pairs) / n
    pe = pa * pb + (1 - pa) * (1 - pb)
    return (po - pe) / (1 - pe) if pe < 1 else 1.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--a", default="grades.json")
    ap.add_argument("--b", default="grades_sonnet.json")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    ga, gb = (json.loads((a.runs / f).read_text()) for f in (a.a, a.b))
    keys = sorted(set(ga) & set(gb))
    checks, kinds, by_cond = [], defaultdict(list), defaultdict(lambda: [[], []])
    same = diff = 0.0
    for k in keys:
        x, y = ga[k], gb[k]
        same += x["score"] == y["score"]
        diff += abs(x["score"] - y["score"])
        c = (x["model_requested"], x["condition"])
        by_cond[c][0].append(x["score"])
        by_cond[c][1].append(y["score"])
        for name in set(x["checks"]) & set(y["checks"]):
            p = (bool(x["checks"][name].get("pass")), bool(y["checks"][name].get("pass")))
            checks.append(p)
            kinds["number" if "value" in x["checks"][name] else "judged"].append(p)
    res = {"n_sessions": len(keys), "n_checks": len(checks),
           "check_agreement": sum(p == q for p, q in checks) / len(checks), "check_kappa": kappa(checks),
           "by_kind": {k: {"n": len(v), "agreement": sum(p == q for p, q in v) / len(v), "kappa": kappa(v)}
                       for k, v in kinds.items()},
           "session_same_score": same / len(keys), "session_mean_abs_diff": diff / len(keys),
           "by_condition": {f"{m} {c}": {"n": len(s[0]), "a": sum(s[0]) / len(s[0]), "b": sum(s[1]) / len(s[1])}
                            for (m, c), s in sorted(by_cond.items())}}
    print(json.dumps(res, indent=1))
    if a.out:
        a.out.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
