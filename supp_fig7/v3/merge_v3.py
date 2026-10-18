"""Merge the families of version 3 (v3/questions/<family>.yaml, v3/results/<family>.json) into questions_v3.yaml and
answers_v3.json, renaming ids that collide (the bulk family's H1-H4 become BH1-BH4: H1 / H2 are Part A's HiRES
questions).  Questions marked ``drop: true`` are left out.

    python benchmarks/fig6/agent_bench/merge_v3.py
"""

from __future__ import annotations

import json
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
RENAME = {"bulk": {"H1": "BH1", "H2": "BH2", "H3": "BH3", "H4": "BH4"}}
ORDER = ["imgfmt", "arcfish", "schic", "bulk"]


def main():
    qs, ans = [], {}
    for fam in ORDER:
        qf, rf = HERE / "v3" / "questions" / f"{fam}.yaml", HERE / "v3" / "results" / f"{fam}.json"
        if not qf.exists():
            print("missing", qf)
            continue
        d = yaml.safe_load(qf.read_text())
        items = d["questions"] if isinstance(d, dict) else d
        res = json.loads(rf.read_text()) if rf.exists() else {}
        for q in items:
            if q.get("drop"):
                continue
            old = q["id"]
            new = RENAME.get(fam, {}).get(old, old)
            q = {**q, "id": new, "family": fam}
            if old != new:
                q["id_in_family"] = old
            qs.append(q)
            if old not in res:
                print(f"{fam}/{old}: no reference")
                continue
            ans[new] = {**res[old], "family": fam}
    ids = [q["id"] for q in qs]
    assert len(ids) == len(set(ids)), "duplicate ids"
    (HERE / "questions_v3.yaml").write_text(yaml.safe_dump({"questions": qs}, sort_keys=False, width=110,
                                                           allow_unicode=True))
    (HERE / "answers_v3.json").write_text(json.dumps(ans, indent=1))
    print(len(qs), "questions:", ", ".join(ids))


if __name__ == "__main__":
    main()
