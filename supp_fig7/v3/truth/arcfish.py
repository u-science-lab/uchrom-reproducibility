"""Reference answers of T1-T5 (v3/questions/arcfish.yaml) from the ArcFISH parity runs of the U-Chrom repository
(``benchmarks/arcfish_parity``: the same inputs through U-Chrom and through the authors' package, a21d889).

    python v3/truth/arcfish.py --runs RUN_DIR [RUN_DIR ...] --out v3/results/arcfish.json

``values`` are U-Chrom's, ``verify.arcfish_a21d889`` the authors'; a question is kept only when they agree
(``agree``).  Datasets: T1 takei2025_chr15_1mb, T2 / T3 takei2021_chr6, T4 su2020_chr2p, T5 bintu2018.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _find(runs, name, suffix):
    for r in runs:
        p = r / f"{name}.{suffix}"
        if p.exists():
            return p
    raise FileNotFoundError(f"{name}.{suffix} in {runs}")


def _tads(runs, name, kind="tads"):
    ours = pd.read_parquet(_find(runs, name, f"ours.{kind}.parquet"))
    theirs = pd.read_csv(_find(runs, name, f"upstream.{kind}.csv"))
    o0 = ours[ours["level"] == 0].sort_values("start1")
    t0 = theirs[theirs["level"] == 0].sort_values("s1")
    # interior boundaries: the start of every level-0 domain but the first
    ob = [int(x) for x in o0["start1"].to_numpy()[1:]]
    tb = [int(x) for x in t0["s1"].to_numpy()[1:]]
    return {"n_domains": len(o0), "n_boundaries": len(ob), "boundaries": ob}, \
           {"n_domains": len(t0), "n_boundaries": len(tb), "boundaries": tb}


def t1(runs):
    return _tads(runs, "takei2025_chr15_1mb")


def t2(runs):
    return _tads(runs, "takei2021_chr6")


def t5(runs):
    return _tads(runs, "bintu2018", kind="ins")


def t3(runs):
    out = []
    for side in ("ours", "upstream"):
        z = np.load(_find(runs, "takei2021_chr6", f"{side}.npz"))
        final, pval, st = np.triu(z["loop_final"].astype(bool), 1), z["loop_pval"], z["starts"]
        ij = np.argwhere(final)
        best = min(ij, key=lambda k: (pval[k[0], k[1]], k[0], k[1]))
        out.append({"n_loops": int(len(ij)), "anchor1": int(st[best[0]]), "anchor2": int(st[best[1]]),
                    "loops": [[int(st[i]), int(st[j]), float(f"{pval[i, j]:.6g}")] for i, j in ij]})
    return out[0], out[1]


def t4(runs):
    ours = pd.read_parquet(_find(runs, "su2020_chr2p", "ours.cpmt.parquet"))
    theirs = pd.read_csv(_find(runs, "su2020_chr2p", "upstream.cpmt.csv"))

    def summary(labels):
        labels = np.asarray(labels)
        sizes = sorted(int((labels == v).sum()) for v in np.unique(labels))
        runs_ = int(1 + (labels[1:] != labels[:-1]).sum())
        return {"n_smaller": sizes[0], "n_larger": sizes[-1], "n_runs": runs_, "n_loci": int(len(labels))}

    return summary(ours["cluster"].to_numpy()), summary(theirs["cpmt"].to_numpy())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    res = {}
    for qid, fn in (("T1", t1), ("T2", t2), ("T3", t3), ("T4", t4), ("T5", t5)):
        ours, theirs = fn(a.runs)
        res[qid] = {"values": ours, "verify": {"arcfish_a21d889": theirs}, "truth": f"arcfish_{qid.lower()}",
                    "agree": ours == theirs}
        print(qid, "agree" if ours == theirs else "DIFFER", json.dumps(ours)[:200])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
