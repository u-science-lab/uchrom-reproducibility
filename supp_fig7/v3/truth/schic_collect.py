"""Merge the U-Chrom references and the independent implementations of the schic family into results/schic.json
(the format of answers_v2.json: {id: {"values", "verify", "truth", "runtime_s"}}).

    python schic_collect.py --work results/schic_work --out results/schic.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def _load(p: Path) -> dict:
    return json.loads(p.read_text()) if p.exists() else {}


def _k30(r: dict) -> dict:
    """R4 runs made with k = 10 as the main value (baseline, variants): report k = 30 as ``purity`` like the graded run."""
    r = {k: v for k, v in r.items() if k not in ("runtime_s_by_chrom", "wall_s")}
    if r.get("k") == 10 and "purity_k30" in r:
        k10, k30 = r.pop("purity"), r.pop("purity_k30")
        pooled = r.pop("purity_L23IT_pooled", None)
        r = {**r, "k": 30, "purity": k30, "purity_k10": k10, "purity_L23IT_pooled_k10": pooled}
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    uc = {**_load(a.work / "schic_uchrom.json"), **_load(a.work / "schic_uchrom_r7.json")}
    sc = _load(a.work / "schic_ref_schicluster.json")
    npy = _load(a.work / "schic_ref_numpy.json")
    sn = _load(a.work / "schic_ref_snaphic.json")
    out = {}
    for q, rec in uc.items():
        vals = dict(rec["values"])
        rt = vals.pop("runtime_s", {})
        entry = {"values": vals, "verify": {}, "truth": rec["truth"],
                 "runtime_s": dict(rt) if isinstance(rt, dict) else {"uchrom": rt}}
        if q in ("R1", "R2"):
            for tag, label in (("", "nuc_dynamics"), ("_def", "nuc_dynamics"), ("_t100", "nuc_dynamics_temps100"),
                               ("_t50d50", "nuc_dynamics_temps50_dyns50")):
                r = _load(a.work / f"schic_ref_nucdyn{tag}.json").get(q + tag)
                if not r or not any("rho" in x for x in r["runs"]) or label in entry["verify"]:
                    continue                          # "_def": the default protocol run in a job of its own
                entry["verify"][label] = {"rho": r.get("rho"), "rho_by_seed": [x.get("rho") for x in r["runs"]],
                                          "seeds": [x.get("seed") for x in r["runs"]], "extra_args": r.get("extra_args", ""),
                                          "commit": r.get("nuc_dynamics_commit"), "n_contacts": r.get("n_contacts")}
                entry["runtime_s"][label + "_per_run"] = [x.get("seconds") for x in r["runs"]]
        if q in ("R3", "R4") and q in sc:
            r = {k: v for k, v in sc[q].items() if k not in ("runtime_s", "runtime_s_by_chrom", "wall_s")}
            entry["verify"]["schicluster"] = r
            entry["runtime_s"]["schicluster"] = sc[q].get("runtime_s")
        if q == "R4" and "R4" in sc and "purity" in vals and "purity" in sc["R4"]:
            # floor: U-Chrom's purity (k = 30) minus 0.01, below every correct variant (reports/schic.md)
            vals["purity_floor"] = round(vals["purity"] - 0.01, 3)
        if q == "R4" and "R4v" in sc:
            entry["verify"]["schicluster_variants"] = {name: _k30(r) for name, r in sc["R4v"].items()
                                                       if isinstance(r, dict)}
        if q == "R4" and "R4b" in sc:
            entry["verify"]["baseline_no_imputation"] = _k30(sc["R4b"])
        if q == "R5":
            if q in npy:
                entry["verify"]["schicluster+cooltools"] = {k: v for k, v in npy[q].items()
                                                            if k not in ("runtime_s_e1", "wall_s")}
                entry["runtime_s"]["cooltools_e1"] = npy[q].get("runtime_s_e1")
            if q in sc:
                entry["runtime_s"]["schicluster_scores"] = sc[q].get("runtime_s")
        if q == "R6" and q in npy:
            entry["verify"]["numpy"] = {k: v for k, v in npy[q].items() if k not in ("published", "runtime_s", "wall_s")}
            entry["verify"]["published_values"] = npy[q].get("published")
            entry["runtime_s"]["numpy"] = npy[q].get("runtime_s")
        if q == "R7" and sn:
            entry["verify"]["snaphic"] = {k: v for k, v in sn.items() if k not in ("runtime_s", "stderr", "stdout",
                                                                                 "log_tail")}
            entry["runtime_s"]["snaphic"] = sn.get("runtime_s")
            first = _load(a.work / "schic_ref_snaphic_newest_pandas3.json")
            if first:
                entry["verify"]["snaphic_newest_packages"] = {"exit": first.get("exit"),
                                                              "runtime_s_until_failure": first.get("runtime_s"),
                                                              "error": (first.get("stderr") or "").strip()[-200:]}
        if q == "R7":
            entry["status"] = "dropped: the authors' SnapHiC does not finish on Python 3.12 (reports/schic.md)"
        out[q] = entry
    a.out.write_text(json.dumps(out, indent=1, default=str) + "\n")
    print(json.dumps({q: {"values": {k: v for k, v in e["values"].items() if not isinstance(v, (dict, list))},
                          "verify": list(e["verify"])} for q, e in out.items()}, indent=1, default=str))


if __name__ == "__main__":
    main()
