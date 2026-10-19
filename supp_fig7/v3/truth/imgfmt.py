"""Reference answers of family ``imgfmt`` (version 3): T6 (Part B, imaging methods; T7 dropped), F1-F4 (Part A, reading the
open store) and the anchors A1, A6, H1, S3.

The implementations run in different environments (``v3/sherlock/imgfmt_truth.sbatch``), each writing one JSON into
a work directory:

* U-Chrom (the agents' ``venv-uchrom``, commit fb43ac4): ``imgfmt_t6.py --impl uchrom``, ``imgfmt_t7.py --impl
  uchrom_*``, ``imgfmt_fmt.py --impl chromdata|uchrom``;
* independent: the authors' SnapFISH-IMPUTE (``venv-sfimpute``: its functions in one process, and its MPI script
  ``run_impute.py``), the authors' jie (``venv-jie``), plain baselines (``venv-imgfmt``), and the raw Zarr / Parquet
  readers written from ``spec.md`` without chromdata (``venv-plain``, the agents' no-U-Chrom environment).

One function per question (the ``truth:`` of ``questions/imgfmt.yaml``) turns those JSONs into the entry of
``results/imgfmt.json`` (``{"values", "verify", "truth", "runtime_s"}``, as ``answers_v3.json``):

    python imgfmt.py --work WORK --anchors ../answers_anchors.json --out results/imgfmt.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

#: grading parameters set from the measurements (see reports/imgfmt.md)
T7_FLOOR_MARGIN = 0.10          # median ARI floor = min(U-Chrom, authors' jie) - margin


def _load(work: Path, name: str) -> dict:
    p = work / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else {}


def _pick(d: dict, keys) -> dict:
    return {k: d.get(k) for k in keys}


def _rel(a, b):
    try:
        return abs(a - b) / abs(b) if b else (0.0 if a == b else float("inf"))
    except TypeError:
        return None


# ------------------------------------------------------------------------------------------------------- T6
def t6_snapfish_holdout(work: Path) -> dict:
    keys = ["rmse_nm", "median_error_nm", "mean_error_nm", "n_holdout", "n_not_imputed"]
    u, f, m, lin = (_load(work, f"t6_{n}") for n in ("uchrom", "sfi_func", "sfi_mpi", "linear"))
    values = _pick(u, keys)
    values["linear_interpolation_rmse_nm"] = lin.get("rmse_nm")
    values["linear_interpolation_median_error_nm"] = lin.get("median_error_nm")
    verify = {"sfimpute_0.0.2_functions": {**_pick(f, keys), "n_iterations": f.get("n_iterations")},
              "sfimpute_0.0.2_run_impute_mpi_resize": _pick(m, keys)}
    agree = {k: {"functions": _rel(f.get(k), u.get(k)), "mpi_resize": _rel(m.get(k), u.get(k))}
             for k in ("rmse_nm", "median_error_nm")}
    return {"values": values, "verify": verify, "agreement_rel": agree, "truth": "t6_snapfish_holdout",
            "runtime_s": {"uchrom": u.get("seconds"), "sfimpute_functions": f.get("seconds"),
                          "sfimpute_mpi": m.get("seconds"), "processes": u.get("jobs")}}


# ------------------------------------------------------------------------------------------------------- T7
def t7_jie_assignment(work: Path, tag: str = "chr3_chr7_chr11") -> dict:
    keys = ["assigned", "purity", "median_ari", "mean_ari", "n_scored_spots", "n_fibres", "n_groups_ari"]
    runs = {n: _load(work, f"t7_{tag}_{n}") for n in
            ("uchrom_default", "uchrom_tuned", "jie_default", "jie_nb", "kmeans2", "dbscan1000", "nnchain", "all_one")}
    u = runs["uchrom_tuned"]
    jie_best = max((runs["jie_default"], runs["jie_nb"]), key=lambda r: r.get("median_ari") or -1)
    floor = None
    if u and jie_best:
        floor = round(min(u["median_ari"], jie_best["median_ari"]) - T7_FLOOR_MARGIN, 2)
    values = {**_pick(u, keys), "median_ari_floor": floor}
    verify = {name: {**_pick(r, keys), "seconds": r.get("seconds")} for name, r in runs.items() if r}
    return {"values": values, "verify": verify, "truth": "t7_jie_assignment",
            "runtime_s": {n: r.get("seconds") for n, r in runs.items() if r}}


# ---------------------------------------------------------------------------------------------------- F1-F4
def _fmt(work: Path, q: str, a: str, b: str, truth: str, keys=None) -> dict:
    x, y = _load(work, f"{q}_{a}"), _load(work, f"{q}_{b}")
    vx, vy = x.get("values", {}), y.get("values", {})
    if keys:
        vx_k, vy_k = _pick(vx, keys), _pick(vy, keys)
    else:
        vx_k, vy_k = vx, vy
    same = {k: (vx_k.get(k) == vy_k.get(k)) if not isinstance(vx_k.get(k), float)
            else (_rel(vy_k.get(k), vx_k.get(k)) or 0) < 1e-9 for k in vx_k}
    return {"values": vx, "verify": {f"raw_{b}": vy_k}, "agree": same, "truth": truth,
            "runtime_s": {a: x.get("seconds"), b: y.get("seconds")},
            "peak_rss_mb": {a: x.get("peak_rss_mb"), b: y.get("peak_rss_mb")}}


def f1_embedded_contacts(work: Path) -> dict:
    return _fmt(work, "f1", "chromdata", "raw", "f1_embedded_contacts",
                ["total_cis", "pixel", "row_sum", "n_pixels", "cell_index"])


def f2_cell_outlines(work: Path) -> dict:
    return _fmt(work, "f2", "chromdata", "raw", "f2_cell_outlines",
                ["n_outlines", "area", "largest_cell", "largest_area", "unit", "n_cells_in_store"])


def f3_provenance(work: Path) -> dict:
    return _fmt(work, "f3", "uchrom", "raw", "f3_provenance",
                ["function", "stored_d_short", "stored_n_pairs_short", "recomputed_d_short",
                 "recomputed_n_pairs_short", "recomputed_d_long", "recomputed_n_pairs_long"])


def f4_large_store_http(work: Path) -> dict:
    out = _fmt(work, "f4", "chromdata", "raw", "f4_large_store_http")
    keep = ["n_traces", "n_cells", "median_rg_all", "largest_median_type_100cells", "smallest_median_type_100cells"]
    v = out["values"]
    # the graded per-type medians under plain names (median_rg[early neurons] -> median_rg_early_neurons)
    per_type = {"median_rg_" + k[len("median_rg["):-1].replace(" ", "_"): v[k] for k in v if k.startswith("median_rg[")}
    out["values"] = {**_pick(v, keep), **per_type, "by_type": v.get("by_type")}
    return out


# ---------------------------------------------------------------------------------------------------- anchors
def anchors(path: Path) -> dict:
    a = json.loads(path.read_text())
    return {k: {**a[k], "anchor": "answers_anchors.json, unchanged"} for k in ("A1", "A6", "H1", "S3")}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--anchors", type=Path, required=True, help="answers_anchors.json")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    # T7 was dropped from version 3 (plain k-means separates the homologs better than the authors' jie: the question
    # would not need the method); t7_jie_assignment stays for the record of the exploration runs.
    res = {"T6": t6_snapfish_holdout(a.work), "F1": f1_embedded_contacts(a.work),
           "F2": f2_cell_outlines(a.work), "F3": f3_provenance(a.work), "F4": f4_large_store_http(a.work),
           **anchors(a.anchors)}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str)[:6000])


if __name__ == "__main__":
    main()
