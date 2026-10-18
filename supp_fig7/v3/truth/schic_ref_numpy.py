"""The independent references of R5 (E1 with cooltools) and R6 (numpy from the paper's definition, and the authors'
published per-cell values), without U-Chrom (environment venv-verify: cooltools 0.7.1, cooler, numpy, pandas).

    python schic_ref_numpy.py --data /data --data3 /data_v3 --shared /shared --work DIR --out schic_ref_numpy.json
        [--only R5,R6]

* R5: the cis E1 of ``contacts_all.mcool::/resolutions/1000000`` (stored weights) with ``cooltools.eigs_cis``
  (ignore_diags 2, clip 99.9, view chr1-19), each chromosome's E1 flipped to correlate positively (Pearson) with the
  CpG ratio; Pearson r with the mean scHiCluster compartment score written by schic_ref_schicluster.py (bins with
  CpG ratio > 0 and an E1).
* R6: Nagano et al. 2017 (schic2, ``analyzeScHiC.r`` / ``ordering.r`` / ``hap_serum_2i_params.r``): each fend pair
  once, fends placed with ``GATC.fends`` (mm9), cis contacts of chr1-19 and X, distance d = |coord2 - coord1|;
  among the contacts with d > 2^14.5: 'near' 2^14.5 < d <= 2^21, 'mitotic band' 2^21 < d <= 2^23.5 (R ``cut``,
  right-closed); the 'newer' phase gates.  Also the published values of the GEO feature table for the same cells.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

AUTOSOMES = [f"chr{i}" for i in range(1, 20)]
CHROMS6 = [f"chr{i}" for i in range(1, 20)] + ["chrX"]


def r5(a) -> dict:
    import bioframe
    import cooler
    import cooltools

    t0 = time.perf_counter()
    clr = cooler.Cooler(f"{a.data / 'schicar_mop/contacts_all.mcool'}::/resolutions/1000000")
    cpg = pd.read_csv(a.data3 / "mm10_cpg/mm10_cpg_1Mb.tsv", sep="\t")
    bins = clr.bins()[:][["chrom", "start", "end"]]
    bins["chrom"] = bins["chrom"].astype(str)
    track = bins.merge(cpg, on=["chrom", "start", "end"], how="left")
    track = track[track["chrom"].isin(AUTOSOMES)].reset_index(drop=True)
    view = bioframe.make_viewframe([(c, 0, int(clr.chromsizes[c]), c) for c in AUTOSOMES])
    _, vecs = cooltools.eigs_cis(clr, phasing_track=track[["chrom", "start", "end", "cpg_ratio"]], view_df=view,
                                 n_eigs=3, ignore_diags=2, clip_percentile=99.9)
    vecs["chrom"] = vecs["chrom"].astype(str)
    vecs = vecs[vecs["chrom"].isin(AUTOSOMES)].merge(track, on=["chrom", "start", "end"], how="left")
    e1 = vecs["E1"].to_numpy(float).copy()
    for c in AUTOSOMES:                                  # orient each chromosome by Pearson with the CpG ratio
        m = (vecs["chrom"] == c).to_numpy() & np.isfinite(e1) & (vecs["cpg_ratio"].to_numpy() > 0)
        if m.sum() > 2 and np.corrcoef(e1[m], vecs["cpg_ratio"].to_numpy()[m])[0, 1] < 0:
            e1[(vecs["chrom"] == c).to_numpy()] *= -1
    vecs["E1"] = e1
    secs = time.perf_counter() - t0
    vecs[["chrom", "start", "end", "E1", "cpg_ratio"]].to_csv(a.work / "schic_R5_cooltools_e1.tsv", sep="\t", index=False)
    out = {"cooltools": cooltools.__version__, "runtime_s_e1": secs}
    mean = a.work / "schic_R5_schicluster_mean.tsv"
    if mean.exists():
        m = pd.read_csv(mean, sep="\t").merge(vecs[["chrom", "start", "end", "E1", "cpg_ratio"]],
                                             on=["chrom", "start", "end"], how="left")
        ok = np.isfinite(m["E1"]) & (m["cpg_ratio"] > 0)
        S, E, C = (m.loc[ok, c].to_numpy(float) for c in ("mean_comp", "E1", "cpg_ratio"))

        def resid(y, x):                                 # least-squares line with intercept
            A = np.column_stack([np.ones_like(x), x])
            return y - A @ np.linalg.lstsq(A, y, rcond=None)[0]

        out["r_partial_cpg"] = float(np.corrcoef(resid(S, C), resid(E, C))[0, 1])
        out["r"] = float(np.corrcoef(S, E)[0, 1])
        out["n_bins"] = int(ok.sum())
        out["r_cpg_e1"] = float(np.corrcoef(C, E)[0, 1])
        out["r_score_cpg"] = float(np.corrcoef(S, C)[0, 1])
    ub = a.work / "schic_R5_uchrom_bins.tsv"
    if ub.exists():                                    # U-Chrom's mean against the cooltools E1, and the two E1s
        u = pd.read_csv(ub, sep="\t")
        u[["chrom", "rng"]] = u["bin"].str.split(":", expand=True)
        u["start"] = u["rng"].str.split("-").str[0].astype(int)
        u = u.merge(vecs[["chrom", "start", "E1"]].rename(columns={"E1": "E1_cooltools"}), on=["chrom", "start"])
        ok = np.isfinite(u["E1"]) & np.isfinite(u["E1_cooltools"])
        out["r_e1_uchrom_vs_cooltools"] = float(np.corrcoef(u.loc[ok, "E1"], u.loc[ok, "E1_cooltools"])[0, 1])
        out["max_abs_e1_diff"] = float((u.loc[ok, "E1"] - u.loc[ok, "E1_cooltools"]).abs().max())
        if mean.exists():
            mm = u.merge(pd.read_csv(mean, sep="\t")[["chrom", "start", "mean_comp"]], on=["chrom", "start"])
            out["max_abs_mean_score_diff_uchrom_vs_schicluster"] = float((mm["mean_scab"] - mm["mean_comp"]).abs().max())
    return out


def _phase(near, mit):
    g = np.full(len(near), "mid-S/G2", dtype=object)
    g[(near > 0.611) & (near <= 0.77)] = "early-S"
    g[near <= 0.611] = "G1"
    g[(mit >= 0.3) & (near <= 0.42)] = "post-M"
    g[(near > 0.42) & (-1.8 * mit - near + 1 < 0)] = "pre-M"
    return g


def _summary(mit, near, prefix=""):
    ph = pd.Series(_phase(near, mit)).value_counts()
    return {f"{prefix}median_f_mitotic": float(np.median(mit)), f"{prefix}median_f_near": float(np.median(near)),
            f"{prefix}n_cells_f_mitotic_ge_0.3": int((mit >= 0.3).sum()),
            **{f"{prefix}n_{p}": int(ph.get(p, 0)) for p in ("post-M", "G1", "early-S", "mid-S/G2", "pre-M")}}


def r6(a) -> dict:
    d = a.data3 / "nagano2017_hap_2i"
    t0 = time.perf_counter()
    cells = pd.read_csv(d / "cells.tsv", sep="\t")["cell"].astype(str).tolist()
    fends = pd.read_csv(d / "GATC.fends", sep="\t", usecols=["fend", "chr", "coord"], dtype={"chr": str})
    fend_ids = fends["fend"].to_numpy(np.int64)
    order = np.argsort(fend_ids)
    fid = fend_ids[order]
    fchr = ("chr" + fends["chr"].astype(str).str.replace("^chr", "", regex=True)).to_numpy()[order]
    fpos = fends["coord"].to_numpy(np.int64)[order]
    keep_chr = set(CHROMS6)
    rows = []
    n_edge = 0
    for cell in cells:
        adj = pd.read_csv(d / "adj" / f"{cell}.adj", sep="\t")
        f1, f2 = adj["fend1"].to_numpy(np.int64), adj["fend2"].to_numpy(np.int64)
        pairs = np.unique(np.stack([np.minimum(f1, f2), np.maximum(f1, f2)], axis=1), axis=0)
        k1, k2 = np.searchsorted(fid, pairs[:, 0]), np.searchsorted(fid, pairs[:, 1])
        k1, k2 = np.clip(k1, 0, len(fid) - 1), np.clip(k2, 0, len(fid) - 1)
        ok = (fid[k1] == pairs[:, 0]) & (fid[k2] == pairs[:, 1])
        c1, c2 = fchr[k1], fchr[k2]
        cis = ok & (c1 == c2) & np.isin(c1, list(keep_chr))
        dist = np.abs(fpos[k2][cis] - fpos[k1][cis])
        valid = dist > 2 ** 14.5
        near = (dist > 2 ** 14.5) & (dist <= 2 ** 21)
        mit = (dist > 2 ** 21) & (dist <= 2 ** 23.5)
        n_edge += int(np.isin(dist, [2 ** 21]).sum())
        rows.append((cell, int(valid.sum()), near.sum() / valid.sum(), mit.sum() / valid.sum()))
    t = pd.DataFrame(rows, columns=["cell", "n_valid", "f_near", "f_mitotic"])
    secs = time.perf_counter() - t0
    out = {"runtime_s": secs, "n_cells": len(t), "contacts_at_exactly_2^21": n_edge,
           **_summary(t["f_mitotic"].to_numpy(), t["f_near"].to_numpy())}
    feat = pd.read_csv(a.shared / "nagano2017_hap/GSE94489_haploids_features_table.txt.gz", sep="\t")
    feat["cell"] = feat["cell_nm"].str.replace("_", ".", regex=False)
    m = t.merge(feat[["cell", "f_near_band", "f_mitotic_band", "group"]], on="cell", how="left")
    dn = (m["f_near"] - m["f_near_band"]).abs()
    dm = (m["f_mitotic"] - m["f_mitotic_band"]).abs()
    out["published"] = {"n_cells_matched": int(m["f_near_band"].notna().sum()),
                        "max_abs_diff_near": float(dn.max()), "median_abs_diff_near": float(dn.median()),
                        "max_abs_diff_mitotic": float(dm.max()), "median_abs_diff_mitotic": float(dm.median()),
                        **_summary(m["f_mitotic_band"].to_numpy(float), m["f_near_band"].to_numpy(float),
                                   prefix="published_"),
                        "published_group_counts": {str(k): int(v) for k, v in m["group"].value_counts().items()}}
    ph_ours = _phase(t["f_near"].to_numpy(), t["f_mitotic"].to_numpy())
    ph_pub = _phase(m["f_near_band"].to_numpy(float), m["f_mitotic_band"].to_numpy(float))
    out["published"]["gate_agreement_with_published_values"] = float((ph_ours == ph_pub).mean())
    t.to_csv(a.work / "schic_R6_numpy_cells.tsv", sep="\t", index=False)
    u = a.work / "schic_R6_uchrom_cells.tsv"
    if u.exists():
        uu = pd.read_csv(u, sep="\t").merge(t, on="cell", suffixes=("_uchrom", "_numpy"))
        out["vs_uchrom"] = {"max_abs_diff_near": float((uu["f_near_uchrom"] - uu["f_near_numpy"]).abs().max()),
                            "max_abs_diff_mitotic": float((uu["f_mitotic_uchrom"] - uu["f_mitotic_numpy"]).abs().max()),
                            "n_cells": int(len(uu))}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--data3", type=Path, required=True)
    ap.add_argument("--shared", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", default="R5,R6")
    a = ap.parse_args()
    out = json.loads(a.out.read_text()) if a.out.exists() else {}
    for q in a.only.split(","):
        t0 = time.perf_counter()
        try:
            v = {"R5": r5, "R6": r6}[q](a)
        except Exception as e:
            import traceback
            traceback.print_exc()
            v = {"error": repr(e)}
        v["wall_s"] = time.perf_counter() - t0
        out[q] = v
        print(q, json.dumps(v, default=str)[:800], flush=True)
        a.out.write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
