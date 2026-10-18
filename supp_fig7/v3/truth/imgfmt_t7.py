"""Question T7 (family imgfmt): raw DNA seqFISH+ detections of Takei et al. 2021 (25-kb loci, replicate 1) assigned to
chromosome copies by a polymer-aware aligner (jie, Jia et al. 2022) and by plain baselines, scored against the
authors' homolog call (``labelID``).

    python imgfmt_t7.py make  --zip DNAseqFISH+.zip --fofct 4DNFIHF3JCBY.csv --out FOLDER [--chroms chr3,chr7,chr11]
    python imgfmt_t7.py run   --impl uchrom_default|uchrom_tuned|jie|kmeans2|dbscan|nnchain|all_one
                              --data FOLDER --out JSON [--jobs 8] [--param k=v ...]

``make`` writes ``detections.tsv`` (what an agent gets: spot_id, cell_id, chrom, locus, start, end, x_nm, y_nm, z_nm)
and ``labels.tsv`` (spot_id, label: the authors' homolog call, -1 = rejected).  ``run`` writes the assignment
(spot_id -> fibre; spots on no fibre left out) and its scores (:func:`score`, the definitions of
``uchrom.im.trace.score_alignment`` re-implemented with pandas / scikit-learn).
"""

from __future__ import annotations

import argparse
import io
import json
import random
import sys
import time
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

MEMBER = "DNAseqFISH+/DNAseqFISH+25kbloci-E14-replicate1.csv"
VOXEL_NM = (103.0, 103.0, 250.0)


def chrom_name(i: int) -> str:
    return "chrX" if int(i) == 20 else f"chr{int(i)}"


def read_fofct(path: Path) -> pd.DataFrame:
    """The 4DN FOF-CT core table: ``##columns=(...)`` names the columns of the rows after the ``#`` header lines."""
    cols, skip = None, 0
    with open(path) as fh:
        for line in fh:
            s = line.strip().strip('"')
            if not s.startswith("#"):
                break
            skip += 1
            if s.startswith("##columns="):
                cols = [c.strip() for c in s.split("=", 1)[1].split(")")[0].strip("(").split(",")]
    df = pd.read_csv(path, skiprows=skip, header=None)
    df = df.iloc[:, :len(cols)]
    df.columns = cols
    return df


def make_input(zpath: Path, fofct: Path, out: Path, chroms) -> dict:
    with zipfile.ZipFile(zpath) as zf, zf.open(MEMBER) as fh:
        raw = pd.read_csv(fh)
    raw["spot_id"] = np.arange(len(raw))
    raw["chrom"] = raw["chromID"].map(chrom_name)
    raw["locus"] = raw["regionID (hyb1-60)"].astype(int)
    raw["cell_id"] = raw["fov"].astype(str) + "_" + raw["cellID"].astype(str)
    # the mm10 locus of every (chromosome, regionID): the traced FOF-CT table of the same experiment holds the same
    # spots (coordinates = pixels x voxel size in um, 3 decimals), with Chrom_Start / Chrom_End
    f = read_fofct(fofct)
    key_f = pd.MultiIndex.from_arrays([f["X"].round(3), f["Y"].round(3), f["Z"].round(3)])
    if key_f.has_duplicates:
        raise ValueError("FOF-CT positions are not unique")
    idx = pd.Series(np.arange(len(f)), index=key_f)
    key_r = pd.MultiIndex.from_arrays([(raw["x"] * VOXEL_NM[0] / 1000).round(3), (raw["y"] * VOXEL_NM[1] / 1000).round(3),
                                       (raw["z"] * VOXEL_NM[2] / 1000).round(3)])
    row = idx.reindex(key_r).to_numpy()
    if np.isnan(row).any():
        raise ValueError(f"{int(np.isnan(row).sum())} detections without a FOF-CT spot")
    m = f.iloc[row.astype(np.int64)]
    same_cell = float(np.mean(m["Cell_ID"].astype(str).to_numpy() == raw["cell_id"].to_numpy()))
    same_label = float(np.mean(m["Trace_ID"].astype(str).str.rsplit("_", n=1).str[-1].to_numpy()
                               == raw["labelID"].astype(str).to_numpy()))
    raw["start"] = m["Chrom_Start"].to_numpy(np.int64)
    raw["end"] = m["Chrom_End"].to_numpy(np.int64)
    loci = raw[["chrom", "locus", "start", "end"]].drop_duplicates()
    if loci.duplicated(["chrom", "locus"]).any():
        raise ValueError("a (chromosome, regionID) maps to two loci")
    sel = raw[raw["chrom"].isin(chroms)].copy()
    sel["x_nm"] = sel["x"] * VOXEL_NM[0]
    sel["y_nm"] = sel["y"] * VOXEL_NM[1]
    sel["z_nm"] = sel["z"] * VOXEL_NM[2]
    out.mkdir(parents=True, exist_ok=True)
    det = sel[["spot_id", "cell_id", "chrom", "locus", "start", "end", "x_nm", "y_nm", "z_nm"]]
    det.to_csv(out / "detections.tsv", sep="\t", index=False, float_format="%.3f")
    sel[["spot_id", "labelID"]].rename(columns={"labelID": "label"}).to_csv(out / "labels.tsv", sep="\t",
                                                                           index=False)
    info = {"member": MEMBER, "chroms": list(chroms), "n_detections": int(len(det)),
            "n_cells": int(det["cell_id"].nunique()), "fofct_match": {"same_cell": same_cell,
                                                                     "label_equals_fofct_homolog": same_label},
            "label_counts": {str(k): int(v) for k, v in sel["labelID"].value_counts().sort_index().items()},
            "loci_per_chrom": {c: int(g["locus"].nunique()) for c, g in det.groupby("chrom")}}
    (out / "build_info.json").write_text(json.dumps(info, indent=1))
    return info


# ----------------------------------------------------------------------------------------------- scoring
def score(assign: pd.DataFrame, det: pd.DataFrame, labels: pd.DataFrame, groups=None) -> dict:
    """``assign``: spot_id, fibre (spots on no fibre absent).  Scored: spots with label >= 0 of the (cell, chrom)
    groups in ``groups`` (default: every group of ``det``).  assigned = share of them on a fibre; purity = per fibre
    (with scored spots), share of its majority label, mean over fibres; median ARI = per (cell, chrom) with more than
    one label, adjusted Rand index of fibre (unassigned = one extra group) vs label, median over groups."""
    from sklearn.metrics import adjusted_rand_score

    df = det[["spot_id", "cell_id", "chrom"]].merge(labels, on="spot_id")
    df = df[df["label"] >= 0]
    if groups is not None:
        g = pd.MultiIndex.from_frame(groups[["cell_id", "chrom"]].astype(str))
        df = df[pd.MultiIndex.from_frame(df[["cell_id", "chrom"]].astype(str)).isin(g)]
    fib = assign.set_index("spot_id")["fibre"].astype(str)
    df = df.assign(fibre=fib.reindex(df["spot_id"].to_numpy()).to_numpy())
    on = df[df["fibre"].notna()]
    purity = on.groupby("fibre")["label"].agg(lambda s: s.value_counts().iloc[0] / len(s))
    aris = []
    for _, g in df.groupby(["cell_id", "chrom"]):
        if g["label"].nunique() > 1:
            aris.append(adjusted_rand_score(g["label"], g["fibre"].fillna("none")))
    per_chrom = {}
    for c, g in df.groupby("chrom"):
        a = [adjusted_rand_score(h["label"], h["fibre"].fillna("none")) for _, h in g.groupby("cell_id")
             if h["label"].nunique() > 1]
        o = g[g["fibre"].notna()]
        per_chrom[c] = {"assigned": float(g["fibre"].notna().mean()),
                        "purity": float(o.groupby("fibre")["label"].agg(
                            lambda s: s.value_counts().iloc[0] / len(s)).mean()) if len(o) else float("nan"),
                        "median_ari": float(np.median(a)) if a else float("nan"), "n_groups_ari": len(a)}
    return {"assigned": float(df["fibre"].notna().mean()), "purity": float(purity.mean()),
            "median_ari": float(np.median(aris)) if aris else float("nan"), "mean_ari": float(np.mean(aris)),
            "n_scored_spots": int(len(df)), "n_fibres": int(on["fibre"].nunique()), "n_groups_ari": len(aris),
            "per_chrom": per_chrom}


# -------------------------------------------------------------------------------------------- U-Chrom
TUNED = dict(max_skip=12, gap_penalty=5.0, min_spots_per_fiber=10, max_fibers_per_chrom=4,
             mean_edge_weight_cutoff=8.0, sister_pair_radius_nm=400.0)


def run_uchrom(det: pd.DataFrame, tuned: bool, params: dict):
    from uchrom.im.trace import SpotAlignerParams, align_spots

    p = SpotAlignerParams(**{**(TUNED if tuned else {}), **params})
    df = det.rename(columns={"x_nm": "x", "y_nm": "y", "z_nm": "z"})[
        ["spot_id", "cell_id", "chrom", "start", "end", "x", "y", "z"]]
    cd = align_spots(df, params=p)
    sp = cd.spots
    return pd.DataFrame({"spot_id": sp["spot_id"].to_numpy(), "fibre": sp["trace_id"].astype(str).to_numpy()}), \
        {"uchrom_params": {k: (v if np.isfinite(v) else str(v)) if isinstance(v, float) else v
                           for k, v in p.__dict__.items()},
         "tau": {c: v.get("tau_nm_per_bp") for c, v in cd.uns.get("jie_polymer", {}).items()}}


# ------------------------------------------------------------------------------------- authors' jie
def _jie_one(args):
    import warnings

    from jie import aligner
    cell_pts, gene_dist, bin_size, kw, seed = args
    random.seed(seed)
    np.random.seed(seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fibres = aligner.find_all_chr(cell_pts, gene_dist, bin_size, **kw)
        except Exception as exc:          # recorded, counted as unassigned
            return [], f"{type(exc).__name__}: {exc}"
    return [f["spot_id"].to_numpy() for f in fibres], None


def fit_nm_per_bp(det: pd.DataFrame, l_p_bp: float = 150.0) -> dict:
    """nm per bp of jie's Gaussian chain from the data, without labels: the median distance from each detection
    to the nearest detection of the next locus of its (cell, chromosome) (median of |r| of a 3-D Gaussian =
    1.5382 sigma, <r^2> = 3 sigma^2 = 2 l_p L c^2 with L the locus spacing in bp)."""
    out = {}
    for chrom, g in det.groupby("chrom"):
        loci = g[["locus", "start"]].drop_duplicates().sort_values("locus")
        step = float(np.median(np.diff(loci["start"].to_numpy())))
        d = []
        for _, c in g.groupby("cell_id"):
            by = {k: v[["x_nm", "y_nm", "z_nm"]].to_numpy() for k, v in c.groupby("locus")}
            for k, a in by.items():
                b = by.get(k + 1)
                if b is None:
                    continue
                d.extend(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min(1))
        med = float(np.median(d))
        r2 = 3 * (med / 1.5382) ** 2
        out[chrom] = {"median_nn_nm": med, "step_bp": step, "nm_per_bp": float(np.sqrt(r2 / (2 * l_p_bp * step)))}
    return out


def run_jie(det: pd.DataFrame, jobs: int, params: dict):
    scale = float(params.pop("nm_scale", 1.0))
    fit = fit_nm_per_bp(det)
    pixel = 103.0
    tasks, keys = [], []
    for (cell, chrom), g in det.groupby(["cell_id", "chrom"], sort=True):
        loci = det.loc[det["chrom"] == chrom, ["locus", "start"]].drop_duplicates().sort_values("locus")
        hyb_of = {l: i for i, l in enumerate(loci["locus"])}
        gene_dist = (loci["start"] - loci["start"].iloc[0]).to_numpy().astype(float)
        bin_size = float(np.median(np.diff(loci["start"].to_numpy())))
        pts = pd.DataFrame({"x_hat": g["x_nm"].to_numpy() / pixel, "y_hat": g["y_nm"].to_numpy() / pixel,
                            "z_hat": g["z_nm"].to_numpy() / pixel, "sig_x": 1.0, "sig_y": 1.0, "sig_z": 1.0,
                            "hyb": g["locus"].map(hyb_of).astype(int).to_numpy(),
                            "spot_id": g["spot_id"].to_numpy()}).sort_values("hyb", kind="mergesort")
        pts = pts.reset_index(drop=True)
        kw = {"nm_per_bp": fit[chrom]["nm_per_bp"] * scale, "pixel_dist": pixel, **params}
        tasks.append((pts, gene_dist, bin_size, kw, len(tasks)))
        keys.append((cell, chrom))
    t0 = time.perf_counter()
    with ProcessPoolExecutor(jobs) as ex:
        res = list(ex.map(_jie_one, tasks, chunksize=4))
    rows, errors = [], []
    for (cell, chrom), (fibres, err) in zip(keys, res):
        if err:
            errors.append(f"{cell} {chrom}: {err}")
        for k, ids in enumerate(fibres):
            rows.append(pd.DataFrame({"spot_id": ids, "fibre": f"{cell}_{chrom}_{k}"}))
    assign = pd.concat(rows) if rows else pd.DataFrame(columns=["spot_id", "fibre"])
    return assign.drop_duplicates("spot_id"), {"nm_per_bp_fit": fit, "nm_scale": scale, "jie_kwargs": params,
                                               "n_errors": len(errors), "errors": errors[:5],
                                               "align_seconds": time.perf_counter() - t0}


# --------------------------------------------------------------------------------------------- baselines
def run_kmeans2(det: pd.DataFrame, params: dict):
    """Two clusters of every (cell, chromosome) by k-means on the 3-D positions (every detection assigned)."""
    from sklearn.cluster import KMeans
    rows = []
    for (cell, chrom), g in det.groupby(["cell_id", "chrom"]):
        xyz = g[["x_nm", "y_nm", "z_nm"]].to_numpy()
        k = min(int(params.get("k", 2)), len(g))
        lab = KMeans(k, n_init=10, random_state=0).fit_predict(xyz) if k > 1 else np.zeros(len(g), int)
        rows.append(pd.DataFrame({"spot_id": g["spot_id"].to_numpy(), "fibre": [f"{cell}_{chrom}_{c}" for c in lab]}))
    return pd.concat(rows), {"k": params.get("k", 2)}


def run_dbscan(det: pd.DataFrame, params: dict):
    from sklearn.cluster import DBSCAN
    eps, ms = float(params.get("eps_nm", 1000.0)), int(params.get("min_samples", 5))
    rows = []
    for (cell, chrom), g in det.groupby(["cell_id", "chrom"]):
        lab = DBSCAN(eps=eps, min_samples=ms).fit_predict(g[["x_nm", "y_nm", "z_nm"]].to_numpy())
        keep = lab >= 0
        rows.append(pd.DataFrame({"spot_id": g["spot_id"].to_numpy()[keep],
                                  "fibre": [f"{cell}_{chrom}_{c}" for c in lab[keep]]}))
    return pd.concat(rows), {"eps_nm": eps, "min_samples": ms}


def run_nnchain(det: pd.DataFrame, params: dict):
    """Greedy nearest-neighbour chaining along the genome: start a chain at the first unused detection of the
    lowest locus; repeatedly step to the nearest unused detection within the next ``max_skip + 1`` loci and
    ``max_step_nm``; chains with fewer than ``min_len`` spots are dropped; repeat until no chain can start."""
    max_skip, max_step, min_len = int(params.get("max_skip", 3)), float(params.get("max_step_nm", 1000.0)), \
        int(params.get("min_len", 5))
    rows = []
    for (cell, chrom), g in det.groupby(["cell_id", "chrom"]):
        g = g.sort_values("locus")
        loc = g["locus"].to_numpy()
        xyz = g[["x_nm", "y_nm", "z_nm"]].to_numpy()
        sid = g["spot_id"].to_numpy()
        used = np.zeros(len(g), bool)
        k = 0
        tried = np.zeros(len(g), bool)
        while True:
            free = np.flatnonzero(~used & ~tried)
            if not len(free):
                break
            cur = free[0]
            chain = [cur]
            used[cur] = True
            while True:
                cand = np.flatnonzero(~used & (loc > loc[cur]) & (loc <= loc[cur] + max_skip + 1))
                if not len(cand):
                    break
                d = np.linalg.norm(xyz[cand] - xyz[cur], axis=1)
                j = int(np.argmin(d))
                if d[j] > max_step:
                    break
                cur = cand[j]
                chain.append(cur)
                used[cur] = True
            if len(chain) >= min_len:
                rows.append(pd.DataFrame({"spot_id": sid[chain], "fibre": f"{cell}_{chrom}_{k}"}))
                k += 1
            else:
                used[chain] = False
                tried[chain[0]] = True
    return pd.concat(rows), {"max_skip": max_skip, "max_step_nm": max_step, "min_len": min_len}


def run_all_one(det: pd.DataFrame, params: dict):
    return pd.DataFrame({"spot_id": det["spot_id"].to_numpy(),
                         "fibre": (det["cell_id"] + "_" + det["chrom"]).to_numpy()}), {}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["make", "run"])
    ap.add_argument("--zip", type=Path)
    ap.add_argument("--fofct", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--chroms", default="chr3,chr7,chr11")
    ap.add_argument("--impl")
    ap.add_argument("--data", type=Path)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--cells", type=int, default=None, help="first N cells (sorted ids) only")
    ap.add_argument("--param", action="append", default=[])
    a = ap.parse_args(argv)
    if a.cmd == "make":
        print(json.dumps(make_input(a.zip, a.fofct, a.out, a.chroms.split(",")), indent=1))
        return
    det = pd.read_csv(a.data / "detections.tsv", sep="\t", dtype={"cell_id": str})
    labels = pd.read_csv(a.data / "labels.tsv", sep="\t")
    chroms = a.chroms.split(",") if a.chroms else None
    if chroms:
        det = det[det["chrom"].isin(chroms)]
    if a.cells:
        keep = sorted(det["cell_id"].unique())[:a.cells]
        det = det[det["cell_id"].isin(keep)]
    params = {}
    for kv in a.param:
        k, v = kv.split("=", 1)
        try:
            params[k] = json.loads(v)
        except ValueError:
            params[k] = v
    t0 = time.perf_counter()
    if a.impl in ("uchrom_default", "uchrom_tuned"):
        assign, info = run_uchrom(det, a.impl == "uchrom_tuned", params)
    elif a.impl == "jie":
        assign, info = run_jie(det, a.jobs, params)
    else:
        assign, info = {"kmeans2": run_kmeans2, "dbscan": run_dbscan, "nnchain": run_nnchain,
                        "all_one": run_all_one}[a.impl](det, params)
    secs = time.perf_counter() - t0
    groups = det[["cell_id", "chrom"]].drop_duplicates()
    sc = score(assign, det, labels, groups)
    res = {"impl": a.impl, "params": params, "seconds": secs, "n_detections": int(len(det)),
           "n_cells": int(det["cell_id"].nunique()), "chroms": sorted(det["chrom"].unique()), **sc, "info": info}
    print(json.dumps({k: v for k, v in res.items() if k != "per_chrom"}, indent=1, default=str), flush=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1, default=str))
    assign.to_csv(a.out.with_suffix(".assign.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
