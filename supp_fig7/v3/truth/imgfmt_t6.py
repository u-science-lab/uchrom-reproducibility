"""Question T6 (family imgfmt): SnapFISH-IMPUTE on a fixed hold-out of observed loci of the Huang et al. 2021 Sox2
tracing, three ways — U-Chrom (``uchrom.im.impute.impute_snapfish``), the authors' package called function by function
(single process, as ``sfimpute.impute.generate_single_reg_target`` + ``recover_3d_coor_single_chr``) and the authors'
MPI script (``run_impute.py`` with ``mpiexec``, which down-scales the distance matrices before the similarity).

    python imgfmt_t6.py make  --src DIR_WITH_HUANG_TABLES --out DATA_FOLDER [--region 129] [--n-hap N] [--frac 0.1]
    python imgfmt_t6.py run   --impl uchrom|sfi_func|sfi_mpi|linear --data DATA_FOLDER --out JSON [--jobs 8]
                              [--run-impute PATH/run_impute.py] [--work DIR]

Each ``run`` writes the imputed coordinates of the held-out loci and the error summary (RMSE, median, mean, not
imputed) to its JSON; the inputs are the data folder's files only.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

COOR = "sox2_coor.tsv"          # region, haploid, pos, x, y, z (nm); NaN = undetected; the held-out loci keep their xyz
ANN = "sox2_ann.tsv"            # region, pos, start, end
HOLD = "holdout.tsv"            # region, haploid, pos: the observed loci to hide


def make_input(src: Path, out: Path, region: str = "129", n_hap: int | None = None, frac: float = 0.1,
               seed: int = 0) -> dict:
    """Data folder of T6 from the registered tables (``huang2021_sox2``): one allele (``region``), its first ``n_hap``
    haploids in file order (all when None), and ``frac`` of its observed loci drawn with
    ``numpy.random.default_rng(seed).choice`` (without replacement) from the observed rows in file order."""
    coor = pd.read_csv(src / "mESC_Sox2_coor_wnan.txt", sep="\t", dtype={"region": str, "haploid": str})
    ann = pd.read_csv(src / "mESC_Sox2_ann.txt", sep="\t", dtype={"region": str})
    c = coor[coor["region"] == region].reset_index(drop=True)
    haps = pd.unique(c["haploid"])
    if n_hap:
        haps = haps[:n_hap]
        c = c[c["haploid"].isin(haps)].reset_index(drop=True)
    a = ann[ann["region"] == region].reset_index(drop=True)
    obs = c[c[["x", "y", "z"]].notna().all(axis=1)]
    n = int(np.floor(len(obs) * frac + 1e-9))
    pick = np.sort(np.random.default_rng(seed).choice(len(obs), size=n, replace=False))
    hold = obs.iloc[pick][["region", "haploid", "pos"]].reset_index(drop=True)
    left = (obs.assign(h=False).set_index(["haploid", "pos"])["h"])
    left.loc[list(zip(hold["haploid"], hold["pos"]))] = True
    kept = (~left).groupby(level=0).sum()
    out.mkdir(parents=True, exist_ok=True)
    c.to_csv(out / COOR, sep="\t", index=False)
    a.to_csv(out / ANN, sep="\t", index=False)
    hold.to_csv(out / HOLD, sep="\t", index=False)
    info = {"region": region, "n_haploids": int(len(haps)), "n_loci": int(len(a)), "n_rows": int(len(c)),
            "n_observed": int(len(obs)), "n_holdout": int(n), "frac": frac, "seed": seed,
            "min_observed_after_holdout": int(kept.min()),
            "loci_never_observed": [int(p) for p in sorted(set(a["pos"]) - set(obs["pos"]))]}
    (out / "build_info.json").write_text(json.dumps(info, indent=1))
    return info


def load(data: Path):
    coor = pd.read_csv(data / COOR, sep="\t", dtype={"region": str, "haploid": str})
    ann = pd.read_csv(data / ANN, sep="\t", dtype={"region": str})
    hold = pd.read_csv(data / HOLD, sep="\t", dtype={"region": str, "haploid": str})
    key = pd.MultiIndex.from_frame(hold[["haploid", "pos"]])
    truth = coor.set_index(["haploid", "pos"]).loc[key, ["x", "y", "z"]].to_numpy(float)
    masked = coor.copy()
    m = pd.MultiIndex.from_frame(coor[["haploid", "pos"]]).isin(key)
    masked.loc[m, ["x", "y", "z"]] = np.nan
    return coor, ann, hold, masked, truth


def summary(est: np.ndarray, truth: np.ndarray) -> dict:
    err = np.linalg.norm(est - truth, axis=1)
    ok = np.isfinite(err)
    return {"rmse_nm": float(np.sqrt(np.mean(err[ok] ** 2))), "median_error_nm": float(np.median(err[ok])),
            "mean_error_nm": float(np.mean(err[ok])), "n_holdout": int(len(err)), "n_not_imputed": int((~ok).sum())}


# ------------------------------------------------------------------------------------------------------ U-Chrom
def run_uchrom(masked: pd.DataFrame, ann: pd.DataFrame, hold: pd.DataFrame, jobs: int):
    from uchrom.im.impute import impute_snapfish
    from uchrom.io import read_snapfish

    cd = read_snapfish(masked, ann)             # each region a chromosome of its own (SnapFISH's view), nm
    t0 = time.perf_counter()
    out = impute_snapfish(cd, n_jobs=jobs, verbose=True)
    secs = time.perf_counter() - t0
    sp = out.spots_with_loci()
    region = str(ann["region"].iloc[0])
    pos_of = dict(zip(ann["start"].astype(np.int64), ann["pos"].astype(int)))
    tab = pd.DataFrame({"trace_id": sp["trace_id"].astype(str).to_numpy(),
                        "pos": [pos_of[int(s)] for s in sp["start"].to_numpy()],
                        "x": out.coords[:, 0], "y": out.coords[:, 1], "z": out.coords[:, 2]})
    tab = tab.drop_duplicates(["trace_id", "pos"]).set_index(["trace_id", "pos"])
    key = pd.MultiIndex.from_arrays([hold["haploid"].astype(str) + "_" + region, hold["pos"].astype(int)])
    return tab.reindex(key)[["x", "y", "z"]].to_numpy(float), secs


# ------------------------------------------------------------------------------- authors' package, one process
def _recover_task(args):
    import warnings

    from sfimpute.impute import recover_3d_coor_single_chr
    with warnings.catch_warnings(), np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        warnings.simplefilter("ignore")
        return recover_3d_coor_single_chr(*args)


def run_sfi_func(masked: pd.DataFrame, ann: pd.DataFrame, hold: pd.DataFrame, jobs: int):
    """The authors' functions in the order of their single-process path: linear interpolation (initial values),
    pairwise distances normalised by genomic distance (Box-Cox per separation), similarity-weighted filling until
    nothing changes (``impute_pdist_one_iter``, as ``generate_single_reg_target``), then L-BFGS-B recovery of every
    haploid (``recover_3d_coor_single_chr``, as ``parallel_single_recover`` but over a process pool)."""
    from concurrent.futures import ProcessPoolExecutor

    from sfimpute.impute import (count_na_in_pdist, impute_pdist_one_iter, interpolate_coors_scipy,
                                 normalize_pdist_by1d, to_dist_df)

    t0 = time.perf_counter()
    coor = masked.copy()
    lnr = interpolate_coors_scipy(ann, coor, "linear")
    np.seterr(divide="ignore", invalid="ignore")
    norm = normalize_pdist_by1d(to_dist_df(coor, ann))
    pdist_arr = norm.to_darr().T
    prev_na, curr_na = None, count_na_in_pdist(pdist_arr)
    n_iter = 0
    while prev_na != curr_na and curr_na != 0:
        P = impute_pdist_one_iter(pdist_arr)
        pdist_arr = np.where(np.isnan(pdist_arr), P, pdist_arr)
        prev_na, curr_na = curr_na, count_na_in_pdist(pdist_arr)
        n_iter += 1
        print(f"iteration {n_iter}: {curr_na} NaN", flush=True)
    norm[norm.vcol()] = pdist_arr.T
    trf = norm[["lmbda", "mu", "var"]].values
    haps = pd.unique(coor["haploid"])
    tasks = []
    for h in haps:
        raw = coor[coor["haploid"] == h].sort_values("pos")[["x", "y", "z"]].to_numpy(float)
        ln = lnr[lnr["haploid"] == h].sort_values("pos")[["x", "y", "z"]].to_numpy(float)
        tasks.append((raw, ln, norm[str(h)].to_numpy(float), trf))
    with ProcessPoolExecutor(jobs) as ex:
        rec = list(ex.map(_recover_task, tasks, chunksize=8))
    secs = time.perf_counter() - t0
    pos = np.sort(ann["pos"].to_numpy())
    tab = pd.concat([pd.DataFrame({"haploid": h, "pos": pos, "x": r[:, 0], "y": r[:, 1], "z": r[:, 2]})
                     for h, r in zip(haps, rec)]).set_index(["haploid", "pos"])
    key = pd.MultiIndex.from_frame(hold[["haploid", "pos"]])
    return tab.reindex(key)[["x", "y", "z"]].to_numpy(float), secs, n_iter


# ------------------------------------------------------------------------------------- authors' MPI script
def run_sfi_mpi(masked: pd.DataFrame, ann: pd.DataFrame, hold: pd.DataFrame, jobs: int, run_impute: Path,
                work: Path):
    work.mkdir(parents=True, exist_ok=True)
    cpath, apath = work / "masked_coor.txt", work / "ann.txt"
    masked.to_csv(cpath, sep="\t", index=False)
    ann.to_csv(apath, sep="\t", index=False)
    outdir = work / "out"
    venv_bin = Path(sys.executable).parent
    cmd = [str(venv_bin / "mpiexec"), "-n", str(jobs), sys.executable, str(run_impute), "-o", str(outdir),
           "-d", str(cpath), "-a", str(apath), "-s", "t6"]
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "OMP_NUM_THREADS": "1"})
    secs = time.perf_counter() - t0
    print(r.stdout[-3000:], r.stderr[-3000:], flush=True)
    if r.returncode:
        raise RuntimeError(f"run_impute.py failed ({r.returncode})")
    rec = pd.read_csv(outdir / "recover_coor_t6.txt", sep="\t", dtype={"region": str, "haploid": str})
    tab = rec.set_index(["haploid", "pos"])
    key = pd.MultiIndex.from_frame(hold[["haploid", "pos"]])
    return tab.reindex(key)[["x", "y", "z"]].to_numpy(float), secs


def run_linear(masked, ann, hold):
    """Baseline: linear interpolation along the genome (the authors' ``interpolate_coors_scipy`` definition)."""
    from scipy import interpolate
    t0 = time.perf_counter()
    rows = []
    x = ann.sort_values("pos")["start"].to_numpy(float)
    pos = ann.sort_values("pos")["pos"].to_numpy()
    for h, df in masked.groupby("haploid", sort=False):
        y = df.sort_values("pos")[["x", "y", "z"]].to_numpy(float)
        ok = ~np.isnan(y[:, 0])
        f = interpolate.interp1d(x[ok], y[ok], axis=0, kind="linear", fill_value="extrapolate")
        v = f(x)
        rows.append(pd.DataFrame({"haploid": h, "pos": pos, "x": v[:, 0], "y": v[:, 1], "z": v[:, 2]}))
    tab = pd.concat(rows).set_index(["haploid", "pos"])
    key = pd.MultiIndex.from_frame(hold[["haploid", "pos"]])
    return tab.reindex(key)[["x", "y", "z"]].to_numpy(float), time.perf_counter() - t0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["make", "run"])
    ap.add_argument("--src", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--region", default="129")
    ap.add_argument("--n-hap", type=int, default=None)
    ap.add_argument("--frac", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--impl", choices=["uchrom", "sfi_func", "sfi_mpi", "linear"])
    ap.add_argument("--data", type=Path)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--run-impute", type=Path)
    ap.add_argument("--work", type=Path)
    a = ap.parse_args(argv)
    if a.cmd == "make":
        print(json.dumps(make_input(a.src, a.out, a.region, a.n_hap, a.frac, a.seed), indent=1))
        return
    coor, ann, hold, masked, truth = load(a.data)
    extra = {}
    if a.impl == "uchrom":
        est, secs = run_uchrom(masked, ann, hold, a.jobs)
    elif a.impl == "sfi_func":
        est, secs, n_iter = run_sfi_func(masked, ann, hold, a.jobs)
        extra["n_iterations"] = n_iter
    elif a.impl == "sfi_mpi":
        est, secs = run_sfi_mpi(masked, ann, hold, a.jobs, a.run_impute, a.work)
    else:
        est, secs = run_linear(masked, ann, hold)
    res = {"impl": a.impl, "data": str(a.data), "jobs": a.jobs, "seconds": secs, **summary(est, truth), **extra,
           "python": sys.version.split()[0]}
    print(json.dumps(res, indent=1), flush=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1))
    np.save(a.out.with_suffix(".npy"), est)


if __name__ == "__main__":
    main()
