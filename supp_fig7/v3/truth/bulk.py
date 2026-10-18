"""Reference answers of the version-3 questions of family ``bulk`` (v3/questions/bulk.yaml): H1-H4, XM1, P1, P2.

Every question is computed twice:

* ``--impl uchrom``: with U-Chrom (the version the agents get, ``venv-uchrom``; functions named as ``truth:`` in the
  YAML);
* ``--impl independent``: without U-Chrom (functions ``<truth>_independent``): the authors' or Open2C's code where
  it exists (H1 ``mustache-hic`` 1.3.3, its command line; H2 cooltools 0.7.1 ``dots``; H3 cooltools ``pileup``;
  H4 scikit-learn's metric MDS) — run in ``venv-hicref`` — and numpy written from the papers for the short
  definitions (XM1; P1 and P2 run in ``venv-uchrom`` but import only ``chromdata``, the store's reader, never
  ``uchrom``: checked at the end);
* ``--impl finalize``: the graded values from both (the H1 band, the floors of H2 / H4, the Monte Carlo
  tolerances of P1), into ``results/bulk.json``.

    python bulk.py --impl uchrom      --data /data --data-v3 /data_v3 --dl /dl --out results/bulk.json [--only H1,..]
    python bulk.py --impl independent --data /data --data-v3 /data_v3 --dl /dl --out results/bulk.json [--only ..]
    python bulk.py --impl finalize    --out results/bulk.json

Results: ``{id: {"values": ..., "verify": {"uchrom": ..., "<independent>": ...}, "truth": fn, "runtime_s": ...}}``
(the layout of ``answers_v2.json``).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

#: data keys of the family -> (root, relative path); "v2" = the version-2 data root, "v3" = data_v3
DATA = {
    "rao_chr21_10kb": ("v3", "rao2014_imr90_chr21_10kb"),
    "su2020_tsv": ("v3", "su2020_chr21_tsv"),
    "takei_rna": ("v3", "takei2025_fov0_rna/takei2025_fov0_rna_counts.tsv"),
    "melanoma": ("v2", "hicrna_human_melanoma/hicrna_human_melanoma.chromdata.zarr"),
    "takei": ("v2", "takei2025_fov0/takei2025_fov0.chromdata.zarr"),
}
COOL = "IMR90_chr21_10kb_hg19.cool"
LOOPS = "GSE63525_IMR90_HiCCUPS_looplist.txt.gz"
RES = 10_000
TOL_BP = 25_000                    # H2: anchor midpoints within 25 kb
SECTION = "human_melanoma_r1"      # P1
N_RANDOM = 1000                    # P1: randomisations per run
MC_SEEDS = 10                      # P1: runs of N_RANDOM with different seeds (the Monte Carlo spread)
MIN_SPOTS = 50                     # P2: traces with at least this many spots

QUESTIONS = {  # id -> truth function (U-Chrom); the independent one is <name>_independent
    "H1": "mustache_loops_chr21", "H2": "hiccups_vs_published_chr21", "H3": "apa_published_chr21",
    "H4": "mds_vs_imaging_chr21", "XM1": "imaging_frequency_vs_hic_chr21", "P1": "msr_copy_number_share",
    "P2": "partial_spearman_rna_rg",
}


class Roots:
    def __init__(self, data: Path, data_v3: Path, dl: Path):
        self.data, self.v3, self.dl = Path(data), Path(data_v3), Path(dl)

    def __call__(self, key: str) -> Path:
        root, rel = DATA[key]
        return (self.data if root == "v2" else self.v3) / rel


# ================================================================================================ shared helpers
def published_loops(path: Path, chrom: str = "21") -> pd.DataFrame:
    """The authors' HiCCUPS loops of one chromosome (GSE63525 columns x1, x2, y1, y2), named as the map (chr21)."""
    t = pd.read_csv(path, sep="\t")
    t = t[t["chr1"].astype(str) == chrom]
    return pd.DataFrame({"chrom1": f"chr{chrom}", "start1": t["x1"].to_numpy(np.int64), "end1": t["x2"].to_numpy(np.int64),
                         "chrom2": f"chr{chrom}", "start2": t["y1"].to_numpy(np.int64),
                         "end2": t["y2"].to_numpy(np.int64)}).reset_index(drop=True)


def anchor_mids(t: pd.DataFrame, res: int = RES):
    e1 = t["end1"] if "end1" in t else t["start1"] + res
    e2 = t["end2"] if "end2" in t else t["start2"] + res
    return ((t["start1"].to_numpy(np.int64) + np.asarray(e1, np.int64)) / 2,
            (t["start2"].to_numpy(np.int64) + np.asarray(e2, np.int64)) / 2)


def matched(a: pd.DataFrame, b: pd.DataFrame, tol: float = TOL_BP) -> np.ndarray:
    """For every loop of ``a``: is there a loop of ``b`` with both anchor midpoints within ``tol`` (Chebyshev)."""
    if not len(a) or not len(b):
        return np.zeros(len(a), bool)
    a1, a2 = anchor_mids(a)
    b1, b2 = anchor_mids(b)
    d = np.maximum(np.abs(a1[:, None] - b1[None, :]), np.abs(a2[:, None] - b2[None, :]))
    return (d <= tol).any(1)


def vs_published(calls: pd.DataFrame, pub: pd.DataFrame) -> dict:
    rec, prec = matched(pub, calls), matched(calls, pub)
    return {"n_loops": int(len(calls)), "recall": float(rec.mean()), "precision": float(prec.mean()) if len(calls)
            else float("nan"), "n_published": int(len(pub)), "n_published_recovered": int(rec.sum())}


def apa_from_window(M: np.ndarray, corner: int = 6) -> dict:
    """Centre over the mean of the lower-left ``corner x corner`` square (rows = first anchor + 5..+10, columns =
    second anchor -10..-5), NaN left out."""
    L = M.shape[0]
    f = L // 2
    ll = M[L - corner:, :corner]
    return {"apa": float(M[f, f] / np.nanmean(ll)), "centre": float(M[f, f]), "lower_left_mean": float(np.nanmean(ll))}


def su_tables(folder: Path):
    """Su 2020 chr21 without U-Chrom: (xyz (copies, loci, 3) in nm, locus labels, Hi-C counts (loci x loci))."""
    hic = pd.read_csv(folder / "Hi-C_contacts_chromosome21.tsv", sep="\t", index_col=0)
    labels = [str(x) for x in hic.index]
    t = pd.read_csv(folder / "chromosome21.tsv", sep="\t",
                    usecols=["Z(nm)", "X(nm)", "Y(nm)", "Genomic coordinate", "Chromosome copy number"])
    pos = {lab: i for i, lab in enumerate(labels)}
    li = t["Genomic coordinate"].map(pos)
    if li.isna().any():
        raise ValueError("loci of chromosome21.tsv missing from the Hi-C table")
    copies, ci = np.unique(t["Chromosome copy number"].to_numpy(), return_inverse=True)
    xyz = np.full((len(copies), len(labels), 3), np.nan)
    xyz[ci, li.to_numpy(int)] = t[["X(nm)", "Y(nm)", "Z(nm)"]].to_numpy(float)
    return xyz, labels, hic.to_numpy(float)


def label_mids(labels) -> np.ndarray:
    se = np.array([[int(v) for v in lab.split(":")[1].split("-")] for lab in labels], float)
    return se.mean(axis=1)


def per_pair_imaging(xyz: np.ndarray, threshold_nm: float = 500.0):
    """Median distance (nm) and contact frequency (< threshold) of every locus pair over the copies where both are
    detected; (n, n) arrays, upper triangle filled (i < j)."""
    n = xyz.shape[1]
    med = np.full((n, n), np.nan)
    freq = np.full((n, n), np.nan)
    nobs = np.zeros((n, n), np.int64)
    for i in range(n - 1):
        d = np.sqrt(((xyz[:, i:i + 1, :] - xyz[:, i + 1:, :]) ** 2).sum(-1))          # copies x (n - i - 1)
        ok = np.isfinite(d)
        k = ok.sum(0)
        nobs[i, i + 1:] = k
        with np.errstate(invalid="ignore", divide="ignore"):
            med[i, i + 1:] = np.nanmedian(np.where(ok, d, np.nan), axis=0)
            freq[i, i + 1:] = np.where(k > 0, ((d < threshold_nm) & ok).sum(0) / np.maximum(k, 1), np.nan)
    return med, freq, nobs


def spearman(a, b) -> float:
    from scipy import stats

    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return float(stats.spearmanr(a[ok], b[ok]).statistic)


def decay_prior_distances(C: np.ndarray, weight: float = 0.05, alpha: float = 4.0) -> np.ndarray:
    """H4 as the question states it: c' = (1 - w) c + w E(|i - j|), E = mean of the non-zero counts of each index
    separation; d = c'^(-1/alpha); diagonal 0."""
    n = len(C)
    E = np.zeros(n)
    for k in range(1, n):
        dg = np.diagonal(C, k)
        nz = dg[dg > 0]
        E[k] = nz.mean() if len(nz) else 0.0
    sep = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    Cp = (1 - weight) * C + weight * E[sep]
    with np.errstate(divide="ignore"):
        D = np.where(Cp > 0, Cp ** (-1.0 / alpha), 0.0)
    np.fill_diagonal(D, 0.0)
    return D


def ranks_partial(x, y, groups) -> dict:
    """Spearman and partial Spearman given groups: ranks (average ties), minus their group mean, Pearson."""
    x, y = pd.Series(np.asarray(x, float)), pd.Series(np.asarray(y, float))
    g = pd.Series(np.asarray(groups)).astype(str)
    rx, ry = x.rank(method="average"), y.rank(method="average")
    ex = rx - rx.groupby(g).transform("mean")
    ey = ry - ry.groupby(g).transform("mean")
    return {"spearman": float(np.corrcoef(rx, ry)[0, 1]), "partial_spearman": float(np.corrcoef(ex, ey)[0, 1]),
            "n_cells": int(len(x)), "n_types": int(g.nunique())}


# ================================================================================================ H1 Mustache
def mustache_loops_chr21(roots: Roots) -> dict:
    """Mustache, U-Chrom's re-implementation, on the balanced 10-kb map: defaults of the question."""
    from uchrom.strc.loop.mustache import MustacheParams, call_loops_mustache

    d = roots("rao_chr21_10kb")
    p = MustacheParams(sigma0=1.6, octaves=2, scales=10, fdr=0.1, sparsity=0.88, max_dist=2_000_000)
    t = time.perf_counter()
    calls = call_loops_mustache(str(d / COOL), chrom="chr21", params=p, key_added=None)
    sec = time.perf_counter() - t
    out = {"n_loops": int(len(calls)), "seconds": sec, **{f"{k}_25kb": v for k, v in
                                                         vs_published(calls, published_loops(d / LOOPS)).items()
                                                         if k in ("recall", "precision")}}
    calls[["chrom1", "start1", "start2", "fdr"]].to_csv(roots.v3 / "bulk_work" / "H1_uchrom_calls.tsv", sep="\t",
                                                        index=False)
    # what an agent gets with the reference command line's own default FDR (-pt 0.2)
    out["n_loops_fdr0.2"] = int(len(call_loops_mustache(str(d / COOL), chrom="chr21", key_added=None,
                                                         params=MustacheParams(fdr=0.2, max_dist=2_000_000))))
    return out


def mustache_loops_chr21_independent(roots: Roots) -> dict:
    """mustache-hic 1.3.3, its command line, on the same .cool (balanced: the stored weights)."""
    import mustache

    d = roots("rao_chr21_10kb")
    exe = str(Path(sys.executable).parent / "mustache")
    work = roots.v3 / "bulk_work"
    out = {"implementation": f"mustache-hic {getattr(mustache, '__version__', '1.3.3')} (command line)"}
    for pt, tag in ((0.1, ""), (0.2, "_fdr0.2")):
        dst = work / f"H1_mustache_pt{pt}.tsv"
        cmd = [exe, "-f", str(d / COOL), "-r", str(RES), "-ch", "chr21", "-pt", str(pt), "-st", "0.88", "-sz", "1.6",
               "-oc", "2", "-i", "10", "-d", "2000000", "-p", "8", "-o", str(dst)]
        t = time.perf_counter()
        r = subprocess.run(cmd, capture_output=True, text=True)
        sec = time.perf_counter() - t
        if r.returncode:
            raise RuntimeError(r.stderr[-3000:])
        calls = pd.read_csv(dst, sep="\t")
        uniq = calls.drop_duplicates(["BIN1_START", "BIN2_START"])
        if tag:
            out["n_loops" + tag] = int(len(uniq))
            continue
        tab = pd.DataFrame({"chrom1": "chr21", "start1": uniq["BIN1_START"], "end1": uniq["BIN1_END"],
                            "start2": uniq["BIN2_START"], "end2": uniq["BIN2_END"]})
        sc = vs_published(tab, published_loops(d / LOOPS))
        out.update({"n_loops": int(len(uniq)), "n_rows_written": int(len(calls)), "seconds": sec,
                    "command": " ".join(cmd[1:]).replace(str(d) + "/", ""),
                    "recall_25kb": sc["recall"], "precision_25kb": sc["precision"]})
        # agreement with U-Chrom's calls (summit bins within 0 / 2 bins)
        u = work / "H1_uchrom_calls.tsv"
        if u.exists():
            ours = pd.read_csv(u, sep="\t")
            A = np.c_[uniq["BIN1_START"], uniq["BIN2_START"]] // RES
            B = np.c_[ours["start1"], ours["start2"]] // RES
            dd = np.maximum(np.abs(A[:, None, 0] - B[None, :, 0]), np.abs(A[:, None, 1] - B[None, :, 1]))
            out.update({"ref_in_uchrom_exact": float((dd.min(1) == 0).mean()),
                        "ref_in_uchrom_2bins": float((dd.min(1) <= 2).mean()),
                        "uchrom_in_ref_2bins": float((dd.min(0) <= 2).mean())})
    return out


# ================================================================================================ H2 HiCCUPS
def _empty_cd():
    from chromdata import ChromData

    return ChromData(np.zeros((0, 3)), pd.DataFrame({"chrom": pd.Series(dtype=str), "start": pd.Series(dtype=np.int64),
                                                     "end": pd.Series(dtype=np.int64),
                                                     "trace_id": pd.Series(dtype=str)}))


def hiccups_vs_published_chr21(roots: Roots) -> dict:
    """HiCCUPS (U-Chrom; cooltools' formulation, HiCCUPS' clustering) on the 10-kb map, loops <= 10 Mb."""
    from uchrom.strc.loop.hiccups import HiccupsParams, call_loops_hiccups

    d = roots("rao_chr21_10kb")
    pub = published_loops(d / LOOPS)
    out = {}
    for md, tag in ((10_000_000, ""), (2_000_000, "_2mb")):
        p = HiccupsParams(p=2, w=5, max_dist=md, fdr=0.1, cluster_radius=20_000)
        t = time.perf_counter()
        loops = call_loops_hiccups(_empty_cd(), contacts=str(d / COOL), chrom="chr21", params=p, key_added=None)
        sec = time.perf_counter() - t
        sc = vs_published(loops, pub)
        out.update({f"{k}{tag}": v for k, v in sc.items()})
        out[f"seconds{tag}"] = sec
    return out


def hiccups_vs_published_chr21_independent(roots: Roots) -> dict:
    """cooltools 0.7.1 ``dots`` (Open2C 2024) with its defaults at 10 kb (w = 5, p = 2, FDR 0.1 per λ-chunk,
    20-kb clustering, HiCCUPS filter), expected_cis with two diagonals ignored."""
    import logging
    import warnings

    import cooler
    import cooltools

    logging.disable(logging.WARNING)
    warnings.filterwarnings("ignore")
    d = roots("rao_chr21_10kb")
    pub = published_loops(d / LOOPS)
    clr = cooler.Cooler(str(d / COOL))
    view = pd.DataFrame({"chrom": ["chr21"], "start": [0], "end": [int(clr.chromsizes["chr21"])], "name": ["chr21"]})
    out = {"implementation": f"cooltools {cooltools.__version__} dots"}
    for md, tag in ((10_000_000, ""), (2_000_000, "_2mb")):
        t = time.perf_counter()
        exp = cooltools.expected_cis(clr, view_df=view, ignore_diags=2, nproc=1)
        dots = cooltools.dots(clr, expected=exp, view_df=view, max_loci_separation=md, lambda_bin_fdr=0.1,
                              clustering_radius=20_000, nproc=8)
        sec = time.perf_counter() - t
        sc = vs_published(dots, pub)
        out.update({f"{k}{tag}": v for k, v in sc.items()})
        out[f"seconds{tag}"] = sec
    return out


# ================================================================================================ H3 APA
def _apa_loops(pub: pd.DataFrame) -> pd.DataFrame:
    """The question's loops: anchor bins of the midpoints, j - i >= 22."""
    m1, m2 = anchor_mids(pub)
    i, j = (m1 // RES).astype(np.int64), (m2 // RES).astype(np.int64)
    return pub[(j - i) >= 22].reset_index(drop=True)


def apa_published_chr21(roots: Roots) -> dict:
    """``uchrom.fea.pileup`` (APA) with the question's definition: mean of the O/E windows (aggregate="mean")."""
    from uchrom.fea import pileup

    d = roots("rao_chr21_10kb")
    pub = published_loops(d / LOOPS)
    t = time.perf_counter()
    r = pileup(str(d / COOL), pub, flank=100_000, expected=True, aggregate="mean", corner=6, min_separation=22,
               ignore_diags=2, key_added=None)
    sec = time.perf_counter() - t
    rs = pileup(str(d / COOL), pub, flank=100_000, expected=True, aggregate="sum", corner=6, min_separation=22,
                ignore_diags=2, key_added=None)
    return {"n_loops": int(r["n"]), "apa": float(r["apa"]), "centre": float(r["center"]),
            "lower_left_mean": float(r["lower_left"]), "seconds": sec, "n_published": int(len(pub)),
            "n_loops_question_rule": int(len(_apa_loops(pub))), "apa_pooled_sum_aggregate": float(rs["apa"])}


def apa_published_chr21_independent(roots: Roots) -> dict:
    """cooltools 0.7.1 ``pileup`` (O/E snippets, expected_cis 'balanced.avg'), mean over loops ignoring NaN."""
    import logging
    import warnings

    import cooler
    import cooltools

    logging.disable(logging.WARNING)
    warnings.filterwarnings("ignore")
    d = roots("rao_chr21_10kb")
    loops = _apa_loops(published_loops(d / LOOPS))
    clr = cooler.Cooler(str(d / COOL))
    view = pd.DataFrame({"chrom": ["chr21"], "start": [0], "end": [int(clr.chromsizes["chr21"])], "name": ["chr21"]})
    t = time.perf_counter()
    exp = cooltools.expected_cis(clr, view_df=view, ignore_diags=2, nproc=1)
    stack = np.asarray(cooltools.pileup(clr, loops, view_df=view, expected_df=exp, expected_value_col="balanced.avg",
                                        flank=100_000, nproc=1), float)
    if stack.shape[0] != len(loops):                    # older layout: features on the last axis
        stack = np.moveaxis(stack, -1, 0)
    M = np.nanmean(stack, axis=0)
    sec = time.perf_counter() - t
    return {"implementation": f"cooltools {cooltools.__version__} pileup", "n_loops": int(len(loops)),
            **apa_from_window(M), "window_shape": list(M.shape), "seconds": sec}


# ================================================================================================ H4 MDS
def mds_vs_imaging_chr21(roots: Roots) -> dict:
    """``uchrom.recon.bulk.mds.reconstruct_mds`` (SMACOF from a classical-MDS start; distance-decay prior 0.05,
    alpha 4) on the raw Su 2020 Hi-C, against ``uchrom.fea.distance_map`` medians of the tracing (U-Chrom readers)."""
    from scipy import stats

    from uchrom.fea.distance import distance_map
    from uchrom.io import read_su2020_hic, read_su2020_tracing, write_locus_cool
    from uchrom.recon.bulk.mds import reconstruct_mds

    d = roots("su2020_tsv")
    loci, counts = read_su2020_hic(d / "Hi-C_contacts_chromosome21.tsv")
    cd = read_su2020_tracing(d / "chromosome21.tsv", extras=False)
    dm = distance_map(cd, "chr21", stat="median")
    img = pd.DataFrame(dm["matrix"], index=dm["bins"]["start"].to_numpy(), columns=dm["bins"]["start"].to_numpy())
    out, runs = {}, {"": dict(balance=False, alpha=4.0, weight=0.05),
                     "_alpha3": dict(balance=False, alpha=3.0, weight=0.05),
                     "_zeros_left_out": dict(balance=False, alpha=4.0, weight=0.0),
                     "_ice_balanced": dict(balance=True, alpha=4.0, weight=0.05)}
    with tempfile.TemporaryDirectory() as tmp:
        cool = write_locus_cool(loci, counts, Path(tmp) / "su2020_chr21.cool", assembly="hg38", balance=True)
        for tag, kw in runs.items():
            t = time.perf_counter()
            st = reconstruct_mds(cool, chrom="chr21", method="smacof", n_iter=1000, tol=1e-6, device="cpu",
                                 key_added=None, **kw)
            sec = time.perf_counter() - t
            sp = st.spots_with_loci()
            xyz = np.asarray(st.coords, float)
            starts = sp["start"].to_numpy(np.int64)
            D = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1))
            I = img.reindex(index=starts, columns=starts).to_numpy(float)
            iu = np.triu_indices(len(starts), 1)
            ok = np.isfinite(I[iu])
            out[f"rho{tag}"] = float(stats.spearmanr(D[iu][ok], I[iu][ok]).statistic)
            out[f"n_pairs{tag}"] = int(ok.sum())
            out[f"seconds{tag}"] = sec
    out["n_loci"] = int(len(loci))
    return out


def mds_vs_imaging_chr21_independent(roots: Roots) -> dict:
    """numpy: the prior and the distances as stated; scikit-learn's metric MDS (SMACOF, several random starts);
    imaging medians in numpy."""
    import sklearn
    from sklearn.manifold import MDS

    d = roots("su2020_tsv")
    xyz, labels, C = su_tables(d)
    med, _, nobs = per_pair_imaging(xyz)
    D = decay_prior_distances(C)
    iu = np.triu_indices(len(D), 1)
    out = {"implementation": f"scikit-learn {sklearn.__version__} MDS (metric SMACOF)"}
    for seed in (0, 1, 2):
        kw = dict(n_components=3, n_init=4, max_iter=3000, eps=1e-9, random_state=seed)
        t = time.perf_counter()
        try:                                          # scikit-learn >= 1.8 names
            mds = MDS(metric="precomputed", init="random", **kw)
            X = mds.fit_transform(D)
        except TypeError:
            mds = MDS(dissimilarity="precomputed", metric=True, **kw)
            X = mds.fit_transform(D)
        sec = time.perf_counter() - t
        S = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
        out[f"rho_seed{seed}"] = spearman(S[iu], med[iu])
        out[f"seconds_seed{seed}"] = sec
    out["rho"] = float(np.mean([out[f"rho_seed{s}"] for s in (0, 1, 2)]))
    out["seconds"] = float(np.mean([out[f"seconds_seed{s}"] for s in (0, 1, 2)]))
    out["n_pairs"] = int(np.isfinite(med[iu]).sum())
    out["rho_target_distances_vs_imaging"] = spearman(D[iu], med[iu])   # the input itself (no embedding)
    return out


# ================================================================================================ XM1
def imaging_frequency_vs_hic_chr21(roots: Roots) -> dict:
    """``uchrom.fea.distance.contact_frequency`` (0.5 µm) on the U-Chrom reader's traces against the Hi-C."""
    from uchrom.fea.distance import contact_frequency
    from uchrom.io import read_su2020_hic, read_su2020_tracing

    d = roots("su2020_tsv")
    t = time.perf_counter()
    loci, counts = read_su2020_hic(d / "Hi-C_contacts_chromosome21.tsv")
    cd = read_su2020_tracing(d / "chromosome21.tsv", extras=False)
    freq, bins, n_traces = contact_frequency(cd, 0.5, "chr21")
    starts = np.array([b[0] for b in bins], np.int64)
    F = pd.DataFrame(freq, index=starts, columns=starts).reindex(index=loci["start"].to_numpy(),
                                                                 columns=loci["start"].to_numpy()).to_numpy()
    iu = np.triu_indices(len(loci), 1)
    mids = ((loci["start"] + loci["end"]) / 2).to_numpy(float)
    near = np.abs(mids[iu[0]] - mids[iu[1]]) <= 2_000_000
    out = {"rho_all": spearman(F[iu], counts[iu]), "rho_le_2mb": spearman(F[iu][near], counts[iu][near]),
           "n_pairs_all": int(np.isfinite(F[iu]).sum()), "n_pairs_le_2mb": int(near.sum()), "n_copies": int(n_traces),
           "seconds": time.perf_counter() - t}
    return out


def imaging_frequency_vs_hic_chr21_independent(roots: Roots) -> dict:
    t = time.perf_counter()
    xyz, labels, C = su_tables(roots("su2020_tsv"))
    _, freq, nobs = per_pair_imaging(xyz, 500.0)
    iu = np.triu_indices(len(labels), 1)
    mids = label_mids(labels)
    near = np.abs(mids[iu[0]] - mids[iu[1]]) <= 2_000_000
    return {"implementation": "numpy", "rho_all": spearman(freq[iu], C[iu]),
            "rho_le_2mb": spearman(freq[iu][near], C[iu][near]), "n_pairs_all": int(np.isfinite(freq[iu]).sum()),
            "n_pairs_le_2mb": int(near.sum()), "n_copies": int(xyz.shape[0]), "seconds": time.perf_counter() - t}


# ================================================================================================ P1 MSR
# --- the benchmark's Moran spectral randomisation, copied from benchmarks/spatial_variable/spatial_variable.py of
#     branch Nanguage/spatial-nulls (PR #127, commit 52c3f77): _msr_basis (its graph given here as a matrix) and msr
def pr127_msr_basis(A: np.ndarray):
    n = len(A)
    r = A.sum(axis=1)
    big = 10.0 * max(1.0, float(r.max()))
    M = A - (r[:, None] + r[None, :]) / n + r.sum() / n ** 2      # H A H, H = I - 11'/n
    M -= big / n                                                # the constant: eigenvalue -big, alone
    lam, V = np.linalg.eigh(M)
    del M
    keep = np.ones(n, bool)
    keep[int(np.argmin(lam))] = False
    return V[:, keep], lam[keep]


def pr127_msr(x: np.ndarray, V: np.ndarray, rng, band: int = 0) -> np.ndarray:
    x = np.asarray(x, float)
    x = x - x.mean(axis=0)
    a = V.T @ x
    if band <= 0:
        return V @ (a * rng.choice([-1.0, 1.0], size=a.shape))
    b = np.empty_like(a)
    for lo in range(0, len(a), band):
        blk = a[lo:lo + band]
        z = rng.standard_normal(blk.shape)
        b[lo:lo + band] = z * (np.linalg.norm(blk, axis=0) / np.linalg.norm(z, axis=0))
    return V @ b


def _r2(y: np.ndarray, X: np.ndarray) -> float:
    Xi = np.c_[np.ones(len(y)), X]
    beta = np.linalg.lstsq(Xi, y, rcond=None)[0]
    res = y - Xi @ beta
    return float(1 - (res ** 2).sum() / ((y - y.mean()) ** 2).sum())


def _null_summary(obs: float, null: np.ndarray) -> dict:
    return {"null_median": float(np.median(null)), "null_q95": float(np.quantile(null, 0.95)),
            "p": float((1 + (null >= obs).sum()) / (len(null) + 1))}


def _melanoma_legs(store: Path, key: str, cells) -> tuple:
    """Per spot: legs per chromosome and the short / long counts of the question (1-Mb maps, autosomes; bin offsets 1
    and 19-53, i.e. separations of 1 Mb and 19-53 Mb counted as offset x bin size) — from the exported per-spot
    coolers (chromdata + cooler only)."""
    import cooler
    from chromdata.embedded import export_scool

    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR", "/tmp")) as tmp:
        sc = Path(tmp) / "r1.scool"
        export_scool(str(store), key, sc, cells=list(cells), trans=True)
        names = list(cooler.fileops.list_coolers(str(sc)))
        rows, legs, short, long_ = [], {}, {}, {}
        first = cooler.Cooler(f"{sc}::{names[0]}")
        bins = first.bins()[["chrom", "start", "end"]][:]
        chrom = bins["chrom"].astype(str).to_numpy()
        auto = np.array([c.replace("chr", "").isdigit() for c in chrom])
        for nm in names:
            clr = cooler.Cooler(f"{sc}::{nm}")
            px = clr.pixels()[:]
            b1, b2 = px["bin1_id"].to_numpy(), px["bin2_id"].to_numpy()
            cnt = px["count"].to_numpy(float)
            ca, cb = chrom[b1], chrom[b2]
            cid = nm.rsplit("/", 1)[-1]
            legs[cid] = (pd.Series(cnt, index=ca).groupby(level=0).sum()
                         .add(pd.Series(cnt, index=cb).groupby(level=0).sum(), fill_value=0))
            cis = (ca == cb) & auto[b1]
            off = np.abs(b2.astype(np.int64) - b1.astype(np.int64))      # bin offset (1-Mb bins)
            short[cid] = float(cnt[cis & (off == 1)].sum())
            long_[cid] = float(cnt[cis & (off >= 19) & (off <= 53)].sum())
    L = pd.DataFrame(legs).T.fillna(0.0)
    return L, pd.Series(short), pd.Series(long_)


def _melanoma_legs_fast(store: Path, key: str, cells) -> tuple:
    """The same three tables as :func:`_melanoma_legs`, read straight from the store's embedded per-spot map (one pass
    over each partition, chromdata's reader; no cooler files): the independent path, and the read an agent can do."""
    from chromdata.embedded import _open

    e = _open(str(store), key)
    r = e.resolutions[0]
    pos = {c: i for i, c in enumerate(e.cells)}
    want = [str(c) for c in cells if str(c) in pos]
    order = np.unique([pos[c] for c in want]).astype(np.int64)
    row_of = np.full(len(e.cells), -1, np.int64)
    row_of[order] = np.arange(len(order))
    chroms = list(e.chroms)
    firsts = np.array([e.first_bin(c, r) for c in chroms], np.int64)
    legs = np.zeros((len(order), len(chroms)))
    short = np.zeros(len(order))
    long_ = np.zeros(len(order))
    parts = [(f"cis/{c}", ci) for ci, c in enumerate(chroms) if f"cis/{c}" in e.g[f"r{r}"]]
    if "trans" in e.g[f"r{r}"]:
        parts.append(("trans", -1))
    for part, ci in parts:
        off = e._offsets(r, part, "cell_offset")
        n = off[order + 1] - off[order]
        cell_rows = np.repeat(row_of[order], n)                      # the cell of every pixel, in read order
        k = 0
        auto = ci >= 0 and chroms[ci].lower().replace("chr", "").isdigit()
        for d in e.cell_pixel_blocks(part, order, r):
            u = cell_rows[k:k + len(d)]
            k += len(d)
            b1, b2, cnt = d[:, 0].astype(np.int64), d[:, 1].astype(np.int64), d[:, 2].astype(np.float64)
            if ci >= 0:                                              # cis: two legs on this chromosome
                legs[:, ci] += 2 * np.bincount(u, weights=cnt, minlength=len(order))
                if auto:
                    o = np.abs(b2 - b1)
                    short += np.bincount(u[o == 1], weights=cnt[o == 1], minlength=len(order))
                    m = (o >= 19) & (o <= 53)
                    long_ += np.bincount(u[m], weights=cnt[m], minlength=len(order))
            else:                                                    # trans: one leg to each chromosome
                c1 = np.searchsorted(firsts, b1, side="right") - 1
                c2 = np.searchsorted(firsts, b2, side="right") - 1
                np.add.at(legs, (u, c1), cnt)
                np.add.at(legs, (u, c2), cnt)
        if k != len(cell_rows):
            raise RuntimeError(f"{part}: read {k} pixels, expected {len(cell_rows)}")
    ids = [e.cells[i] for i in order]
    return (pd.DataFrame(legs, index=ids, columns=chroms), pd.Series(short, index=ids), pd.Series(long_, index=ids))


def _section_spots(cd) -> pd.DataFrame:
    cells = cd.cells.copy()
    cells.index = cells.index.astype(str)
    return cells[cells["section"].astype(str) == SECTION]


def _section_key(cd) -> str:
    """The per-spot map of the section (keys of uns['linked_scool'], chromdata's reader)."""
    from chromdata.cdata import _linked_records

    keys = [k for k in _linked_records(cd.uns.get("linked_scool")) if SECTION in k]
    if len(keys) != 1:
        raise SystemExit(f"per-spot map of {SECTION}: {keys}")
    return keys[0]


def _cn_pcs(L: pd.DataFrame, k: int = 5) -> np.ndarray:
    L = L.drop(columns=[c for c in L.columns if str(c) in ("chrY", "Y")], errors="ignore")
    share = L.div(L.sum(axis=1), axis=0)
    prof = np.log2(share / share.median()).replace([np.inf, -np.inf], np.nan)
    Z = ((prof - prof.mean()) / prof.std(ddof=1)).fillna(0.0).to_numpy(float)
    u, s, _ = np.linalg.svd(Z, full_matrices=False)
    return u[:, :k] * s[:k], list(L.columns)


def msr_copy_number_share(roots: Roots) -> dict:
    """y = U-Chrom's separation_profile (log2 short / long); copy-number PCs from the contact legs (as version 2's
    S3); spot graph = U-Chrom's spatial_neighbors (radius 1.5 on array_col / array_row = the 8 adjacent spots);
    the null = PR #127's MSR code (bands of 20), 1,000 draws x MC_SEEDS seeds."""
    from chromdata import ChromData

    from uchrom.fea import separation_profile
    from uchrom.spatial import spatial_neighbors

    store = roots("melanoma")
    t0 = time.perf_counter()
    cd = ChromData.read(str(store), backed=True)
    g = _section_spots(cd)
    key = _section_key(cd)
    sp = separation_profile(cd, source="contacts", contacts=key, cells=list(g.index), key_added=None)
    sp.index = sp.index.astype(str)
    L, _, _ = _melanoma_legs(store, key, g.index)
    L.index = L.index.astype(str)
    ids = [c for c in g.index if c in L.index and np.isfinite(sp.loc[c, "log2_short_long"])]
    y = sp.loc[ids, "log2_short_long"].to_numpy(float)
    C, chroms = _cn_pcs(L.loc[ids])
    obs = _r2(y, C)
    t_feat = time.perf_counter() - t0
    gg = g.loc[ids]
    cdq = ChromData(coords=np.zeros((0, 3)), spots=pd.DataFrame({"chrom": pd.Series([], dtype=str), "start": [],
                                                                 "end": [], "trace_id": []}),
                    cells=pd.DataFrame({"cell_id": ids}))
    cdq.set_cell_positions(gg[["array_col", "array_row"]].to_numpy(float), cell_ids=ids, unit="spot", frame="tissue")
    graph = spatial_neighbors(cdq, method="radius", radius=1.5)
    if list(map(str, graph.cell_ids)) != ids:
        raise SystemExit("graph order")
    A = graph.weights("b").toarray()
    t1 = time.perf_counter()
    V, lam = pr127_msr_basis(A)
    t_basis = time.perf_counter() - t1
    runs = []
    for seed in range(MC_SEEDS):
        rng = np.random.default_rng(seed)
        null = np.array([_r2(y, pr127_msr(C, V, rng, band=20)) for _ in range(N_RANDOM)])
        runs.append(_null_summary(obs, null))
    t_all = time.perf_counter() - t0
    return {"r2_observed": obs, **{k: runs[0][k] for k in runs[0]}, "mc": runs, "n_spots": int(len(ids)),
            "n_links": int(A.sum()), "chroms": chroms, "median_y": float(np.median(y)),
            "seconds_features": t_feat, "seconds_basis": t_basis,
            "seconds_one_run_of_1000": (t_all - t_feat - t_basis) / MC_SEEDS, "seconds": t_all}


def msr_copy_number_share_independent(roots: Roots) -> dict:
    """numpy written from Wagner & Dray 2015: y and the legs from the per-spot pixels (one pass over the embedded map),
    PCs by SVD, the basis
    from an orthonormal complement of the constant (Q' A Q), bands of 20 from the smallest eigenvalue."""
    from chromdata import ChromData

    store = roots("melanoma")
    t0 = time.perf_counter()
    cd = ChromData.read(str(store), backed=True)
    g = _section_spots(cd)
    key = _section_key(cd)
    L, short, long_ = _melanoma_legs_fast(store, key, g.index)
    L.index = L.index.astype(str)
    short.index, long_.index = short.index.astype(str), long_.index.astype(str)
    ids = [c for c in g.index if c in L.index]
    y = np.log2((short.loc[ids].to_numpy() + 0.5) / (long_.loc[ids].to_numpy() + 0.5))
    C, chroms = _cn_pcs(L.loc[ids])
    obs = _r2(y, C)
    t_feat = time.perf_counter() - t0
    r = g.loc[ids, "array_row"].to_numpy(int)
    c = g.loc[ids, "array_col"].to_numpy(int)
    A = ((np.abs(r[:, None] - r[None, :]) <= 1) & (np.abs(c[:, None] - c[None, :]) <= 1)).astype(float)
    np.fill_diagonal(A, 0.0)
    n = len(A)
    t1 = time.perf_counter()
    # orthonormal basis of the complement of the constant: QR of [1 | I], columns 2..n
    Q = np.linalg.qr(np.c_[np.ones(n), np.eye(n)[:, :n - 1]])[0][:, 1:]
    lam, W = np.linalg.eigh(Q.T @ A @ Q)                      # ascending
    V = Q @ W
    del Q, W
    t_basis = time.perf_counter() - t1
    Cc = C - C.mean(axis=0)
    a = V.T @ Cc                                              # coefficients, (n - 1) x 5
    bands = [slice(lo, min(lo + 20, len(a))) for lo in range(0, len(a), 20)]
    norms = np.stack([np.linalg.norm(a[b], axis=0) for b in bands])          # bands x 5
    runs = []
    for seed in range(MC_SEEDS):
        rng = np.random.default_rng(10_000 + seed)
        null = np.empty(N_RANDOM)
        for k in range(N_RANDOM):
            z = rng.standard_normal(a.shape)
            for bi, b in enumerate(bands):
                z[b] *= norms[bi] / np.linalg.norm(z[b], axis=0)
            null[k] = _r2(y, V @ z)
        runs.append(_null_summary(obs, null))
    t_all = time.perf_counter() - t0
    return {"implementation": "numpy", "r2_observed": obs, **{k: runs[0][k] for k in runs[0]}, "mc": runs,
            "n_spots": int(n), "n_links": int(A.sum()), "chroms": chroms, "median_y": float(np.median(y)),
            "seconds_features": t_feat, "seconds_basis": t_basis,
            "seconds_one_run_of_1000": (t_all - t_feat - t_basis) / MC_SEEDS, "seconds": t_all}


# ================================================================================================ P2
def _takei_cells(cd) -> pd.DataFrame:
    c = cd.cells.copy()
    c.index = c.index.astype(str)
    return c


def partial_spearman_rna_rg(roots: Roots) -> dict:
    """Rg per trace = uchrom.fea.distance.radius_of_gyration; RNA totals = uchrom.io.read_seqfish_transcriptome on
    the Zenodo zip (and checked against the question's table); partial Spearman as stated (scipy ranks)."""
    from scipy import stats

    from chromdata import ChromData
    from uchrom.fea.distance import radius_of_gyration
    from uchrom.io.seqfish_multiomics import read_seqfish_transcriptome

    t = time.perf_counter()
    cd = ChromData.read(str(roots("takei")))
    cells = _takei_cells(cd)
    rg = radius_of_gyration(cd)
    n = cd.spots.groupby("trace_id", observed=True).size()
    cell_of = cd.spots.groupby("trace_id", observed=True)["cell_id"].first().astype(str)
    keep = n.index[n >= MIN_SPOTS]
    med = rg.reindex(keep).groupby(cell_of.reindex(keep).to_numpy()).median()
    ad = read_seqfish_transcriptome(roots.dl / "takei2025_rna" / "transcriptomic-data.zip", dataset="cerebellum",
                                    assay="mrna", replicate=1, observed_cell_ids=list(cells.index))
    tot = pd.Series(np.asarray(ad.obs["n_transcripts"], float), index=ad.obs_names.astype(str))
    table = pd.read_csv(roots("takei_rna"), sep="\t", index_col=0)
    table.index = table.index.astype(str)
    same = bool(np.array_equal(table.sum(axis=1).reindex(tot.index).to_numpy(float), tot.to_numpy()))
    ids = [c for c in cells.index if c in med.index]
    x, yv = tot.reindex(ids).fillna(0).to_numpy(float), med.reindex(ids).to_numpy(float)
    rx, ry = stats.rankdata(x), stats.rankdata(yv)
    T = pd.get_dummies(cells.loc[ids, "cell_type"].astype(str)).to_numpy(float)          # indicators (no intercept)
    ex = rx - T @ np.linalg.lstsq(T, rx, rcond=None)[0]
    ey = ry - T @ np.linalg.lstsq(T, ry, rcond=None)[0]
    return {"partial_spearman": float(stats.pearsonr(ex, ey).statistic),
            "spearman": float(stats.spearmanr(x, yv).statistic), "n_cells": int(len(ids)),
            "n_traces_used": int(len(keep)), "n_types": int(T.shape[1]), "rna_table_equals_uchrom_reader": same,
            "median_total_transcripts": float(np.median(x)), "median_cell_rg_um": float(np.median(yv)),
            "seconds": time.perf_counter() - t}


def partial_spearman_rna_rg_independent(roots: Roots) -> dict:
    from chromdata import ChromData

    t = time.perf_counter()
    cd = ChromData.read(str(roots("takei")))
    cells = _takei_cells(cd)
    sp = cd.spots[["trace_id", "cell_id"]].astype(str).reset_index(drop=True)
    xyz = pd.DataFrame(np.asarray(cd.coords, float), columns=["x", "y", "z"])
    df = pd.concat([sp, xyz], axis=1)
    df = df[np.isfinite(xyz.to_numpy()).all(1)]
    grp = df.groupby("trace_id")
    cen = grp[["x", "y", "z"]].transform("mean")
    df = df.assign(r2=((df[["x", "y", "z"]] - cen) ** 2).sum(axis=1))
    per = df.groupby("trace_id").agg(n=("r2", "size"), ms=("r2", "mean"), cell=("cell_id", "first"))
    per = per[per["n"] >= MIN_SPOTS]
    med = np.sqrt(per["ms"]).groupby(per["cell"]).median()
    table = pd.read_csv(roots("takei_rna"), sep="\t", index_col=0)
    table.index = table.index.astype(str)
    tot = table.sum(axis=1)
    ids = [c for c in cells.index if c in med.index]
    res = ranks_partial(tot.reindex(ids).fillna(0).to_numpy(float), med.reindex(ids).to_numpy(float),
                        cells.loc[ids, "cell_type"].astype(str).to_numpy())
    return {"implementation": "numpy / pandas", **res, "n_traces_used": int(len(per)), "seconds": time.perf_counter() - t}


# ================================================================================================ finalize
def finalize(out: dict) -> dict:
    """Graded values from the two implementations (rules in v3/reports/bulk.md)."""
    def both(q):
        v = out.get(q, {}).get("verify", {})
        return v.get("uchrom"), next((x for k, x in v.items() if k != "uchrom"), None)

    u, i = both("H1")
    if u and i:
        lo, hi = sorted((u["n_loops"], i["n_loops"]))
        mid = (lo + hi) / 2
        tol = float(np.ceil(((hi - mid) / mid + 0.01) * 100) / 100)     # covers both, + 1 % of the midpoint
        out["H1"]["values"] = {"n_loops": mid, "n_loops_uchrom": u["n_loops"], "n_loops_mustache_hic": i["n_loops"],
                               "band": [mid * (1 - tol), mid * (1 + tol)]}
        out["H1"]["grade"] = {"numbers": ["n_loops"], "tol": tol}
    u, i = both("H2")
    if u and i:
        rf = float(np.floor((min(u["recall"], i["recall"]) - 0.05) * 100) / 100)
        pf = float(np.floor((min(u["precision"], i["precision"]) - 0.05) * 100) / 100)
        out["H2"]["values"] = {"recall": u["recall"], "precision": u["precision"], "n_loops": u["n_loops"],
                               "recall_floor": rf, "precision_floor": pf}
    u, i = both("H3")
    if u and i:
        out["H3"]["values"] = {"n_loops": u["n_loops"], "apa": u["apa"]}
    u, i = both("H4")
    if u and i:
        fl = float(np.floor((min(u["rho"], i["rho"]) - 0.05) * 100) / 100)
        out["H4"]["values"] = {"rho": u["rho"], "rho_floor": fl}
    u, i = both("XM1")
    if u and i:
        out["XM1"]["values"] = {"rho_all": u["rho_all"], "rho_le_2mb": u["rho_le_2mb"]}
    u, i = both("P1")
    if u and i:
        runs = u["mc"] + i["mc"]
        val = {"r2_observed": u["r2_observed"]}
        tols = {}
        for k in ("null_median", "null_q95", "p"):
            xs = np.array([r[k] for r in runs])
            val[k] = float(np.mean(xs))
            tols[k] = float(np.max(np.abs(xs - xs.mean())) / abs(xs.mean())) if xs.mean() else 0.0
        val["mc_relative_spread"] = tols
        out["P1"]["values"] = val
    u, i = both("P2")
    if u and i:
        out["P2"]["values"] = {"partial_spearman": u["partial_spearman"], "spearman": u["spearman"]}
    return out


# ================================================================================================ main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--impl", choices=("uchrom", "independent", "finalize"), required=True)
    ap.add_argument("--data", type=Path, default=Path("/data"))
    ap.add_argument("--data-v3", type=Path, default=Path("/data_v3"))
    ap.add_argument("--dl", type=Path, default=Path("/dl"))
    ap.add_argument("--only", default="")
    ap.add_argument("--out", type=Path, default=HERE.parent / "results" / "bulk.json")
    a = ap.parse_args()
    out = json.loads(a.out.read_text()) if a.out.exists() else {}
    if a.impl == "finalize":
        a.out.write_text(json.dumps(finalize(out), indent=1, default=float))
        return
    roots = Roots(a.data, a.data_v3, a.dl)
    (roots.v3 / "bulk_work").mkdir(parents=True, exist_ok=True)
    only = [q for q in (a.only.split(",") if a.only else QUESTIONS)]
    for q in only:
        name = QUESTIONS[q] + ("" if a.impl == "uchrom" else "_independent")
        t = time.perf_counter()
        try:
            vals = globals()[name](roots)
        except Exception as e:                                        # keep the other questions going
            import traceback

            traceback.print_exc()
            vals = {"error": f"{type(e).__name__}: {e}"}
        sec = time.perf_counter() - t
        rec = out.setdefault(q, {"values": {}, "verify": {}, "truth": QUESTIONS[q], "runtime_s": {}})
        impl = "uchrom" if a.impl == "uchrom" else vals.get("implementation", "independent").split(" (")[0]
        rec["verify"][impl] = vals
        rec["runtime_s"][impl] = round(sec, 2)
        print(q, impl, json.dumps(vals, default=float)[:1500], flush=True)
        a.out.write_text(json.dumps(out, indent=1, default=float))
    if a.impl == "independent":
        print("uchrom imported:", "uchrom" in sys.modules)
        if "uchrom" in sys.modules:
            raise SystemExit("the independent path imported uchrom")


if __name__ == "__main__":
    main()
