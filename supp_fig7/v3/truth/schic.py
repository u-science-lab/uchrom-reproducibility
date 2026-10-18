"""Reference answers of the schic family of version 3 (questions/schic.yaml), computed with U-Chrom (the version the
agents get: commit fb43ac4, environment venv-uchrom).

    python schic.py --data /data --data3 /data_v3 --out schic_uchrom.json [--only R1,R3]

One function per question (``truth:`` in schic.yaml).  Each returns the graded values plus diagnostics; the runtime
of the U-Chrom computation is under ``runtime_s``.  The independent implementations are in ``schic_ref_*.py``
(authors' code in their own environments, numpy / cooltools); ``schic_collect.py`` merges both into
``results/schic.json``.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

DATA = {"stevens": "stevens2017_cell1", "schicar": "schicar_mop/schicar_mop.chromdata.zarr",
        "scool": "schicar_mop/contacts_cells_1Mb.scool", "mcool": "schicar_mop/contacts_all.mcool"}
DATA3 = {"stevens2": "stevens2017_cell2", "cpg": "mm10_cpg/mm10_cpg_1Mb.tsv", "nagano2i": "nagano2017_hap_2i"}

#: R3: the two scHiCAR cells (chosen by schic_explore.py: the cell with the median chr1 contacts of Oligo and of L6CT)
R3_CELLS = ("GTGAAATCTACGAGATAC", "TAGACCGGTATGTACCCG")
#: R4 / R5: autosomes
AUTOSOMES = [f"chr{i}" for i in range(1, 20)]
#: R6: Nagano et al. 2017 (schic2 analyzeScHiC.r / hap_serum_2i_params.r)
NEAR = (2 ** 14.5, 2 ** 21)
MITOTIC = (2 ** 21, 2 ** 23.5)
NAGANO_CHROMS = [f"chr{i}" for i in range(1, 20)] + ["chrX"]


# ------------------------------------------------------------------ R1 / R2: NucDynamics vs the published structure
def _pdb_model(path, model: int = 1, chain: str = "a") -> pd.DataFrame:
    """Particles of one model and chain of a Nuc3D PDB export: x, y, z and the sequence position (last column)."""
    import gzip

    rows, m = [], 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("MODEL"):
                m = int(line.split()[1])
            elif m == model and line.startswith(("HETATM", "ATOM")) and line[21] == chain:
                rows.append((float(line[30:38]), float(line[38:46]), float(line[46:54]), int(line.split()[-1])))
            elif m > model:
                break
    return pd.DataFrame(rows, columns=["x", "y", "z", "pos"])


def _binned_distances(xyz: np.ndarray, pos: np.ndarray, res: int = 1_000_000) -> pd.Series:
    """Pairwise distances of the particles averaged into bins of ``res``; Series indexed by (bin i, bin j), i < j."""
    b = pos // res
    df = pd.DataFrame(xyz, columns=["x", "y", "z"]).assign(b=b).groupby("b").mean()
    d = np.sqrt(((df.to_numpy()[:, None, :] - df.to_numpy()[None, :, :]) ** 2).sum(-1))
    i, j = np.triu_indices(len(df), 1)
    return pd.Series(d[i, j], index=pd.MultiIndex.from_arrays([df.index[i], df.index[j]]))


def _nucdyn_vs_published(folder: Path, seeds=(0, 1)):
    from uchrom.recon.sc.nucdyn import reconstruct_nucdyn

    pairs = next(folder.glob("*contact_pairs.txt.gz"))
    pdb = next(folder.glob("*.pdb.gz"))
    pub = {m: _pdb_model(pdb, m) for m in (1, 2)}
    ref = _binned_distances(pub[1][["x", "y", "z"]].to_numpy(), pub[1]["pos"].to_numpy())
    a, b = _binned_distances(pub[2][["x", "y", "z"]].to_numpy(), pub[2]["pos"].to_numpy()).align(ref, join="inner")
    rho_pub = float(st.spearmanr(a, b).statistic)
    rhos, secs = [], []
    for seed in seeds:
        t0 = time.perf_counter()
        cd = reconstruct_nucdyn(str(pairs), n_models=1, seed=seed, particle_sizes=(8e6, 4e6, 2e6, 1e6), device="cpu",
                                key_added=None)
        secs.append(time.perf_counter() - t0)
        sp = cd.spots_with_loci()
        sel = (sp["chrom"].astype(str) == "chr1").to_numpy()
        dist = _binned_distances(np.asarray(cd.coords)[sel], sp.loc[sel, "start"].to_numpy(np.int64))
        x, y = dist.align(ref, join="inner")
        rhos.append(float(st.spearmanr(x, y).statistic))
    floor = round(min(min(rhos), rho_pub) - 0.05, 3)
    return {"rho": float(np.mean(rhos)), "rho_uchrom_by_seed": rhos, "rho_published_model1_vs_model2": rho_pub,
            "rho_floor": floor, "n_bins_published": int(pub[1]["pos"].floordiv(1_000_000).nunique()),
            "n_pairs_compared": int(len(x)),
            "runtime_s": {"uchrom_per_run": secs}}


def nucdyn_cell1(a):
    """R1 (= v2 M3): Stevens 2017 Cell 1."""
    return _nucdyn_vs_published(a.data / DATA["stevens"])


def nucdyn_cell2(a):
    """R2: Stevens 2017 Cell 2, the same protocol and floor rule."""
    return _nucdyn_vs_published(a.data3 / DATA3["stevens2"])


# ------------------------------------------------------------------ R3: scHiCluster RWR imputation of two cells
def _cis_upper(scool: str, cell: str, chrom: str):
    """The cell's pixels of ``chrom`` as stored (upper triangle, diagonal included), dense ``(n, n)``."""
    import cooler

    clr = cooler.Cooler(f"{scool}::/cells/{cell}")
    px = clr.matrix(balance=False, as_pixels=True, join=False).fetch(chrom)
    lo = clr.offset(chrom)
    n = clr.extent(chrom)[1] - lo
    u = np.zeros((n, n))
    np.add.at(u, (px["bin1_id"].to_numpy() - lo, px["bin2_id"].to_numpy() - lo), px["count"].to_numpy(float))
    return u


def schicluster_rwr(u: np.ndarray, *, pad: int = 1, rp: float = 0.5, tol: float = 1e-6, max_iter: int = 30):
    """scHiCluster's imputation of one chromosome as its authors' code does it (``schicluster.schicluster
    .impute_cpu``, Zhou et al. 2019), with U-Chrom's convolution and RWR kernels: A = log2(U + U^T + 1) of the stored
    upper triangle U; (2 pad + 1)^2 mean filter with zeros outside; diagonal set to 0, then 1 on the diagonal of
    the bins without any value; random walk with restart from Q = I; returns the symmetrised (Q + Q^T) / 2."""
    from uchrom.fea.impute_rwr import convolve_maps, rwr_dense

    a = np.log2(u + u.T + 1.0)
    b = convolve_maps(a[None], pad)[0]
    np.fill_diagonal(b, 0.0)
    b[np.diag_indices(len(b))] = (b.sum(axis=0) == 0).astype(float)
    return rwr_dense(b[None], rp=rp, tol=tol, max_iter=max_iter, self_loops=False)[0]


def rwr_pair_corr(a):
    """R3: Pearson r of the two cells' imputed chr1 maps over the pixels i < j."""
    from uchrom.emb.hic import schicluster_impute

    scool = str(a.data / DATA["scool"])
    t0 = time.perf_counter()
    us = [_cis_upper(scool, c, "chr1") for c in R3_CELLS]
    qs = [schicluster_rwr(u) for u in us]
    secs = time.perf_counter() - t0
    n = len(us[0])
    iu = np.triu_indices(n, 1)
    r = float(np.corrcoef(qs[0][iu], qs[1][iu])[0, 1])
    out = {"r": r, "cells": list(R3_CELLS), "n_bins": int(n), "n_pixels": int(len(iu[0])),
           "contacts_chr1": [float(u.sum()) for u in us], "runtime_s": {"uchrom": secs}}
    # how much the definition matters (reported, not graded)
    sym = [u + u.T - np.diag(np.diag(u)) for u in us]
    var = {}
    q = [schicluster_impute(m, pad=1, rp=0.5) for m in sym]
    var["uchrom_schicluster_impute_default_counts"] = float(np.corrcoef(q[0][iu], q[1][iu])[0, 1])
    q = [schicluster_impute(np.log2(u + u.T + 1.0), pad=1, rp=0.5) for u in us]
    var["uchrom_schicluster_impute_default_log2"] = float(np.corrcoef(q[0][iu], q[1][iu])[0, 1])
    iu0 = np.triu_indices(n, 0)
    var["with_diagonal"] = float(np.corrcoef(qs[0][iu0], qs[1][iu0])[0, 1])
    var["spearman"] = float(st.spearmanr(qs[0][iu], qs[1][iu]).statistic)
    var["raw_counts_pearson"] = float(np.corrcoef(sym[0][iu], sym[1][iu])[0, 1])
    out["variants"] = var
    return out


# ------------------------------------------------------------------ R4: scHiCluster embedding, kNN purity
def schicluster_knn(a, k: int = 30):
    """R4: embed_contacts (scHiCluster: pad 1, rp 0.5, top 20 %, 20 PCs per chromosome, chr1-19, PCA to 20)
    and the k-nearest-neighbour purity of cells['cell_type'] in the 20 PCs (k = 30: the neighbourhood size at which
    imputation separates most clearly from no imputation; reports/schic.md)."""
    from chromdata import ChromData
    from uchrom.emb.hic import embed_contacts
    from uchrom.emb.metrics import knn_purity

    cd = ChromData.read(str(a.data / DATA["schicar"]))
    t0 = time.perf_counter()
    embed_contacts(cd, chroms=AUTOSOMES, pad=1, rp=0.5, top_pct=20.0, n_per_chrom=20, n_pcs=20, methods=("pca",),
                   random_state=0, key_added="hic3")
    secs = time.perf_counter() - t0
    x = np.asarray(cd.cellm["hic3_pca"], float)
    labels = cd.cells["cell_type"].astype(str).to_numpy()
    ok = np.isfinite(x).all(axis=1)
    out = {"purity": knn_purity(x[ok], labels[ok], k=k), "k": k, "n_cells": int(ok.sum()),
           "n_types": int(len(set(labels[ok]))), "n_dims": int(x.shape[1]), "runtime_s": {"uchrom": secs}}
    out["purity_k_variants"] = {kk: knn_purity(x[ok], labels[ok], k=kk) for kk in (5, 10, 15)}
    pooled = np.where(np.char.startswith(labels.astype(str), "L23IT"), "L23IT", labels)
    out["purity_L23IT_pooled"] = knn_purity(x[ok], pooled[ok], k=k)
    np.save(a.out.parent / "schic_R4_uchrom_pca.npy", x)
    return out


# ------------------------------------------------------------------ R5: scA/B (scHiCluster's compartment score)
def _cpg(a) -> pd.DataFrame:
    t = pd.read_csv(a.data3 / DATA3["cpg"], sep="\t")
    t["chrom"] = t["chrom"].astype(str)
    return t


def _e1(mcool: str, chroms, cpg: pd.DataFrame, resolution=1_000_000) -> pd.DataFrame:
    from uchrom.strc.comp.eig import EigCompartmentParams, call_compartments_eig

    t = call_compartments_eig(mcool, chrom=chroms, resolution=resolution,
                              phasing=cpg.rename(columns={"cpg_ratio": "value"})[["chrom", "start", "end", "value"]],
                              params=EigCompartmentParams(ignore_diags=2, phasing_method="pearson"), key_added=None)
    col = "E1" if "E1" in t else "e1"
    return t[["chrom", "start", "end", col]].rename(columns={col: "E1"})


def scab_vs_e1(a):
    """R5: per cell scHiCluster's compartment score (cis, diagonal removed, rows normalised, x CpG ratio; bins
    without a cis contact take their own CpG ratio), averaged over the cells, against the cis E1 of the 1-Mb
    pseudo-bulk of all cells (cooltools eigs_cis formulation, oriented by the CpG ratio), chr1-19, bins with CpG
    ratio > 0 and an E1.  Graded: the partial correlation given the CpG ratio (``r_partial_cpg``); the plain r is
    nearly that of the CpG ratio alone and is reported only."""
    from chromdata import ChromData
    from uchrom.strc.comp.single_cell import ScABParams, single_cell_compartments

    cd = ChromData.read(str(a.data / DATA["schicar"]))
    cpg = _cpg(a)
    t0 = time.perf_counter()
    frame = single_cell_compartments(cd, track=cpg.rename(columns={"cpg_ratio": "value"})[
        ["chrom", "start", "end", "value"]], chrom=AUTOSOMES, params=ScABParams(contacts="cis", normalize="none"),
        key_added=None)
    secs_scab = time.perf_counter() - t0
    labels = pd.Index(frame.columns)
    key = cpg["chrom"] + ":" + cpg["start"].astype(str) + "-" + cpg["end"].astype(str)
    own = pd.Series(cpg["cpg_ratio"].to_numpy(), index=key).reindex(labels).to_numpy()
    vals = frame.to_numpy(np.float64)
    n_nan = int(np.isnan(vals).sum())
    vals = np.where(np.isnan(vals), own[None, :], vals)                  # scHiCluster: own CpG without contacts
    t1 = time.perf_counter()
    e1 = _e1(str(a.data / DATA["mcool"]), AUTOSOMES, cpg)
    secs_e1 = time.perf_counter() - t1
    ekey = e1["chrom"].astype(str) + ":" + e1["start"].astype(str) + "-" + e1["end"].astype(str)
    E = pd.Series(e1["E1"].to_numpy(float), index=ekey).reindex(labels).to_numpy()
    cells = frame.index.astype(str)
    ctype = cd.cells["cell_type"].astype(str).reindex(cells).to_numpy()
    ok = np.isfinite(E) & (own > 0)

    def r(v):
        m = ok & np.isfinite(v)
        return float(np.corrcoef(v[m], E[m])[0, 1])

    mean_all = vals.mean(axis=0)

    def pr(v):                                         # partial correlation of v and E1 given the CpG ratio
        m = ok & np.isfinite(v)
        r_se, r_sc, r_ec = (float(np.corrcoef(x, y)[0, 1]) for x, y in ((v[m], E[m]), (v[m], own[m]), (E[m], own[m])))
        return (r_se - r_sc * r_ec) / np.sqrt((1 - r_sc ** 2) * (1 - r_ec ** 2)), r_sc, r_ec

    partial, r_sc, r_ec = pr(mean_all)
    out = {"r_partial_cpg": float(partial), "r": r(mean_all), "r_score_cpg": r_sc, "r_cpg_e1": r_ec,
           "n_bins": int(ok.sum()), "n_cells": int(len(cells)),
           "bins_without_contacts_filled": n_nan, "runtime_s": {"uchrom_scab": secs_scab, "uchrom_e1": secs_e1}}
    skipped = np.nanmean(frame.to_numpy(np.float64), axis=0)
    diag = {"r_cpg_e1": r(own), "r_mean_nan_skipped": r(skipped), "r_partial_nan_skipped": float(pr(skipped)[0]),
            "r_partial_median": float(pr(np.median(vals, axis=0))[0]),
            "spearman": float(st.spearmanr(mean_all[ok], E[ok]).statistic)}
    with np.errstate(invalid="ignore"):
        resid = mean_all - np.polyval(np.polyfit(own[ok], mean_all[ok], 1), own)
    diag["r_residual_given_cpg"] = r(resid)
    by_type = {}
    for ty in pd.Series(ctype).value_counts().index[:6]:
        by_type[ty] = {"r": r(vals[ctype == ty].mean(axis=0)), "r_partial": float(pr(vals[ctype == ty].mean(axis=0))[0])}
    from uchrom.strc.comp.single_cell import normalize_scab
    diag["r_partial_rank_normalised_cells"] = float(pr(np.nanmean(normalize_scab(frame.to_numpy(np.float64), "rank"),
                                                                     axis=0))[0])
    diag["r_by_type_vs_all_cells_E1"] = by_type
    out["diagnostics"] = diag
    pd.DataFrame({"bin": labels, "mean_scab": mean_all, "E1": E, "cpg": own}).to_csv(
        a.out.parent / "schic_R5_uchrom_bins.tsv", sep="\t", index=False)
    return out


# ------------------------------------------------------------------ R6: Nagano 2017 cell-cycle fractions
def _fends(path: Path):
    f = pd.read_csv(path, sep="\t", usecols=["fend", "chr", "coord"], dtype={"chr": str})
    f = f.sort_values("fend")
    chrom = "chr" + f["chr"].astype(str).str.replace("^chr", "", regex=True)
    return f["fend"].to_numpy(np.int64), chrom.to_numpy(), f["coord"].to_numpy(np.int64)


def _phase(near: np.ndarray, mit: np.ndarray) -> np.ndarray:
    """The 'newer' gates of schic2 (ordering.r group_cells_by_metrics; hap_serum_2i_params.r)."""
    g = np.full(len(near), "mid-S/G2", dtype=object)
    g[(near > 0.611) & (near <= 0.77)] = "early-S"
    g[near <= 0.611] = "G1"
    g[(mit >= 0.3) & (near <= 0.42)] = "post-M"
    g[(near > 0.42) & (-1.8 * mit - near + 1 < 0)] = "pre-M"
    return g


def nagano_cellcycle(a):
    """R6: per cell the share of cis contacts (distance > 2^14.5 bp) at 2^14.5 < d <= 2^21 ('near') and
    2^21 < d <= 2^23.5 ('mitotic band'); medians over the 395 cells and the cells per phase by the paper's gates.
    U-Chrom: ContactStats (log2 bins of 1/8 octave from 2^10) + contact_features."""
    from chromdata import ChromData
    from uchrom.fea.contact_stats import ContactFeatureParams, contact_features
    from uchrom.io.pairs_stats import ContactStats, distance_edges

    d = a.data3 / DATA3["nagano2i"]
    cells = pd.read_csv(d / "cells.tsv", sep="\t")["cell"].astype(str).tolist()
    t0 = time.perf_counter()
    fid, fchr, fpos = _fends(d / "GATC.fends")
    stats = ContactStats(cells, NAGANO_CHROMS, distance_edges(base=2.0 ** 10, per_octave=8, max_distance=2.0 ** 27.75),
                         orientation=False)
    cidx = {c: i for i, c in enumerate(NAGANO_CHROMS)}
    n_pairs = []
    for i, cell in enumerate(cells):
        adj = pd.read_csv(d / "adj" / f"{cell}.adj", sep="\t")
        f1, f2 = adj["fend1"].to_numpy(np.int64), adj["fend2"].to_numpy(np.int64)
        lo, hi = np.minimum(f1, f2), np.maximum(f1, f2)
        key = np.unique(lo * (1 << 32) + hi)                       # each fend pair once
        lo, hi = key >> 32, key & ((1 << 32) - 1)
        k1 = np.clip(np.searchsorted(fid, lo), 0, len(fid) - 1)
        k2 = np.clip(np.searchsorted(fid, hi), 0, len(fid) - 1)
        known = (fid[k1] == lo) & (fid[k2] == hi)
        c1 = pd.Series(fchr[k1]).map(cidx).fillna(-1).to_numpy(np.int64)
        c2 = pd.Series(fchr[k2]).map(cidx).fillna(-1).to_numpy(np.int64)
        stats.add(np.where(known, i, -1), c1, fpos[k1], c2, fpos[k2])
        n_pairs.append(len(key))
    no_spots = pd.DataFrame({"chrom": pd.Series(dtype=str), "start": pd.Series(dtype=np.int64),
                             "end": pd.Series(dtype=np.int64), "trace_id": pd.Series(dtype=str),
                             "cell_id": pd.Series(dtype=str)})
    cd = ChromData(np.zeros((0, 3)), no_spots, cells=pd.DataFrame(index=pd.Index(cells, name="cell_id")))
    stats.write(cd, key="contacts")
    p = ContactFeatureParams(min_distance=NEAR[0], bands=(MITOTIC, NEAR), slopes=())
    t = contact_features(cd, key="contacts", params=p, key_added=None)
    secs = time.perf_counter() - t0
    fcols = [c for c in t.columns if c.startswith("frac_") and c not in ("frac_short", "frac_long")]
    mit, near = t[fcols[0]].to_numpy(float), t[fcols[1]].to_numpy(float)
    phase = _phase(near, mit)
    counts = pd.Series(phase).value_counts()
    out = {"median_f_mitotic": float(np.median(mit)), "median_f_near": float(np.median(near)),
           "n_cells": len(cells), "n_cells_f_mitotic_ge_0.3": int((mit >= 0.3).sum()),
           **{f"n_{ph}": int(counts.get(ph, 0)) for ph in ("post-M", "G1", "early-S", "mid-S/G2", "pre-M")},
           "columns": fcols, "median_fend_pairs": float(np.median(n_pairs)), "runtime_s": {"uchrom": secs}}
    pd.DataFrame({"cell": cells, "f_mitotic": mit, "f_near": near, "phase": phase}).to_csv(
        a.out.parent / "schic_R6_uchrom_cells.tsv", sep="\t", index=False)
    return out


# ------------------------------------------------------------------ R7: SnapHiC loops, 50 cells, chr19
def snaphic_chr19(a):
    """R7: SnapHiC with its defaults (whole-chromosome RWR, restart 0.05, candidates 100 kb - 2 Mb, the mm10 filter
    regions) on chr19 of the 50 Nagano 2017 diploid cells of data_v3/nagano2017_dip_chr19 (their 10-kb maps)."""
    from chromdata import ChromData
    from uchrom.strc.loop.snaphic import SnapHiCParams, call_loops_snaphic

    d = a.data3 / "nagano2017_dip_chr19"
    cells = (d / "cells.txt").read_text().split()
    cd = ChromData.read(str(a.shared / "nagano2017_dip_serum/nagano2017_dip_serum.chromdata.zarr"))
    t0 = time.perf_counter()
    loops = call_loops_snaphic(cd, cells=cells, chrom="chr19", params=SnapHiCParams(),
                               exclude=str(a.shared / "snaphic/mm10_filter_regions.txt"), key_added=None)
    secs = time.perf_counter() - t0
    loops.to_csv(a.out.parent / "schic_R7_uchrom_loops.tsv", sep="\t", index=False)
    return {"n_loops": int(len(loops)), "n_cells": len(cells), "columns": list(loops.columns),
            "runtime_s": {"uchrom": secs}}


FUNCS = {f.__name__: f for f in (nucdyn_cell1, nucdyn_cell2, rwr_pair_corr, schicluster_knn, scab_vs_e1,
                                 nagano_cellcycle, snaphic_chr19)}
IDS = {"R1": "nucdyn_cell1", "R2": "nucdyn_cell2", "R3": "rwr_pair_corr", "R4": "schicluster_knn",
       "R5": "scab_vs_e1", "R6": "nagano_cellcycle", "R7": "snaphic_chr19"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--data3", type=Path, required=True)
    ap.add_argument("--shared", type=Path, default=Path("/shared"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", default=",".join(IDS))
    ap.add_argument("--cells", default=None, help="R3 cells, comma-separated (overrides R3_CELLS)")
    ap.add_argument("--export-cells", type=Path, default=None,
                    help="write the scHiCAR cells (cell, cell_type) for the independent scripts and stop")
    a = ap.parse_args()
    if a.export_cells:
        from chromdata import ChromData
        cd = ChromData.read(str(a.data / DATA["schicar"]))
        pd.DataFrame({"cell": cd.cells.index.astype(str), "cell_type": cd.cells["cell_type"].astype(str).to_numpy()}
                     ).to_csv(a.export_cells, sep="\t", index=False)
        return
    if a.cells:
        global R3_CELLS
        R3_CELLS = tuple(a.cells.split(","))
    out = json.loads(a.out.read_text()) if a.out.exists() else {}
    for q in a.only.split(","):
        fn = IDS[q]
        t0 = time.perf_counter()
        try:
            vals = FUNCS[fn](a)
        except Exception as e:  # keep the others going; the error is the result
            import traceback
            traceback.print_exc()
            vals = {"error": repr(e)}
        out[q] = {"values": vals, "truth": fn, "wall_s": time.perf_counter() - t0}
        print(q, json.dumps(vals, default=str)[:600], flush=True)
        a.out.write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
