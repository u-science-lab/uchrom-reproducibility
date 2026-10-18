"""The independent references of R3-R5 with the authors' scHiCluster code (no U-Chrom; environment venv-schicluster:
schicluster 1.3.2 from PyPI and the GitHub tag v1.3.4 checked out in src/).

    python schic_ref_schicluster.py --data /data --data3 /data_v3 --cells-tsv schic_cells.tsv --out schic_ref_schicluster.json
        [--only R3,R4,R5] [--r3-cells A,B] [--workers 8]

* R3: ``schicluster.schicluster.impute_cpu`` (Zhou et al. 2019; log2(A + A^T + 1), 3x3 mean filter, diagonal removed,
  RWR from Q = I until the Frobenius change < 1e-6, at most 30 steps) on per-cell text files ``<cell>_chr1.txt``
  (``bin1 bin2 count`` of the stored upper triangle, as the function reads them); symmetrised (Q + Q^T) / 2.
* R4: the same imputation for every cell and chromosome (chr1-19; the function's body, fed with the pixels
  instead of a text file), binarised as ``hicluster_cpu`` does it (each cell's values above their 80th
  percentile over the n x n matrix), PCA with 20 components per chromosome (``hicluster_cpu`` uses
  int(0.2 * min(cells, n^2)) - 1 = 1,061 here; the question states 20), concatenated, PCA to 20; kNN purity.
* R5: ``single_chrom_compartment`` of scHiCluster v1.3.4 (``schicluster/compartment/call_compartment.py``, executed
  from the checkout: PyPI's 1.3.2 lacks the row normalisation) per cell and chromosome on
  ``cooler.matrix(balance=False, sparse=True).fetch(chrom)``, the CpG table of the question; mean over the cells
  (written for schic_ref_numpy.py, which correlates it with the cooltools E1).
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path("/opt/venv/src/scHiCluster")
AUTOSOMES = [f"chr{i}" for i in range(1, 20)]
_G: dict = {}


def _cell_pixels(scool: str, cell: str, chrom: str) -> np.ndarray:
    """``bin1 bin2 count`` (chromosome bins) of the stored upper triangle."""
    import cooler

    clr = cooler.Cooler(f"{scool}::/cells/{cell}")
    px = clr.matrix(balance=False, as_pixels=True, join=False).fetch(chrom)
    lo = clr.offset(chrom)
    return np.column_stack([px["bin1_id"].to_numpy() - lo, px["bin2_id"].to_numpy() - lo,
                            px["count"].to_numpy(float)]).astype(float)


def r3(a) -> dict:
    import cooler
    import schicluster
    import schicluster.schicluster as S

    scool = str(a.data / "schicar_mop/contacts_cells_1Mb.scool")
    cells = a.r3_cells.split(",")
    first = cooler.Cooler(f"{scool}::/cells/{cells[0]}")
    res = int(first.binsize)
    ngene = int(first.chromsizes["chr1"] / res) + 1                      # as hicluster_cpu sizes the matrices
    tmp = Path(tempfile.mkdtemp())
    t0 = time.perf_counter()
    qs = []
    for c in cells:
        np.savetxt(tmp / f"{c}_chr1.txt", _cell_pixels(scool, c, "chr1"), fmt="%d\t%d\t%g")
        name, q = S.impute_cpu([str(tmp / c), "1", ngene, 1, 0.5])
        qs.append(q.reshape(ngene, ngene))
    secs = time.perf_counter() - t0
    iu = np.triu_indices(ngene, 1)
    sym = [(q + q.T) / 2 for q in qs]
    return {"r": float(np.corrcoef(sym[0][iu], sym[1][iu])[0, 1]),
            "r_unsymmetrised_upper": float(np.corrcoef(qs[0][iu], qs[1][iu])[0, 1]),
            "asymmetry_max_abs": float(max(np.abs(q - q.T).max() for q in qs)),
            "n_bins": ngene, "cells": cells, "schicluster": schicluster.__version__, "runtime_s": secs}


def _impute(args):
    """The body of ``schicluster.schicluster.impute_cpu`` with the pixels in place of ``np.loadtxt``."""
    from scipy.sparse import csr_matrix

    import schicluster.schicluster as S

    D, ngene, pad, rp = args[:4]
    log = args[4] if len(args) > 4 else True
    if len(D) == 0:
        D = np.zeros((1, 3))
    A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape=(ngene, ngene)).toarray()
    A = np.log2(A + A.T + 1) if log else A + A.T - np.diag(np.diag(A))   # variant: counts, diagonal once
    A = S.neighbor_ave_cpu(A, pad)
    Q = A if rp == -1 else S.random_walk_cpu(A, rp)            # rp = -1: no random walk (as impute_cpu)
    return Q.reshape(ngene * ngene)


def _init(scool):
    _G["scool"] = scool


def _pixels_of(args):
    cell, chrom = args
    return _cell_pixels(_G["scool"], cell, chrom)


def r4b(a) -> dict:
    """R4 baseline (not a reference): the same pipeline without imputation (pad 0, no random walk: log2 counts
    binarised at their 80th percentile), to see what the floor rejects."""
    return r4(a, pad=0, rp=-1, tag="baseline_")


def r4v(a) -> dict:
    """R4 variants of a correct scHiCluster (not references): the paper's description without the code's log2, and
    the binarisation over the upper triangle (i < j) instead of the whole n x n matrix, to size the floor's margin."""
    return {"no_log2": r4(a, log=False, tag="nolog_"), "upper_triangle": r4(a, tri=True, tag="tri_")}


def r4(a, k: int = 30, n_pcs: int = 20, pad: int = 1, rp: float = 0.5, tag: str = "", log: bool = True,
       tri: bool = False) -> dict:
    import cooler
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    scool = str(a.data / "schicar_mop/contacts_cells_1Mb.scool")
    cells = pd.read_csv(a.cells_tsv, sep="\t", dtype=str)
    names = cells["cell"].tolist()
    first = cooler.Cooler(f"{scool}::/cells/{names[0]}")
    res = int(first.binsize)
    t0 = time.perf_counter()
    blocks = []
    timing = {}
    with ProcessPoolExecutor(a.workers, initializer=_init, initargs=(scool,)) as ex:
        for chrom in AUTOSOMES:
            t1 = time.perf_counter()
            ngene = int(first.chromsizes[chrom] / res) + 1
            pix = list(ex.map(_pixels_of, [(c, chrom) for c in names], chunksize=64))
            q = np.array(list(ex.map(_impute, [(p, ngene, pad, rp, log) for p in pix], chunksize=32)))
            if tri:                                                       # variant: entries i < j only
                iu = np.triu_indices(ngene, 1)
                q = (q.reshape(len(q), ngene, ngene) + q.reshape(len(q), ngene, ngene).transpose(0, 2, 1))[:, iu[0], iu[1]] / 2
            thres = np.percentile(q, 80, axis=1)                          # hicluster_cpu, prct = 20
            q = (q > thres[:, None])
            blocks.append(PCA(n_components=n_pcs, random_state=0).fit_transform(q.astype(np.float64)))
            timing[chrom] = time.perf_counter() - t1
            print(chrom, f"{timing[chrom]:.1f} s", flush=True)
            del q, pix
    feats = np.concatenate(blocks, axis=1)
    x = PCA(n_components=n_pcs, random_state=0).fit_transform(feats)
    secs = time.perf_counter() - t0
    labels = cells["cell_type"].to_numpy()
    out = {"k": k, "n_cells": len(names), "runtime_s": secs, "runtime_s_by_chrom": timing}
    for kk in (k, 5, 10, 15):
        nn = NearestNeighbors(n_neighbors=kk + 1).fit(x)
        idx = nn.kneighbors(x, return_distance=False)[:, 1:]
        out["purity" if kk == k else f"purity_k{kk}"] = float((labels[idx] == labels[:, None]).mean())
    pooled = np.where(pd.Series(labels).str.startswith("L23IT"), "L23IT", labels)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(x)
    idx = nn.kneighbors(x, return_distance=False)[:, 1:]
    out["purity_L23IT_pooled"] = float((pooled[idx] == pooled[:, None]).mean())
    np.save(a.out.parent / f"schic_R4_schicluster_{tag}pca.npy", x)
    return out


def _compartment_function():
    """``single_chrom_compartment`` (+ ``compartment_strength``) of the v1.3.4 checkout, executed as written."""
    import scipy.sparse as sp

    src = (SRC / "schicluster/compartment/call_compartment.py").read_text()
    keep = [n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef)
            and n.name in ("single_chrom_compartment", "compartment_strength")]
    ns = {"np": np, "pd": pd, "diags": sp.diags, "csr_matrix": sp.csr_matrix}
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(SRC / "call_compartment.py"), "exec"), ns)
    return ns["single_chrom_compartment"]


def _comp_cell(cell):
    import cooler

    func = _G.get("func") or _G.setdefault("func", _compartment_function())
    clr = cooler.Cooler(f"{_G['scool']}::/cells/{cell}")
    out = []
    for chrom in AUTOSOMES:
        m = clr.matrix(balance=False, sparse=True).fetch(chrom)
        comp, _ = func(m, _G["cpg"][chrom])
        out.append(np.asarray(comp, float))
    return np.concatenate(out)


def _init5(scool, cpg):
    _G["scool"] = scool
    _G["cpg"] = cpg


def r5(a) -> dict:
    scool = str(a.data / "schicar_mop/contacts_cells_1Mb.scool")
    cells = pd.read_csv(a.cells_tsv, sep="\t", dtype=str)["cell"].tolist()
    cpg = pd.read_csv(a.data3 / "mm10_cpg/mm10_cpg_1Mb.tsv", sep="\t")
    by = {c: cpg.loc[cpg["chrom"] == c, "cpg_ratio"].reset_index(drop=True) for c in AUTOSOMES}
    t0 = time.perf_counter()
    with ProcessPoolExecutor(a.workers, initializer=_init5, initargs=(scool, by)) as ex:
        vals = np.array(list(ex.map(_comp_cell, cells, chunksize=32)))
    secs = time.perf_counter() - t0
    sel = cpg[cpg["chrom"].isin(AUTOSOMES)]
    sel = pd.concat([cpg[cpg["chrom"] == c] for c in AUTOSOMES])
    out = pd.DataFrame({"chrom": sel["chrom"].to_numpy(), "start": sel["start"].to_numpy(),
                        "end": sel["end"].to_numpy(), "mean_comp": vals.mean(axis=0)})
    out.to_csv(a.out.parent / "schic_R5_schicluster_mean.tsv", sep="\t", index=False)
    commit = (SRC / ".git/HEAD").read_text().strip() if (SRC / ".git/HEAD").exists() else None
    return {"n_cells": len(cells), "n_bins": int(vals.shape[1]), "runtime_s": secs, "checkout_head": commit}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--data3", type=Path, required=True)
    ap.add_argument("--cells-tsv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", default="R3,R4,R5")
    ap.add_argument("--r3-cells", default="GTGAAATCTACGAGATAC,TAGACCGGTATGTACCCG")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    out = json.loads(a.out.read_text()) if a.out.exists() else {}
    for q in a.only.split(","):
        t0 = time.perf_counter()
        try:
            v = {"R3": r3, "R4": r4, "R4b": r4b, "R4v": r4v, "R5": r5}[q](a)
        except Exception as e:
            import traceback
            traceback.print_exc()
            v = {"error": repr(e)}
        v["wall_s"] = time.perf_counter() - t0
        out[q] = v
        print(q, json.dumps(v, default=str)[:600], flush=True)
        a.out.write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
