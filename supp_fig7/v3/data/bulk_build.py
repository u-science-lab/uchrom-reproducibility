"""Data folders of the version-3 questions of family ``bulk`` (H1-H4, XM1, P2), built from registered sources
without U-Chrom's code (cooler / pandas / pyarrow / numpy only; run in the
benchmark's ``venv-uchrom``, which has pyarrow, but never imports ``uchrom``):

``rao2014_imr90_chr21_10kb/``  (data key ``rao_chr21_10kb``; questions H1, H2, H3)
    ``IMR90_chr21_10kb_hg19.cool``: the Rao et al. 2014 IMR90 chr21 5-kb slice of the GEO .hic (registered source
    ``rao2014_imr90_chr21``) summed into 10-kb bins (``cooler.coarsen_cooler``, factor 2) and ICE-balanced as
    ``cooler balance --cis-only`` with its defaults (weights in the column ``weight``);
    ``GSE63525_IMR90_HiCCUPS_looplist.txt.gz``: the authors' HiCCUPS loops (source ``rao2014_annotations``).
``su2020_chr21_tsv/``  (data key ``su2020_tsv``; questions H4, XM1)
    ``chromosome21.tsv``, ``Hi-C_contacts_chromosome21.tsv``, ``README_August_2020.txt``: the Su et al. 2020 files
    as published (source ``su2020_chr21``, md5-checked).
``takei2025_fov0_rna/``  (data key ``takei_rna``; question P2)
    ``takei2025_fov0_rna_counts.tsv``: RNA seqFISH transcripts per cell and gene of the 47 cells of the Takei 2025
    FOV-0 store (rows: ``cell_id`` as in the store, ``rep_fov_cellID``; columns: genes), counted from the mRNA spot
    tables of cerebellum replicate 1 in the Zenodo ``transcriptomic-data.zip`` (source ``takei2025_rna``,
    md5-checked; barcoded and non-barcoded mRNA panels).

    python bulk_build.py --shared /shared --data /data --out /data_v3 --dl /dl [--only rao,su,takei]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

TAKEI_ZIP = "https://zenodo.org/records/7693825/files/transcriptomic-data.zip"
TAKEI_MD5 = "c29594d9ed0f4dfea67c0140c56050a1"
SU_MD5 = {"chromosome21.tsv": "170d9d8b081edaa5f4ccaad8aad25eff",
          "Hi-C_contacts_chromosome21.tsv": "7842ef3149e139fe08c6de7df13b3084",
          "README_August_2020.txt": "a93629f38fcc11acdb46ced5fc74bee1"}
LOOPS_MD5 = {"GSE63525_IMR90_HiCCUPS_looplist.txt.gz": "326c6dbafcb860ad1e726c544f71a5a6"}


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def copy_checked(src: Path, dst: Path, digest: str):
    if dst.exists() and md5(dst) == digest:
        return
    got = md5(src)
    if got != digest:
        raise SystemExit(f"{src}: md5 {got} != {digest}")
    shutil.copyfile(src, f"{dst}.part")
    os.replace(f"{dst}.part", dst)


def build_rao(shared: Path, out: Path) -> dict:
    import cooler

    d = out / "rao2014_imr90_chr21_10kb"
    d.mkdir(parents=True, exist_ok=True)
    src = shared / "rao2014_imr90" / "IMR90_chr21_5kb_hg19.cool"
    dst = d / "IMR90_chr21_10kb_hg19.cool"
    info = {"source_5kb": str(src), "source_5kb_md5": md5(src)}
    if not dst.exists():
        tmp = d / "IMR90_chr21_10kb_hg19.cool.part"
        if tmp.exists():
            tmp.unlink()
        cooler.coarsen_cooler(str(src), str(tmp), factor=2, chunksize=10_000_000)
        clr = cooler.Cooler(str(tmp))
        _, stats = cooler.balance_cooler(clr, cis_only=True, store=True)   # = cooler balance --cis-only
        info["balance_converged"] = bool(np.all(stats["converged"]))
        info["balance_var"] = np.asarray(stats["var"], float).ravel().tolist()
        os.replace(tmp, dst)
    clr = cooler.Cooler(str(dst))
    w = clr.bins()[:]["weight"].to_numpy()
    info.update({"cool": str(dst), "binsize": int(clr.binsize), "chromnames": list(clr.chromnames),
                 "n_bins": int(clr.info["nbins"]), "n_pixels": int(clr.info["nnz"]),
                 "sum_counts": int(clr.pixels()[:]["count"].sum()), "n_bins_masked": int(np.isnan(w).sum()),
                 "balance_attrs": {k: (v if isinstance(v, (int, float, str, bool)) else str(v))
                                   for k, v in clr.open("r")["bins/weight"].attrs.items()}})
    for name, digest in LOOPS_MD5.items():
        copy_checked(shared / "rao2014_annotations" / name, d / name, digest)
    loops = pd.read_csv(d / "GSE63525_IMR90_HiCCUPS_looplist.txt.gz", sep="\t")
    info["published_loops_chr21"] = int((loops["chr1"].astype(str) == "21").sum())
    return info


def build_su(shared: Path, out: Path) -> dict:
    d = out / "su2020_chr21_tsv"
    d.mkdir(parents=True, exist_ok=True)
    for name, digest in SU_MD5.items():
        copy_checked(shared / "su2020_imr90" / name, d / name, digest)
    return {"files": sorted(SU_MD5), "md5": SU_MD5}


def _download(url: str, dst: Path, digest: str):
    if dst.exists() and md5(dst) == digest:
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(f"{dst}.part")
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, length=1 << 22)
    got = md5(tmp)
    if got != digest:
        raise SystemExit(f"{url}: md5 {got} != {digest}")
    os.replace(tmp, dst)


def build_takei(data: Path, out: Path, dl: Path) -> dict:
    """Counts of RNA seqFISH spots per cell and gene for the FOV-0 store's 47 cells (replicate 1, fov 0)."""
    zpath = dl / "takei2025_rna" / "transcriptomic-data.zip"
    t = time.time()
    _download(TAKEI_ZIP, zpath, TAKEI_MD5)
    info = {"zip": str(zpath), "zip_md5": TAKEI_MD5, "download_s": round(time.time() - t, 1)}
    cells = pd.read_parquet(data / "takei2025_fov0" / "takei2025_fov0.chromdata.zarr" / "tables" / "cells.parquet")
    ids = [str(c) for c in (cells["cell_id"] if "cell_id" in cells.columns else cells.index)]
    info["n_store_cells"] = len(ids)
    with zipfile.ZipFile(zpath) as zf:
        members = [m for m in zf.namelist() if m.lower().endswith(".csv")]
        info["zip_csv_members"] = members
        use = [m for m in members if "cerebellum" in m.lower() and "replicate1" in m.lower() and "mrna" in m.lower()]
        info["members_used"] = use
        counts = None
        for m in use:
            with zf.open(m) as fh:
                for chunk in pd.read_csv(fh, chunksize=1_000_000, low_memory=False):
                    if "rep" in chunk.columns and not (chunk["rep"].astype(int) == 1).all():
                        chunk = chunk[chunk["rep"].astype(int) == 1]
                    cid = "1_" + chunk["fov"].astype(int).astype(str) + "_" + chunk["cellID"].astype(int).astype(str)
                    keep = cid.isin(ids)
                    if not keep.any():
                        continue
                    g = pd.DataFrame({"cell_id": cid[keep].to_numpy(), "gene": chunk.loc[keep, "geneID"].astype(str)
                                      .to_numpy()}).groupby(["cell_id", "gene"]).size()
                    counts = g if counts is None else counts.add(g, fill_value=0)
    if counts is None:
        raise SystemExit("no RNA spot of the FOV-0 cells found")
    wide = counts.astype(np.int64).unstack("gene", fill_value=0).reindex(ids, fill_value=0)
    wide = wide[sorted(wide.columns)]
    wide.index.name = "cell_id"
    d = out / "takei2025_fov0_rna"
    d.mkdir(parents=True, exist_ok=True)
    wide.to_csv(d / "takei2025_fov0_rna_counts.tsv", sep="\t")
    tot = wide.sum(axis=1)
    info.update({"n_cells": int(len(wide)), "n_genes": int(wide.shape[1]), "n_cells_without_rna": int((tot == 0).sum()),
                 "total_transcripts_median": float(tot.median()), "total_transcripts_min": int(tot.min()),
                 "total_transcripts_max": int(tot.max()), "table": str(d / "takei2025_fov0_rna_counts.tsv"),
                 "table_md5": md5(d / "takei2025_fov0_rna_counts.tsv")})
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--shared", type=Path, required=True, help="the shared data directory ($SCRATCH/uchrom-nd/data)")
    ap.add_argument("--data", type=Path, required=True, help="the shared data root (Takei FOV-0 store)")
    ap.add_argument("--out", type=Path, required=True, help="the version-3 data root (data_v3)")
    ap.add_argument("--dl", type=Path, required=True, help="download cache (the RNA zip)")
    ap.add_argument("--only", default="rao,su,takei")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    report = {}
    for part in a.only.split(","):
        t = time.time()
        if part == "rao":
            report[part] = build_rao(a.shared, a.out)
        elif part == "su":
            report[part] = build_su(a.shared, a.out)
        elif part == "takei":
            report[part] = build_takei(a.data, a.out, a.dl)
        report[part]["seconds"] = round(time.time() - t, 1)
        print(part, json.dumps(report[part], default=str)[:2000], flush=True)
    dst = a.out / "bulk_build_report.json"
    old = json.loads(dst.read_text()) if dst.exists() else {}
    dst.write_text(json.dumps({**old, **report}, indent=1, default=str))
    if "uchrom" in sys.modules:
        raise SystemExit("bulk_build.py imported uchrom")


if __name__ == "__main__":
    main()
