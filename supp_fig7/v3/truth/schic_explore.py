"""Inspection of the inputs of the schic family (prints only; used to fix the questions' details)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(sys.argv[1])          # the shared data root


def main():
    import cooler
    from chromdata import ChromData

    cd = ChromData.read(str(root / "schicar_mop/schicar_mop.chromdata.zarr"))
    c = cd.cells
    print("cells", c.shape, list(c.columns)[:60])
    print(c["cell_type"].astype(str).value_counts().to_string())
    for col in c.columns:
        if c[col].dtype == object or str(c[col].dtype) == "category":
            if c[col].nunique() < 40:
                print("categorical column", col, c[col].nunique())
    print("uns linked_scool", cd.uns.get("linked_scool"))
    print("uns linked_cool keys", list((cd.uns.get("linked_cool") or {}).keys())[:30])
    scool = root / "schicar_mop/contacts_cells_1Mb.scool"
    names = cooler.fileops.list_scool_cells(str(scool))
    print("scool cells", len(names), names[:3])
    clr = cooler.Cooler(f"{scool}::{names[0]}")
    print("scool binsize", clr.binsize, "chroms", list(clr.chromnames), "n_bins", clr.info["nbins"])
    ids = set(map(str, c.index))
    print("scool names in cells:", sum(n.split("/")[-1] in ids for n in names))
    mc = root / "schicar_mop/contacts_all.mcool"
    print("mcool resolutions", cooler.fileops.list_coolers(str(mc)))
    for t in ("Oligo", "L6CT"):
        p = root / f"schicar_mop/contacts_{t}.cool"
        if p.exists():
            cc = cooler.Cooler(str(p))
            print(t, "cool binsize", cc.binsize, "weights", "weight" in cc.bins().columns)
    # chr1 cis contacts per cell (off-diagonal and all)
    import h5py
    rows = []
    with h5py.File(scool, "r") as f:
        chrom = f["bins/chrom"][:]
        chr1 = np.flatnonzero(chrom == 0)
        lo, hi = int(chr1[0]), int(chr1[-1]) + 1
        for n in names:
            g = f["cells"][n.split("/")[-1]]
            off = g["indexes/bin1_offset"]
            a, b = int(off[lo]), int(off[hi])
            b1 = g["pixels/bin1_id"][a:b]
            b2 = g["pixels/bin2_id"][a:b]
            cnt = g["pixels/count"][a:b]
            m = b2 < hi
            rows.append((n.split("/")[-1], int(cnt[m].sum()), int(cnt[m & (b1 != b2)].sum()), int(m.sum())))
    t = pd.DataFrame(rows, columns=["cell", "chr1_cis", "chr1_cis_offdiag", "chr1_pixels"]).set_index("cell")
    t["cell_type"] = c["cell_type"].astype(str).reindex(t.index)
    print(t.describe().to_string())
    print(t.groupby("cell_type")["chr1_cis"].median().to_string())
    for ty in ("Oligo", "L6CT", "L23IT.1", "Astro"):
        s = t[t["cell_type"] == ty].sort_values("chr1_cis")
        if len(s):
            mid = s.iloc[len(s) // 2]
            print("median cell of", ty, s.index[len(s) // 2], mid.to_dict())
    if "n_contacts_hic" in c:
        print("n_contacts_hic", c["n_contacts_hic"].describe().to_string())


if __name__ == "__main__":
    main()
