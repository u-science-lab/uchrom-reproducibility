"""Reference answers of the four anchors of Part A (A1, A6, H1, S3: reading the data), the functions that computed
them, unchanged (run with the agents' first environment, U-Chrom fb43ac4, on a compute cluster).

    python supp_fig7/v3/truth/anchors.py --data-root <agentbench/data> [--out supp_fig7/v3/answers_anchors.json]

Writes ``{id: {"values": ..., "truth": fn}}``; ``truth/imgfmt.py`` copies them into ``results/imgfmt.json``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from chromdata import ChromData

HERE = Path(__file__).resolve().parent
DATA = {"takei": "takei2025_fov0/takei2025_fov0.chromdata.zarr", "schicar": "schicar_mop/schicar_mop.chromdata.zarr",
        "melanoma": "hicrna_human_melanoma/hicrna_human_melanoma.chromdata.zarr"}
ATLAS_HIRES = "https://uchrom-atlas-r2.u-science.org/hires_brain.chromdata.zarr"
ANCHORS = {"A1": "counts", "A6": "schicar_top_type", "H1": "atlas_hires_counts", "S3": "melanoma_dosage"}
_CACHE: dict = {}


def load(root: Path, key: str) -> ChromData:
    if key not in _CACHE:
        _CACHE[key] = ChromData.read(str(root / DATA[key]))
    return _CACHE[key]

def schicar_path(root, rel):
    return root / "schicar_mop" / rel

def _section_map(cd, section: str) -> str:
    from uchrom.fea.separation import _per_cell_map_keys

    keys = [k for k in _per_cell_map_keys(cd) if section in k]
    if len(keys) != 1:
        raise SystemExit(f"per-spot map of {section}: {keys}")
    return keys[0]

def counts(root):
    cd = load(root, "takei")
    return {"n_cells": int(cd.n_cells), "n_traces": int(cd.n_traces), "n_spots": int(cd.n_spots)}

def schicar_top_type(root):
    import h5py
    cd = load(root, "schicar")
    c = cd.cells
    vc = c["cell_type"].astype(str).value_counts()
    top = str(vc.index[0])
    sel = c["cell_type"].astype(str) == top
    # cross-check: contacts actually stored in the linked per-cell .scool (1 Mb, symmetric-upper, all chroms)
    path = schicar_path(root, cd.uns["linked_scool"]["per_cell"]["path"])
    with h5py.File(path, "r") as f:
        scool_sum = pd.Series({k: int(f["cells"][k].attrs["sum"]) for k in c.index[sel].astype(str)})
    return {"cell_type": top, "n": int(vc.iloc[0]), "median_contacts": float(c.loc[sel, "n_contacts_hic"].median()),
            "median_scool_pixel_sum": float(scool_sum.median()),
            "note": "median_contacts = cells['n_contacts_hic'] (per-cell Hi-C contact count stored with the cell "
                    "metadata); the per-cell .scool pixel sums (1 Mb, all chromosomes) are given for context "
                    "and differ because the .scool keeps only binned cis+trans contacts on the 21 chromosomes."}

def melanoma_dosage(root, min_contacts: int = 10_000):
    """S3: per-spot chromosome shares of the contact legs from the 1-Mb per-spot maps (cis 2 legs, trans 1 + 1)."""
    import tempfile

    import cooler
    from chromdata.embedded import export_scool

    store = root / DATA["melanoma"]
    cd = ChromData.read(str(store), backed=True)
    key = _section_map(cd, "human_melanoma_r1")
    with tempfile.TemporaryDirectory() as tmp:
        sc = Path(tmp) / "r1.scool"
        export_scool(str(store), key, sc, trans=True)
        names = list(cooler.fileops.list_coolers(str(sc)))
        legs = {}
        for n in names:
            clr = cooler.Cooler(f"{sc}::{n}")
            bins = clr.bins()[["chrom"]][:]
            px = clr.pixels()[:]
            ca = bins["chrom"].to_numpy()[px["bin1_id"].to_numpy()]
            cb = bins["chrom"].to_numpy()[px["bin2_id"].to_numpy()]
            cnt = px["count"].to_numpy(float)
            tot = pd.Series(cnt, index=ca).groupby(level=0).sum().add(pd.Series(cnt, index=cb).groupby(level=0).sum(),
                                                                         fill_value=0)
            if cnt.sum() >= min_contacts:
                legs[n.rsplit("/", 1)[-1]] = tot
    L = pd.DataFrame(legs).T.fillna(0.0)
    L = L.drop(columns=[c for c in L.columns if str(c) in ("chrY", "Y")], errors="ignore")
    share = L.div(L.sum(axis=1), axis=0)
    prof = np.log2(share / share.median())
    sd = prof.std().sort_values(ascending=False)
    return {"top_chrom": str(sd.index[0]), "sd_top": float(sd.iloc[0]), "n_spots": int(len(L)),
            "sd_by_chrom": {str(k): float(x) for k, x in sd.head(6).items()}}

def atlas_hires_counts(root):
    cd = ChromData.read(ATLAS_HIRES, backed=True)
    per_cell = cd.spots["cell_id"].astype(str).value_counts()
    return {"n_cells": int(cd.n_cells), "n_traces": int(cd.n_traces), "median_spots_per_cell": float(per_cell.median())}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=HERE.parent / "answers_anchors.json")
    a = ap.parse_args()
    out = {}
    for qid, fn in ANCHORS.items():
        out[qid] = {"values": globals()[fn](a.data_root), "truth": fn}
        print(qid, json.dumps(out[qid]["values"], default=str)[:300])
    a.out.write_text(json.dumps(out, indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
