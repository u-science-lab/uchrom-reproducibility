"""Data folders of the version-3 questions of family ``imgfmt`` (T6, T7, F2, F3), from registered sources
(``v3/data/imgfmt.md`` describes each folder; F1, F4 read public atlas stores over HTTP and need no folder).

    python imgfmt_build.py sox2     --dl DL --out DATA_V3/imgfmt_sox2 [--n-hap N]                    (any env)
    python imgfmt_build.py seqfish  --shared SHARED --out DATA_V3/imgfmt_seqfish                     (any env)
    python imgfmt_build.py outlines --ifzip ZIP --shared SHARED --fov 4 --out OUTLINES.json          (venv with shapely)
    python imgfmt_build.py store    --shared SHARED --outlines OUTLINES.json --fov 4 --out STORE     (U-Chrom env)

``sox2``     Huang et al. 2021 Sox2 tables of SnapFISH-IMPUTE (source ``huang2021_sox2``, downloaded + md5-checked
             here) -> the 129 allele, a fixed hold-out of 10 % of its observed loci (``imgfmt_t6.make_input``).
``seqfish``  Takei et al. 2021 raw DNA seqFISH+ detections, 25-kb loci, replicate 1 (source ``seqfish``,
             DNAseqFISH+.zip) with the mm10 locus of every probe from the traced FOF-CT table of the same experiment
             (source ``takei``, 4DNFIHF3JCBY) -> detections + the authors' homolog labels (``imgfmt_t7.make_input``).
``outlines`` Nuclear outlines of the cells of one field of view of the same experiment from the authors' DAPI
             segmentation (Zenodo 3735329 ``IF-binned1x1x1-E14-replicate1-Pos<fov>.zip``: per cell, the voxels of the
             segmented nucleus, ``...-Ch4-Pos<fov>-Cell<k>-...csv``): the x-y footprint of the voxels (union over z),
             holes filled, as the polygon of its pixel edges in um (103 nm pixels); the largest part when the footprint
             is not connected.  Checked against the DNA detections of the same cells (share of spots inside).
``store``    The traced FOF-CT table (4DNFIHF3JCBY) + its cell table (4DNFIFINA2U9, source ``takei_tables``),
             cells of that field of view only, the outlines as ``cell_shapes['nucleus']`` and two results of
             ``uchrom.fea.separation_profile`` with provenance (questions F2, F3) -> ``takei2021_fov<k>.chromdata.zarr``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "truth"))

SOX2_COMMIT = "7f7f1a703aaf6db98147183588a9da6ee1d1c7be"
SOX2_MD5 = {"mESC_Sox2_coor_wnan.txt": "ab391637b75e412eec39d4082e6e3f8d",
            "mESC_Sox2_ann.txt": "c28b7b7e0f590a424f255388997f54bc"}
SEQFISH_MD5 = "c4e7b23df669bbc7ef9b1436712b4a2c"
TAKEI_MD5 = {"4DNFIHF3JCBY.csv": "c2f5a8db66bfc5ba9cbe6870b415aaa1", "4DNFIFINA2U9.csv": "799c25d9bf53ab0104911b47946d9920"}
IF_ZIP = {4: ("IF-binned1x1x1-E14-replicate1-Pos4.zip", "21e25a076802e8a89ce5dace95f5c6ce", 2471856846)}
PIXEL_UM = 0.103
#: the results stored in the F2/F3 store: key -> SeparationParams fields (bp, um)
RESULTS = {"separation": dict(short=(50_000.0, 150_000.0), long=(500_000.0, 1_000_000.0), stat="median",
                              duplicates="mean"),
           "separation_mean": dict(short=(25_000.0, 100_000.0), long=(750_000.0, 1_500_000.0), stat="mean",
                                   duplicates="mean")}


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def check_md5(path: Path, want: str):
    got = md5(path)
    if got != want:
        raise IOError(f"md5 mismatch for {path}: {got} != {want}")
    print(f"md5 ok {path.name} {got}", flush=True)


# ---------------------------------------------------------------------------------------------------- sox2
def build_sox2(dl: Path, out: Path, n_hap):
    import imgfmt_t6
    src = dl / "huang2021_sox2"
    src.mkdir(parents=True, exist_ok=True)
    for name, want in SOX2_MD5.items():
        p = src / name
        if not p.exists():
            url = f"https://raw.githubusercontent.com/hyuyu104/SnapFISH-IMPUTE/{SOX2_COMMIT}/data/{name}"
            urllib.request.urlretrieve(url, p)
        check_md5(p, want)
    info = imgfmt_t6.make_input(src, out, region="129", n_hap=n_hap, frac=0.1, seed=0)
    print(json.dumps(info, indent=1))


# ------------------------------------------------------------------------------------------------- seqfish
def build_seqfish(shared: Path, out: Path, chroms):
    import imgfmt_t7
    check_md5(shared / "DNAseqFISH+.zip", SEQFISH_MD5)
    check_md5(shared / "4DNFIHF3JCBY.csv", TAKEI_MD5["4DNFIHF3JCBY.csv"])
    info = imgfmt_t7.make_input(shared / "DNAseqFISH+.zip", shared / "4DNFIHF3JCBY.csv", out, chroms)
    print(json.dumps(info, indent=1))


# ------------------------------------------------------------------------------------------------ outlines
class HTTPRangeFile:
    """A read-only, seekable file over HTTP range requests (Zenodo serves the IF zips at ~0.3 MB/s, so only the
    central directory and the members we need are fetched: ~60 MB instead of 2.47 GB).  Every member read through
    ``zipfile`` is checked against the CRC-32 of the archive's central directory."""

    def __init__(self, url: str, size: int, block: int = 1 << 20):
        self.url, self.size, self.pos, self.block = url, size, 0, block
        self.fetched = 0
        self._cache = {}

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def _get(self, start: int, end: int) -> bytes:
        for attempt in range(6):
            try:
                req = urllib.request.Request(self.url, headers={"Range": f"bytes={start}-{end - 1}",
                                                                "User-Agent": "curl/8.4.0"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    data = r.read()
                if len(data) != end - start:
                    raise IOError(f"short read {len(data)} != {end - start}")
                self.fetched += len(data)
                return data
            except Exception as exc:          # transient network errors: retry
                print(f"range {start}-{end}: {exc}; retry {attempt + 1}", flush=True)
                time.sleep(5 * (attempt + 1))
        raise IOError(f"range {start}-{end} failed")

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        out = bytearray()
        while n > 0:
            b0 = (self.pos // self.block) * self.block
            if b0 not in self._cache:
                if len(self._cache) > 64:
                    self._cache.clear()
                self._cache[b0] = self._get(b0, min(b0 + self.block, self.size))
            chunk = self._cache[b0][self.pos - b0:self.pos - b0 + n]
            out += chunk
            self.pos += len(chunk)
            n -= len(chunk)
        return bytes(out)


def build_outlines(ifzip, shared: Path, fov: int, out: Path):
    import shapely
    from shapely.geometry import MultiPolygon, Polygon, box
    from shapely.ops import unary_union

    name, want, size = IF_ZIP[fov]
    if str(ifzip).startswith("http"):
        src = HTTPRangeFile(str(ifzip), size)
        integrity = "HTTP range reads of the members; each checked against the archive's CRC-32 by zipfile"
    else:
        check_md5(Path(ifzip), want)
        src = ifzip
        integrity = f"md5 {want}"
    rows, rings = [], {}
    with zipfile.ZipFile(src) as zf:
        members = [m for m in zf.namelist() if re.search(rf"-Ch4-Pos{fov}-Cell\d+-", m)]
        for m in sorted(members, key=lambda s: int(re.search(r"Cell(\d+)-", s).group(1))):
            k = int(re.search(r"Cell(\d+)-", m).group(1))
            with zf.open(m) as fh:
                vox = pd.read_csv(fh, usecols=["hybID", "xIdx", "yIdx", "zIdx"])
            vox = vox[vox["hybID"] == vox["hybID"].min()]
            px = np.unique(vox[["xIdx", "yIdx"]].to_numpy().round().astype(np.int64), axis=0)
            # pixel i (1-based, MATLAB) covers [i - 1, i) in pixel units
            geom = unary_union([box(x - 1, y - 1, x, y) for x, y in px])
            parts = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
            main = max(parts, key=lambda g: g.area)
            filled = Polygon(main.exterior)
            ring = np.asarray(filled.exterior.coords) * PIXEL_UM
            cid = f"{fov}_{k}"
            rings[cid] = ring.tolist()
            rows.append({"cell_id": cid, "n_voxels": int(len(vox)), "n_pixels_xy": int(len(px)),
                         "n_parts": len(parts), "area_px": float(filled.area), "area_um2": float(filled.area) * PIXEL_UM ** 2,
                         "n_vertices": int(len(ring)), "member": m})
    # check against the DNA detections of the same field of view (25-kb replicate 1, pixels)
    with zipfile.ZipFile(shared / "DNAseqFISH+.zip") as zf, \
            zf.open("DNAseqFISH+/DNAseqFISH+25kbloci-E14-replicate1.csv") as fh:
        det = pd.read_csv(fh, usecols=["fov", "cellID", "x", "y"])
    det = det[det["fov"] == fov]
    check = {}
    for orient in ("xy", "yx"):
        inside, n = 0, 0
        for r in rows:
            k = int(r["cell_id"].split("_")[1])
            d = det[det["cellID"] == k]
            if not len(d):
                continue
            poly = Polygon(np.asarray(rings[r["cell_id"]]) / PIXEL_UM)
            px, py = (d["x"], d["y"]) if orient == "xy" else (d["y"], d["x"])
            inside += int(shapely.contains_xy(poly, px.to_numpy(), py.to_numpy()).sum())
            n += len(d)
        check[orient] = {"spots": n, "inside_own_outline": inside / n if n else None}
    det_cells = sorted(det["cellID"].unique().tolist())
    info = {"zip": name, "zenodo_md5": want, "size": size, "integrity": integrity,
            "bytes_fetched": getattr(src, "fetched", None), "fov": fov, "n_outlines": len(rows), "cells": rows,
            "dna_cells_in_fov": det_cells, "dna_check": check}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"info": info, "rings_um": rings}, indent=1))
    print(json.dumps({k: v for k, v in info.items() if k != "cells"}, indent=1))
    print(pd.DataFrame(rows).drop(columns="member").to_string())


# --------------------------------------------------------------------------------------------------- store
def build_store(shared: Path, outlines: Path, fov: int, out: Path):
    from chromdata import ChromData
    from uchrom.fea import SeparationParams, separation_profile

    check_md5(shared / "4DNFIHF3JCBY.csv", TAKEI_MD5["4DNFIHF3JCBY.csv"])
    check_md5(shared / "4DNFIFINA2U9.csv", TAKEI_MD5["4DNFIFINA2U9.csv"])
    cd = ChromData.from_fofct(shared / "4DNFIHF3JCBY.csv", cell_table=shared / "4DNFIFINA2U9.csv")
    keep = [c for c in cd.cells.index.astype(str) if c.startswith(f"{fov}_")]
    sub = cd[np.flatnonzero(cd.spots["cell_id"].astype(str).isin(keep).to_numpy())]
    doc = json.loads(outlines.read_text())
    rings = doc["rings_um"]
    ids = [c for c in rings if c in set(keep)]
    sub.set_cell_shapes([np.asarray(rings[c]) for c in ids], key="nucleus", cell_ids=ids, unit="um", frame="fov",
                        in_coords_frame=bool(doc["info"]["dna_check"]["xy"]["inside_own_outline"] > 0.9),
                        status="measured",
                        source="Takei et al. 2021, Zenodo 3735329 " + doc["info"]["zip"] + ": DAPI segmentation "
                               "voxels of each nucleus (Ch4), x-y footprint, holes filled, pixel-edge polygon",
                        method="union of the 103-nm pixels of the segmented voxels over z")
    for key, kw in RESULTS.items():
        separation_profile(sub, params=SeparationParams(**kw), key_added=key)
    sub.uns["source"] = (f"Takei et al. 2021, Nature 590:344: DNA seqFISH+ of E14 mouse ES cells, replicate 1, field of "
                         f"view {fov} (4DN 4DNFIHF3JCBY core table + 4DNFIFINA2U9 cell table); nuclear outlines "
                         f"from the authors' DAPI segmentation (Zenodo 3735329)")
    if out.exists():
        import shutil
        shutil.rmtree(out)
    sub.write(out)
    back = ChromData.read(out)
    info = {"store": str(out), "n_cells": int(back.n_cells), "n_traces": int(back.n_traces),
            "n_spots": int(back.n_spots), "n_outlines": int(len(back.cell_shapes["nucleus"])),
            "results": {k: {"function": back.results.record(k).function, "params": back.results.record(k).params}
                        for k in back.results.keys()},
            "dup_spots_per_trace_locus": int(back.spots.duplicated(["trace_id", "bin_id"]).sum()),
            "cells_columns": [str(c) for c in back.cells.columns]}
    # the outline areas against the authors' nuclear area of the cell table (4DNFIFINA2U9 "area(um2)")
    areas = {r["cell_id"]: r["area_um2"] for r in doc["info"]["cells"]}
    col = next((c for c in back.cells.columns if "area" in str(c).lower() and "cyto" not in str(c).lower()), None)
    if col is not None:
        a = pd.DataFrame({"outline": pd.Series(areas), "authors": back.cells[col].reindex(list(areas))})
        info["area_check"] = {"column": str(col), "pearson": float(a.corr().iloc[0, 1]),
                              "median_ratio_outline_over_authors": float((a["outline"] / a["authors"]).median())}
    print(json.dumps(info, indent=1, default=str))
    (out.parent / "build_info.json").write_text(json.dumps(info, indent=1, default=str))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["sox2", "seqfish", "outlines", "store"])
    ap.add_argument("--dl", type=Path)
    ap.add_argument("--shared", type=Path)
    ap.add_argument("--ifzip", type=str)
    ap.add_argument("--outlines", type=Path)
    ap.add_argument("--fov", type=int, default=4)
    ap.add_argument("--n-hap", type=int, default=None)
    ap.add_argument("--chroms", default="chr3,chr7,chr11")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    t0 = time.perf_counter()
    if a.cmd == "sox2":
        build_sox2(a.dl, a.out, a.n_hap)
    elif a.cmd == "seqfish":
        build_seqfish(a.shared, a.out, a.chroms.split(","))
    elif a.cmd == "outlines":
        build_outlines(a.ifzip, a.shared, a.fov, a.out)
    else:
        build_store(a.shared, a.outlines, a.fov, a.out)
    print(f"{a.cmd}: {time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()
