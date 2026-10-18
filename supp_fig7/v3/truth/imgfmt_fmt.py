"""Questions F1-F4 (family imgfmt, Part A: reading the open store), each computed twice: ``chromdata`` (the U-Chrom
environment: chromdata / uchrom APIs) and ``raw`` (no chromdata: zarr / pyarrow / fsspec / numpy only, following
``packages/chromdata/spec.md``; run in the agents' plain environment).

    python imgfmt_fmt.py f1 --impl chromdata|raw --url URL --cell CELL --chrom chr7 --bin1 30 --bin2 31
    python imgfmt_fmt.py f2 --impl chromdata|raw --store STORE --key nucleus --cell CELL
    python imgfmt_fmt.py f3 --impl uchrom|raw   --store STORE --key separation --cell CELL --chrom chr7
    python imgfmt_fmt.py f4 --impl chromdata|raw --url URL --chrom chr19
    (every command: --out JSON)

Each prints and writes a JSON with the values, the seconds and the peak RSS of the process.
"""

from __future__ import annotations

import argparse
import json
import resource
import struct
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def peak_rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0      # Linux: kB


# ================================================================================================ raw helpers
def _store(url_or_path: str):
    import zarr
    if "://" in url_or_path:
        return zarr.storage.FsspecStore.from_url(url_or_path, read_only=True)
    return zarr.storage.LocalStore(url_or_path, read_only=True)


def _open_file(url_or_path: str, rel: str):
    if "://" in url_or_path:
        import fsspec
        return fsspec.open(f"{url_or_path.rstrip('/')}/{rel}", "rb").open()
    return open(Path(url_or_path) / rel, "rb")


def _read_parquet(url_or_path: str, rel: str, columns=None) -> pd.DataFrame:
    import pyarrow.parquet as pq
    with _open_file(url_or_path, rel) as fh:
        return pq.read_table(fh, columns=columns).to_pandas()


def _root_attrs(url_or_path: str) -> dict:
    with _open_file(url_or_path, "zarr.json") as fh:
        return json.load(fh)


# ============================================================================== F1: embedded per-cell contacts
def f1_stats(pix: np.ndarray, b1: int, b2: int) -> dict:
    """pix: (n, 3) chromosome-local bin1, bin2, count (upper triangle, bin1 <= bin2)."""
    total = int(pix[:, 2].sum())
    pixel = int(pix[(pix[:, 0] == b1) & (pix[:, 1] == b2), 2].sum())
    row = int(pix[(pix[:, 0] == b1) | (pix[:, 1] == b1), 2].sum())      # symmetric row, diagonal once
    diag = int(pix[pix[:, 0] == pix[:, 1], 2].sum())
    return {"total_cis": total, "pixel": pixel, "row_sum": row, "n_pixels": int(len(pix)), "diagonal": diag,
            "upper_triangle_ok": bool((pix[:, 0] <= pix[:, 1]).all())}


def f1_chromdata(url, cell, chrom, b1, b2, key="per_cell"):
    from chromdata.embedded import open_all
    e = open_all(url)[key]
    idx = e.cells.index(cell)
    pix = e.cell_pixels(f"cis/{chrom}", [idx])
    out = f1_stats(pix, b1, b2)
    out["cell_index"] = idx
    out["resolution"] = e.resolutions[0]
    return out


def f1_pick(url, cell_type, key="per_cell"):
    """Exploration: the cell of ``cell_type`` whose n_contacts is closest to the type's median; where it sits in
    the cells table and in the embedded map's cell list."""
    import zarr
    cells = _read_parquet(url, "tables/cells.parquet")
    if "cell_id" in cells.columns:
        cells = cells.set_index("cell_id")
    cells.index = cells.index.astype(str)
    g = zarr.open_group(_store(url), path=f"embedded/contacts/{key}", mode="r")
    emb = [str(c) for c in g["cells"][:]]
    sub = cells[cells["cell_type"] == cell_type]
    med = sub["n_contacts"].median()
    pick = (sub["n_contacts"] - med).abs().sort_values(kind="mergesort").index[0]
    same_order = emb == list(cells.index)
    return {"cell": pick, "n_contacts": int(sub.loc[pick, "n_contacts"]), "type_median": float(med),
            "n_type": int(len(sub)), "row_in_cells_table": int(list(cells.index).index(pick)),
            "index_in_embedded_cells": emb.index(pick), "same_order": same_order, "n_embedded": len(emb),
            "first_cells_table": list(cells.index[:3]), "first_embedded": emb[:3]}


def f1_raw(url, cell, chrom, b1, b2, key="per_cell"):
    import zarr
    st = _store(url)
    g = zarr.open_group(st, path=f"embedded/contacts/{key}", mode="r")
    meta = dict(g.attrs)["uchrom_contacts"]
    res = int(meta["resolutions"][0])
    cells = [str(c) for c in g["cells"][:]]
    idx = cells.index(cell)
    off = g[f"r{res}/cis/{chrom}/cell_offset"][idx:idx + 2].astype(np.int64)
    pix = g[f"r{res}/cis/{chrom}/pix"][:, int(off[0]):int(off[1])].T.astype(np.int64)
    out = f1_stats(pix, b1, b2)
    out["cell_index"] = idx
    out["resolution"] = res
    return out


# ================================================================================================ F2: outlines
def wkb_polygon_rings(buf: bytes):
    """Rings (lists of (x, y[, z]) arrays) of a WKB Polygon / Polygon Z (ISO 1003 or EWKB flag), no library."""
    endian = "<" if buf[0] == 1 else ">"
    (gtype,) = struct.unpack(endian + "I", buf[1:5])
    has_z = bool(gtype & 0x80000000) or (gtype // 1000) in (1, 3)
    base = (gtype & 0xFFFF) % 1000
    if base != 3:
        raise ValueError(f"not a polygon (WKB type {gtype})")
    dim = 3 if has_z else 2
    pos = 5
    (n_rings,) = struct.unpack(endian + "I", buf[pos:pos + 4]); pos += 4
    rings = []
    for _ in range(n_rings):
        (n,) = struct.unpack(endian + "I", buf[pos:pos + 4]); pos += 4
        arr = np.frombuffer(buf[pos:pos + 8 * dim * n], dtype=endian + "f8").reshape(n, dim)
        pos += 8 * dim * n
        rings.append(arr)
    return rings


def shoelace(ring: np.ndarray) -> float:
    x, y = ring[:, 0], ring[:, 1]
    return float(abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / 2.0)


def polygon_area(rings) -> float:
    return shoelace(rings[0]) - sum(shoelace(r) for r in rings[1:])


def f2_summary(areas: pd.Series, cell: str, unit) -> dict:
    return {"n_outlines": int(len(areas)), "area": float(areas.loc[cell]), "unit": unit,
            "largest_cell": str(areas.idxmax()), "largest_area": float(areas.max()),
            "median_area": float(areas.median())}


def f2_chromdata(store, key, cell):
    from chromdata import ChromData
    from chromdata.cellspatial import polygon_rings
    cd = ChromData.read(store, backed=True)
    shp = cd.cell_shapes[key]
    areas = pd.Series({str(c): polygon_area([np.asarray(r) for r in polygon_rings(g)])
                       for c, g in zip(shp["cell_id"], shp["geometry"])})
    rec = (cd.uns.get("cell_shapes") or {}).get(key, {})
    out = f2_summary(areas, cell, rec.get("unit"))
    out["n_cells_in_store"] = int(cd.n_cells)
    return out


def f2_raw(store, key, cell):
    tab = _read_parquet(store, f"tables/cell_shapes/{key}.parquet")
    areas = pd.Series({str(c): polygon_area(wkb_polygon_rings(bytes(g)))
                       for c, g in zip(tab["cell_id"], tab["geometry"])})
    # the unit: uns record (spec: uns/ attrs {values: JSON})
    with _open_file(store, "uns/zarr.json") as fh:
        uns = json.load(fh)["attributes"]["values"]
    rec = (uns.get("cell_shapes") or {}).get(key, {})
    cells = _read_parquet(store, "tables/cells.parquet")
    out = f2_summary(areas, cell, rec.get("unit") if isinstance(rec, dict) else None)
    out["n_cells_in_store"] = int(len(cells))
    return out


# ============================================================================================ F3: provenance
def f3_uchrom(store, key, cell, chrom):
    from chromdata import ChromData
    from uchrom.fea import SeparationParams, separation_profile
    cd = ChromData.read(store)
    rec = cd.results.record(key)
    tab = rec.value.set_index("cell_id") if "cell_id" in rec.value.columns else rec.value
    p = dict(rec.params)
    fields = set(SeparationParams.__dataclass_fields__)
    kw = {k: (tuple(v) if k in ("short", "long") else v) for k, v in p.items() if k in fields}
    kw["chroms"] = [chrom]
    t = separation_profile(cd, cells=[cell], params=SeparationParams(**kw), key_added=None)
    return {"function": rec.function, "params": p, "uchrom_version": rec.uchrom_version,
            "stored_d_short": float(tab.loc[cell, "d_short"]), "stored_d_long": float(tab.loc[cell, "d_long"]),
            "stored_n_pairs_short": int(tab.loc[cell, "n_pairs_short"]),
            "recomputed_d_short": float(t.loc[cell, "d_short"]), "recomputed_d_long": float(t.loc[cell, "d_long"]),
            "recomputed_n_pairs_short": int(t.loc[cell, "n_pairs_short"]),
            "recomputed_n_pairs_long": int(t.loc[cell, "n_pairs_long"])}


def _categories(store, rel):
    df = _read_parquet(store, rel)
    return df.iloc[:, 0].astype(str).to_numpy()


def separation_numpy(coords: pd.DataFrame, bins: pd.DataFrame, p: dict, chroms) -> dict:
    """coords: bin_id, trace (any id), x, y, z of ONE cell; bins indexed by bin_id with chrom, start, end."""
    mid = ((bins["start"] + bins["end"]) / 2.0)
    df = coords.copy()
    df["chrom"] = bins.loc[df["bin_id"], "chrom"].astype(str).to_numpy()
    df = df[df["chrom"].isin(chroms)]
    dup = p.get("duplicates", "last")
    if dup in ("mean", "median"):
        df = df.groupby(["trace", "bin_id"], as_index=False)[["x", "y", "z"]].agg(dup)
    elif dup in ("first", "last"):
        df = df.drop_duplicates(["trace", "bin_id"], keep=dup)
    df = df[np.isfinite(df[["x", "y", "z"]]).all(axis=1)]
    out = {}
    for name in ("short", "long"):
        lo, hi = (float(v) for v in p[name])
        d_all = []
        for _, g in df.groupby("trace"):
            xyz = g[["x", "y", "z"]].to_numpy(float)
            m = mid.loc[g["bin_id"]].to_numpy(float)
            i, j = np.triu_indices(len(g), 1)
            s = np.abs(m[j] - m[i])
            sel = (s >= lo) & (s < hi)
            d_all.append(np.linalg.norm(xyz[i[sel]] - xyz[j[sel]], axis=1))
        d = np.concatenate(d_all) if d_all else np.zeros(0)
        stat = np.median if p.get("stat", "median") == "median" else np.mean
        out[f"d_{name}"] = float(stat(d)) if len(d) else float("nan")
        out[f"n_pairs_{name}"] = int(len(d))
    return out


def f3_raw(store, key, cell, chrom):
    import zarr
    g = zarr.open_group(_store(store), path="results", mode="r")
    node = dict(g.attrs)["keys"][key]
    attrs = dict(g[node].attrs)
    tab = _read_parquet(store, f"results/{node}/table.parquet")
    tab = tab.set_index("cell_id") if "cell_id" in tab.columns else tab
    p = attrs["_params"]
    bins = _read_parquet(store, "tables/bins.parquet")
    if "bin_id" in bins.columns:
        bins = bins.set_index("bin_id")
    cells_cat = _categories(store, "tables/categories/spots.cell_id.parquet")
    code = int(np.flatnonzero(cells_cat == cell)[0])
    tcols = ["bin_id", "trace_id", "cell_id", "x", "y", "z"]
    co = _read_parquet(store, f"tables/coords/chrom={chrom}/part-0.parquet", columns=tcols)
    co = co[co["cell_id"] == code].rename(columns={"trace_id": "trace"})
    rec = separation_numpy(co, bins, p, [chrom])
    return {"function": attrs.get("_function"), "params": p, "uchrom_version": attrs.get("_uchrom_version"),
            "stored_d_short": float(tab.loc[cell, "d_short"]), "stored_d_long": float(tab.loc[cell, "d_long"]),
            "stored_n_pairs_short": int(tab.loc[cell, "n_pairs_short"]),
            "recomputed_d_short": rec["d_short"], "recomputed_d_long": rec["d_long"],
            "recomputed_n_pairs_short": rec["n_pairs_short"], "recomputed_n_pairs_long": rec["n_pairs_long"]}


# ====================================================================================== F4: large store, HTTP
def rg_by_trace(trace: np.ndarray, xyz: np.ndarray) -> pd.Series:
    df = pd.DataFrame({"t": trace, "x": xyz[:, 0], "y": xyz[:, 1], "z": xyz[:, 2]})
    df = df[np.isfinite(df[["x", "y", "z"]]).all(axis=1)]
    c = df.groupby("t")[["x", "y", "z"]].transform("mean")
    sq = ((df[["x", "y", "z"]] - c) ** 2).sum(axis=1)
    return np.sqrt(sq.groupby(df["t"]).mean())


def f4_summary(rg: pd.Series, cell_of_trace: pd.Series, cell_type: pd.Series, types) -> dict:
    df = pd.DataFrame({"rg": rg, "cell": cell_of_trace.reindex(rg.index).astype(str)})
    df["type"] = cell_type.reindex(df["cell"]).to_numpy()
    by = df.groupby("type")["rg"].agg(["median", "size"])
    by_cells = df.groupby("type")["cell"].nunique()
    out = {"n_traces": int(len(df)), "n_cells": int(df["cell"].nunique()),
           "median_rg_all": float(df["rg"].median()),
           "by_type": {str(t): {"median_rg": float(r["median"]), "n_traces": int(r["size"]),
                                "n_cells": int(by_cells[t])} for t, r in by.iterrows()}}
    big = by[by_cells.reindex(by.index) >= 100]
    out["largest_median_type_100cells"] = str(big["median"].idxmax())
    out["smallest_median_type_100cells"] = str(big["median"].idxmin())
    for t in types:
        out[f"median_rg[{t}]"] = float(by.loc[t, "median"]) if t in by.index else None
    return out


def f4_chromdata(url, chrom, types):
    from chromdata import ChromData
    cd = ChromData.read(url, backed=True)
    sub = cd.get_chrom(chrom, columns="coords")
    sp = sub.spots
    rg = rg_by_trace(sp["trace_id"].astype(str).to_numpy(), np.asarray(sub.coords, float))
    cell_of_trace = sp.groupby(sp["trace_id"].astype(str), observed=True)["cell_id"].first().astype(str)
    ct = cd.cells["cell_type"].copy()
    ct.index = ct.index.astype(str)
    out = f4_summary(rg, cell_of_trace, ct, types)
    out["xyz_unit"] = cd.uns.get("xyz_unit")
    return out


def f4_raw(url, chrom, types):
    tcols = ["trace_id", "cell_id", "x", "y", "z"]
    co = _read_parquet(url, f"tables/coords/chrom={chrom}/part-0.parquet", columns=tcols)
    traces = _categories(url, "tables/categories/spots.trace_id.parquet")
    cells = _categories(url, "tables/categories/spots.cell_id.parquet")
    rg = rg_by_trace(co["trace_id"].to_numpy(), co[["x", "y", "z"]].to_numpy(np.float64))
    cell_of_code = co.groupby("trace_id")["cell_id"].first()
    cell_of_trace = pd.Series(cells[cell_of_code.to_numpy()], index=cell_of_code.index)
    ct = _read_parquet(url, "tables/cells.parquet")
    if "cell_id" in ct.columns:
        ct = ct.set_index("cell_id")
    ct.index = ct.index.astype(str)
    out = f4_summary(rg, cell_of_trace, ct["cell_type"], types)
    out["n_trace_categories"] = int(len(traces))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("q", choices=["f1", "f2", "f3", "f4"])
    ap.add_argument("--impl", required=True)
    ap.add_argument("--url")
    ap.add_argument("--store")
    ap.add_argument("--cell")
    ap.add_argument("--chrom")
    ap.add_argument("--bin1", type=int)
    ap.add_argument("--bin2", type=int)
    ap.add_argument("--key")
    ap.add_argument("--types", default="")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    t0 = time.perf_counter()
    types = [t for t in a.types.split("|") if t]
    if a.q == "f1" and a.impl == "pick":
        val = f1_pick(a.url, types[0], a.key or "per_cell")
    elif a.q == "f1":
        fn = f1_chromdata if a.impl == "chromdata" else f1_raw
        val = fn(a.url, a.cell, a.chrom, a.bin1, a.bin2, a.key or "per_cell")
    elif a.q == "f2":
        val = (f2_chromdata if a.impl == "chromdata" else f2_raw)(a.store, a.key or "nucleus", a.cell)
    elif a.q == "f3":
        val = (f3_uchrom if a.impl == "uchrom" else f3_raw)(a.store, a.key or "separation", a.cell, a.chrom)
    else:
        val = (f4_chromdata if a.impl == "chromdata" else f4_raw)(a.url, a.chrom, types)
    res = {"question": a.q, "impl": a.impl, "values": val, "seconds": time.perf_counter() - t0,
           "peak_rss_mb": peak_rss_mb(), "python": sys.version.split()[0]}
    print(json.dumps(res, indent=1, default=str), flush=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
