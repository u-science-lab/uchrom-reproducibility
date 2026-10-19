"""The data of T1, P1 and P2 in community formats, next to the shared ChromData stores (``v3/data/neutral.md``),
so that a condition without U-Chrom does not have to read a ChromData store by hand.  Each export is checked against
its store value by value: the references do not change.

    python neutral_build.py takei    --store TAKEI_STORE --rna RNA_TSV --out DATA_V3/takei2025_fov0_fofct
    python neutral_build.py melanoma --store MELANOMA_STORE --section human_melanoma_r1 --out DATA_V3/hicrna_melanoma_r1
    (U-Chrom's environment; --results results/bulk.json to compare with the references of P1 / P2)

``takei``    Takei et al. 2025, field of view 0 (the shared store ``takei2025_fov0``) -> a 4DN FOF-CT core table
             (``chromdata.ChromData.to_fofct``: one row per spot, X / Y / Z in um, Trace_ID, Cell_ID, the locus;
             the store's spot tracks left out) and its FOF-CT cell table (every column of the store's ``cells``,
             ``cell_type`` among them).  Checked: the two tables read with pandas alone equal the store spot by spot
             and cell by cell; P2's statistic recomputed from the files (pandas / numpy) equals the reference.
``melanoma`` Guo et al. 2026 spatial Hi-C-RNA, one section (the atlas store ``hicrna_human_melanoma``) -> the
             per-spot 1-Mb maps of the section with the contacts between chromosomes as one ``.scool``
             (``chromdata.embedded.export_scool(..., trans=True)``, the read P1's reference makes) and a spot table
             (cell_id = the scool's cell name, array_row, array_col).  Checked: the legs per chromosome and the short /
             long counts of P1 read from the scool with cooler equal those read from the store's embedded map; P1's
             observed R^2 recomputed from the files equals the reference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "truth"))


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def read_fofct(path: Path, float_precision: str | None = "round_trip") -> tuple[pd.DataFrame, dict]:
    """A FOF-CT table with pandas alone: the ``#`` header lines, the columns of ``##columns=(...)``, the rows (decimals
    parsed exactly: pandas' default parser can be one unit in the last place off)."""
    header, n = {}, 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            n += 1
            header[line.lstrip("#").split("=", 1)[0].split(":", 1)[0].strip()] = line.rstrip("\n")
    cols = header["columns"].split("(", 1)[1].rstrip(")").split(",")
    return pd.read_csv(path, skiprows=n, header=None, names=[c.strip() for c in cols], low_memory=False,
                       float_precision=float_precision), header


def check(name: str, ok: bool, report: dict) -> None:
    report.setdefault("checks", {})[name] = bool(ok)
    print(f"{'ok  ' if ok else 'FAIL'} {name}", flush=True)


# ================================================================================================ takei (T1, P2)
def takei(a) -> dict:
    from chromdata import ChromData

    from bulk import MIN_SPOTS, ranks_partial

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    core, cells_csv = out / "takei2025_fov0_fofct_core.csv", out / "takei2025_fov0_fofct_cells.csv"
    cd = ChromData.read(a.store)
    print(f"store: {len(cd.spots)} spots, {len(cd.cells)} cells; "
          f"spot columns {list(cd.spots.columns)}; spot tracks {list(cd.spot_tracks.columns)}; "
          f"cells columns {list(cd.cells.columns)}; uns xyz_unit={cd.uns.get('xyz_unit')!r} "
          f"genome={cd.uns.get('genome_assembly')!r}", flush=True)
    cd.to_fofct(core, cell_table=cells_csv, spot_tracks=False)
    rep = {"files": {p.name: {"bytes": p.stat().st_size, "md5": md5(p)} for p in (core, cells_csv)}}

    t, head = read_fofct(core)
    c, chead = read_fofct(cells_csv)
    rep["core_header"], rep["cells_header"] = list(head.values()), list(chead.values())
    rep["core_columns"], rep["cells_columns"] = list(t.columns), list(c.columns)
    sp = cd.spots_with_loci()
    xyz = np.asarray(cd.coords)
    rep["store_coords_dtype"] = str(xyz.dtype)
    check("core: one row per spot", len(t) == len(sp), rep)
    check("core: X, Y, Z equal the store", np.array_equal(t[["X", "Y", "Z"]].to_numpy(xyz.dtype), xyz,
                                                         equal_nan=True), rep)
    # what pandas' default (fast) parser makes of the same decimals: within one unit in the last place
    t_fast, _ = read_fofct(core, float_precision=None)
    rep["default_parser_max_abs_diff_um"] = float(np.nanmax(np.abs(t_fast[["X", "Y", "Z"]].to_numpy(float) - xyz)))
    for col, ref in (("Trace_ID", sp["trace_id"]), ("Cell_ID", sp["cell_id"]), ("Chrom", sp["chrom"])):
        check(f"core: {col} equals the store", np.array_equal(t[col].astype(str).to_numpy(),
                                                              ref.astype(str).to_numpy()), rep)
    for col, ref in (("Chrom_Start", sp["start"]), ("Chrom_End", sp["end"])):
        check(f"core: {col} equals the store", np.array_equal(t[col].to_numpy(np.int64),
                                                              ref.to_numpy(np.int64)), rep)
    cells = cd.cells.copy()
    cells.index = cells.index.astype(str)
    check("cells: one row per cell, same ids", list(c["Cell_ID"].astype(str)) == list(cells.index), rep)
    check("cells: cell_type equals the store", list(c["cell_type"].astype(str)) ==
          list(cells["cell_type"].astype(str)), rep)

    # P2 from the files alone (the question's recipe, pandas / numpy)
    df = t[["Trace_ID", "Cell_ID", "X", "Y", "Z"]].astype({"Trace_ID": str, "Cell_ID": str})
    df = df[np.isfinite(df[["X", "Y", "Z"]].to_numpy()).all(1)]
    cen = df.groupby("Trace_ID")[["X", "Y", "Z"]].transform("mean")
    df = df.assign(r2=((df[["X", "Y", "Z"]] - cen) ** 2).sum(axis=1))
    per = df.groupby("Trace_ID").agg(n=("r2", "size"), ms=("r2", "mean"), cell=("Cell_ID", "first"))
    per = per[per["n"] >= MIN_SPOTS]
    med = np.sqrt(per["ms"]).groupby(per["cell"]).median()
    rna = pd.read_csv(a.rna, sep="\t", index_col=0)
    rna.index = rna.index.astype(str)
    ct = c.assign(Cell_ID=c["Cell_ID"].astype(str)).set_index("Cell_ID")["cell_type"].astype(str)
    ids = [i for i in ct.index if i in med.index]
    p2 = ranks_partial(rna.sum(axis=1).reindex(ids).fillna(0).to_numpy(float), med.reindex(ids).to_numpy(float),
                       ct.loc[ids].to_numpy())
    ref = json.load(open(a.results))["P2"]["values"] if a.results else {}
    rep["P2_from_files"], rep["P2_reference"] = p2, ref
    print(f"P2 from the files: {p2}; reference {ref}", flush=True)
    if ref:
        check("P2 from the files equals the reference",
              all(abs(p2[k] - ref[k]) < 1e-12 for k in ("spearman", "partial_spearman")), rep)
    return rep


# ================================================================================================ melanoma (P1)
def scool_legs(path: Path) -> tuple:
    """P1's per-spot tables from a .scool with cooler alone: legs per chromosome (cis 2, trans 1 + 1), the counts at
    bin offset 1 and 19-53 within autosomes."""
    import cooler

    names = list(cooler.fileops.list_coolers(str(path)))
    first = cooler.Cooler(f"{path}::{names[0]}")
    chrom = first.bins()["chrom"][:].astype(str).to_numpy()
    chroms = list(first.chromnames)
    auto = np.array([c.replace("chr", "").isdigit() for c in chrom])
    legs, short, long_ = {}, {}, {}
    for nm in names:
        px = cooler.Cooler(f"{path}::{nm}").pixels()[:]
        b1, b2 = px["bin1_id"].to_numpy(np.int64), px["bin2_id"].to_numpy(np.int64)
        cnt = px["count"].to_numpy(float)
        cid = nm.rsplit("/", 1)[-1]
        legs[cid] = (pd.Series(cnt, index=chrom[b1]).groupby(level=0).sum()
                     .add(pd.Series(cnt, index=chrom[b2]).groupby(level=0).sum(), fill_value=0))
        cis = (chrom[b1] == chrom[b2]) & auto[b1]
        off = b2 - b1
        short[cid] = float(cnt[cis & (off == 1)].sum())
        long_[cid] = float(cnt[cis & (off >= 19) & (off <= 53)].sum())
    L = pd.DataFrame(legs).T.reindex(columns=chroms).fillna(0.0)
    return L, pd.Series(short), pd.Series(long_), first.binsize


def melanoma(a) -> dict:
    from chromdata import ChromData
    from chromdata.embedded import export_scool

    from bulk import _cn_pcs, _melanoma_legs_fast, _r2, _section_key, _section_spots
    import bulk

    bulk.SECTION = a.section
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    sc = out / f"{a.section}.spots_1Mb.scool"
    spots_csv = out / f"{a.section}_spots.csv"
    cd = ChromData.read(a.store, backed=True)
    g = _section_spots(cd)
    key = _section_key(cd)
    print(f"section {a.section}: {len(g)} spots, map {key!r}; cells columns {list(cd.cells.columns)}", flush=True)
    t0 = time.perf_counter()
    if sc.exists():
        sc.unlink()
    export_scool(str(a.store), key, sc, cells=list(g.index), trans=True)
    t_export = time.perf_counter() - t0
    tab = pd.DataFrame({"cell_id": g.index.astype(str), "array_row": g["array_row"].to_numpy(int),
                        "array_col": g["array_col"].to_numpy(int)})
    tab.to_csv(spots_csv, index=False)
    rep = {"map": key, "n_spots": int(len(g)), "export_s": round(t_export, 1),
           "files": {p.name: {"bytes": p.stat().st_size, "md5": md5(p)} for p in (sc, spots_csv)}}

    t0 = time.perf_counter()
    L, short, long_, binsize = scool_legs(sc)
    rep["read_scool_s"] = round(time.perf_counter() - t0, 1)
    Ls, ss, ls = _melanoma_legs_fast(Path(a.store), key, g.index)
    for x in (Ls, ss, ls):
        x.index = x.index.astype(str)
    check("scool: 1-Mb bins", binsize == 1_000_000, rep)
    check("scool: one map per spot of the section", sorted(L.index) == sorted(g.index.astype(str)), rep)
    Lc = L.reindex(index=Ls.index, columns=Ls.columns).fillna(0.0)
    check("scool: legs per chromosome equal the store's map", np.array_equal(Lc.to_numpy(), Ls.to_numpy()), rep)
    check("scool: offset-1 counts equal the store's map", np.array_equal(short.reindex(ss.index).to_numpy(),
                                                                         ss.to_numpy()), rep)
    check("scool: offset 19-53 counts equal the store's map", np.array_equal(long_.reindex(ls.index).to_numpy(),
                                                                             ls.to_numpy()), rep)
    rep["contacts"] = {"legs": float(Lc.to_numpy().sum()), "trans_spots": int((Lc.sum(axis=1) > 0).sum())}

    # P1's observed statistic from the files alone
    pos = pd.read_csv(spots_csv, dtype={"cell_id": str}).set_index("cell_id")
    ids = [i for i in pos.index if i in L.index]
    y = np.log2((short.loc[ids].to_numpy() + 0.5) / (long_.loc[ids].to_numpy() + 0.5))
    C, _ = _cn_pcs(L.loc[ids])
    r2 = _r2(y, C)
    ref = json.load(open(a.results))["P1"]["values"] if a.results else {}
    rep["P1_r2_from_files"], rep["P1_r2_reference"] = r2, ref.get("r2_observed")
    print(f"P1 observed R^2 from the files: {r2!r}; reference {ref.get('r2_observed')!r}", flush=True)
    if ref:
        check("P1 observed R^2 from the files equals the reference", abs(r2 - ref["r2_observed"]) < 1e-9, rep)
    return rep


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["takei", "melanoma"])
    ap.add_argument("--store", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rna", help="takei: the RNA counts of P2 (takei2025_fov0_rna_counts.tsv)")
    ap.add_argument("--section", default="human_melanoma_r1")
    ap.add_argument("--results", help="results/bulk.json (the P1 / P2 references to compare with)")
    ap.add_argument("--report", help="JSON report (default <out>/../neutral_build_<what>.json)")
    a = ap.parse_args()
    t0 = time.perf_counter()
    rep = {"what": a.what, "store": a.store, **(takei(a) if a.what == "takei" else melanoma(a))}
    rep["seconds"] = round(time.perf_counter() - t0, 1)
    import chromdata
    rep["chromdata"] = getattr(chromdata, "__version__", "?")
    path = Path(a.report) if a.report else Path(a.out).parent / f"neutral_build_{a.what}.json"
    path.write_text(json.dumps(rep, indent=1, default=str))
    bad = [k for k, v in rep.get("checks", {}).items() if not v]
    print(f"report {path}; {len(rep.get('checks', {}))} checks, failed: {bad}", flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
