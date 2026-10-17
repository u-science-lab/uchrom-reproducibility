"""Captures of the U-Chrom web browser for Fig. 2b: melanoma section r1 of Guo et al. 2026 (spatial Hi-C-RNA) in
four linked views, with a patch of the tissue selected.

The browser colours and plots cells by columns of ``cells``.  The atlas store has the spots' contact-distance
histograms and contacts per chromosome (``cellm["contacts.*"]``), not the derived values, so the script first makes
a capture copy of a local atlas store (``--store``: an APFS clone where the file system makes one, else a copy; the
store itself is not changed) and writes into it

* ``short_long_ratio``: the spot's short / long-range contact ratio, destriped, QC spots of the section
  (``results/spatial_variable_melanoma_spots.parquet``, ``benchmarks/spatial_variable``), clipped to its
  1st-99th percentiles;
* ``chr10_share``: log2 of the spot's share of contact legs on chromosome 10 over the section median (a copy-number
  readout: 2 x cis + trans per chromosome, ``cellm["contacts.cis"]`` / ``["contacts.trans"]``), QC spots;
* ``cell_group``: the authors' cell types in six groups (tumour, T / NK, B / plasma, myeloid, stroma, endothelial).

* ``linked_bin_matrices["scAB"]``: each spot's single-cell A/B compartment score per 1-Mb bin (``--scab``: the
  table of ``benchmarks/spatial_variable`` stage compartments, ``uchrom.strc.comp.single_cell_compartments``
  on the per-spot 1-Mb maps with the hg38 CpG track; 90 MB, kept in the run's cache, not in ``results/``), written as
  an AnnData next to the copy, each 1-Mb value on the store's two 500-kb bins it covers (track ``scAB@1Mb``); the
  scores are each spot's ranks 0-1, shown minus 0.5 (> 0: A-like, < 0: B-like).

Views: the tissue and the RNA UMAP (the authors' scVI) coloured by the ratio, following the selection; two
histograms of every spot with the selected patch overlaid -- the ratio and the chromosome 10 share; a genome view of
a window (``--window``, default chromosome 12, whose copy number the patch keeps): the patch's summed per-spot contact
map, the same as log2 over all spots' (each over its total, 5 x 5 pixels; the browser's ``ratio_to``), the mean A/B
score of the patch and of all spots, and the RefSeq genes (UCSC hg38 refGene, registered as ``hg38_refgene``) with a
few named.

The selected patch is the analysis' copy-number patch (``in_patch``: ``benchmarks/spatial_variable`` stage
patch -- the QC spots in the top fifth of the second principal component of the spots' copy-number profiles,
smoothed over the spot grid; defined by copy number alone).  Then, as ``capture_fig2a.py``: the browser on the copy,
headless Chromium (Paper theme, device scale 2), the views opened, the patch selected with ``select_cells`` (as a
lasso), each view saved with ``capture_view`` and the whole page with the browser's frame.

Writes ``captures/ed_melanoma_<view>.png``, ``captures/ed_melanoma_browser.png`` (the page) and ``ed_melanoma_browser.json`` (where
the views sit on it); the copy goes to ``_cache/`` (not tracked).

    python manuscripts/nmeth-bc/figures/capture_ed_melanoma.py --store PATH/hicrna_human_melanoma.chromdata.zarr
"""

from __future__ import annotations

import argparse
import json
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))  # this repository's root (paper_paths)
from paper_paths import REPRO, UCHROM_REPO  # noqa: E402
sys.path.insert(0, str(REPRO / "benchmarks" / "fig6"))   # this repository
from bench_frontend import LAUNCH, free_port, start_server  # noqa: E402

OUT = HERE / "captures"
CACHE = HERE / "_cache"
SPOTS = HERE / "results" / "spatial_variable_melanoma_spots.parquet"
SECTION = "human_melanoma_r1"
#: genes named in the gene track (those in the window): on 12, CDK4 (amplified / mutated in melanoma), MDM2, KRAS;
#: on 9p, the melanocyte markers MLANA / TYRP1, the tumour suppressors CDKN2A / B and PTPRD, JAK2
LABELLED = {"CDK4", "MDM2", "KRAS", "MLANA", "TYRP1", "CDKN2A", "CDKN2B", "PTPRD", "JAK2"}
GROUPS = {"tumour": ["Melanoma_Invasive", "Melanoma_Melanocytic", "Melanoma_NeuralCrest", "Melanoma_Transitory"],
          "T / NK": ["Activated CD8", "Deletional tolerance CD8", "Dysfunc CD8", "Precursor exhaustion", "Tcm", "Tem",
                     "Terminal exhaustion", "Tfh", "Treg", "NK"],
          "B / plasma": ["B cells memory", "B cells naive", "Plasmablasts"],
          "myeloid": ["MDSC", "Macrophage M1", "Macrophage M2", "Monocytes", "cDC", "pDC"],
          "stroma": ["CAFs MSC iCAF-like", "CAFs myCAF-like", "PVL differentiated", "PVL immature"],
          "endothelial": ["Endothelial"]}
#: view title -> open_view arguments; the views after the selection compare it with every spot
VIEWS = {"tissue": {"type": "spatial", "spatial": "image", "region": SECTION, "image_opacity": 0.25},
         "umap": {"type": "embedding", "embedding": "published_scvi_umap"}}
COMPARE = {"hist_ratio": {"type": "stats", "stats": "hist", "field_x": "short_long_ratio", "cells": "all"},
           "hist_chr10": {"type": "stats", "stats": "hist", "field_x": "chr10_share", "cells": "all"}}


def write_scab(table: Path, cd, out: Path) -> dict:
    """The spots' 1-Mb scA/B scores as an AnnData over the store's (500-kb) bins, and its uns record."""
    import anndata as ad

    t = pd.read_parquet(table)
    t = t.set_index("cell_id") if "cell_id" in t.columns else t
    cols = {}
    for c in t.columns:                                          # "chr1:0-1000000"
        chrom, rng = str(c).split(":")
        cols[(chrom, int(rng.split("-")[0]))] = c
    bins = cd.bins.reset_index(drop=True)
    take = [cols.get((str(ch), int(st) // 1_000_000 * 1_000_000)) for ch, st in zip(bins["chrom"], bins["start"])]
    keep = np.flatnonzero([c is not None for c in take])
    X = t[[take[i] for i in keep]].to_numpy(np.float32) - 0.5      # per-spot ranks 0-1: > 0 A-like, < 0 B-like
    var = pd.DataFrame({"bin_id": keep.astype(np.int64), "resolution": np.full(len(keep), 1_000_000, np.int64)},
                       index=[str(i) for i in keep])
    ad.AnnData(X=X, obs=pd.DataFrame(index=t.index.astype(str)), var=var).write_h5ad(out)
    return {"path": out.name, "format": "h5ad", "label": "scA/B",
            "description": "single-cell A/B compartment score of each spot (U-Chrom single_cell_compartments on the "
                           "per-spot 1-Mb maps: the spot's bins ranked 0-1 by the CpG of their contact partners), "
                           "minus 0.5 (> 0: A-like); each 1-Mb value on the two 500-kb bins it covers",
            "cell_axis": "obs index = cell ids", "var_axis": "var bin_id = row of bins"}


def gene_table(chrom: str, start: int, end: int):
    """RefSeq protein-coding genes (NM_ transcripts, the union per gene) in a window, from UCSC hg38 refGene
    (registered dataset ``hg38_refgene``); the genes of ``LABELLED`` carry their names."""
    import gzip
    import re

    import uchrom.datasets as ds
    from chromdata.intervals import IntervalTable

    f = Path(ds.fetch("hg38_refgene", verbose=False))
    f = next(f.rglob("*.gtf.gz")) if f.is_dir() else f
    spans: dict = {}
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                continue
            col = line.split("\t")
            if col[2] != "transcript":
                continue
            s0, e0 = int(col[3]) - 1, int(col[4])
            if e0 <= start or s0 >= end or 'transcript_id "NM_' not in col[8]:
                continue
            name = re.search(r'gene_name "([^"]+)"', col[8]).group(1)
            a, b = spans.get(name, (s0, e0))
            spans[name] = (min(a, s0), max(b, e0))
    t = pd.DataFrame([(chrom, a, b, n if n in LABELLED else "") for n, (a, b) in spans.items()],
                     columns=["chrom", "start", "end", "label"]).sort_values("start", ignore_index=True)
    print(f"genes: {len(t)} in {chrom}:{start}-{end}, labelled {sorted(set(t.label) - {''})}", flush=True)
    return IntervalTable.from_frame(t, kind="segment")


def capture_store(src: Path, scab: Optional[Path] = None, window: Optional[tuple] = None):
    """The capture copy with the three columns, and the ids of the selected patch."""
    from chromdata import ChromData

    dst = CACHE / f"{src.name.split('.')[0]}.fig2.chromdata.zarr"
    if dst.exists():
        shutil.rmtree(dst)
    CACHE.mkdir(parents=True, exist_ok=True)
    if platform.system() == "Darwin":
        subprocess.run(["cp", "-cR", str(src), str(dst)], check=True)      # APFS clone: no data copied
    else:
        shutil.copytree(src, dst)
    cd = ChromData.read(dst)
    ids = pd.Index(np.asarray(cd._cell_ids_axis(), dtype=str))
    spots = pd.read_parquet(SPOTS)
    sec = spots[spots.section == SECTION]
    v = sec["destriped.q30.short_long_ratio"]                       # the shown section's QC spots only
    lo, hi = np.nanquantile(v, [0.01, 0.99])
    cells = cd.cells.copy()
    cells["short_long_ratio"] = v.clip(lo, hi).reindex(ids).to_numpy()
    chroms = cd.uns["contact_stats"]["contacts"]["chroms"]
    legs = 2 * np.asarray(cd.cellm["contacts.cis"], float) + np.asarray(cd.cellm["contacts.trans"], float)
    legs = pd.DataFrame(legs, index=ids, columns=chroms).drop(columns=["chrY"], errors="ignore")
    share = legs["chr10"] / legs.sum(axis=1).replace(0, np.nan)
    share = share.where(ids.isin(sec.index))
    cells["chr10_share"] = np.log2(share / share.median()).to_numpy()
    if scab is not None:
        cd.uns["linked_bin_matrices"] = {"scAB": write_scab(scab, cd, dst.parent / f"{dst.name.split('.')[0]}_scab.h5ad")}
    if window is not None:
        cd.intervals["genes"] = gene_table(*window)
    group = {ct: g for g, cts in GROUPS.items() for ct in cts}
    cells["cell_group"] = pd.Categorical(cells["cell_type"].map(group), categories=list(GROUPS))
    cd.cells = cells
    cd.write(dst)                               # the store it was read from: keeps its embedded copies
    patch = list(sec.index[sec["in_patch"].astype(bool)])          # the analysis' patch (stage patch)
    print(f"{dst}: short / long clipped to {lo:.2f}-{hi:.2f}; patch of {len(patch)} of {len(sec)} QC spots",
          flush=True)
    return dst, patch


def capture(store: Path, patch: list, size: int, debug: bool, genome: bool = False, window: tuple = ("chr12", 0, 133_275_309)) -> None:
    import requests
    from playwright.sync_api import sync_playwright

    import uchrom_browser as ub

    port = free_port()
    proc = start_server(port)
    base = f"http://127.0.0.1:{port}"
    try:
        requests.post(base + "/api/datasets", json={"path": str(store.resolve())}, timeout=1800).raise_for_status()
        with sync_playwright() as p:
            browser = p.chromium.launch(**LAUNCH)
            page = browser.new_page(viewport={"width": 2 * size, "height": int(2 * size * 0.86)}, device_scale_factor=2)
            page.add_init_script("localStorage.setItem('uchrom.browser.theme', 'paper')")
            page.goto(base + "/")
            ws = ub.connect(base)
            for _ in range(120):                                   # the page connects and opens its default views
                try:
                    if ws.describe().get("views"):
                        break
                except ub.client.BrowserError:
                    pass
                page.wait_for_timeout(500)
            page.wait_for_timeout(2000)
            for side in ("left side bar", "inspector"):
                b = page.locator(f'button[title^="Hide the {side}"]')
                if b.count():
                    b.first.click()
            panes = {title: ws.open_view(title=title, **args)["pane"] for title, args in VIEWS.items()}
            for v in ws.describe()["views"]:
                if v["pane"] not in panes.values():
                    ws.close_view(pane=v["pane"])
            ws.arrange(layout="grid")
            ws.color_by(field="short_long_ratio")
            page.wait_for_timeout(3000)
            ref = ws.select_cells(cells=patch, name="copy-number patch")
            ref = ref.get("group", ref) if isinstance(ref, dict) else ref
            for title, args in COMPARE.items():                    # all spots, the patch overlaid
                panes[title] = ws.open_view(title=title, compare_groups=[ref], **args)["pane"]
            if genome:
                rows = [{"kind": "contacts", "source": f"per_spot.{SECTION}"},
                        {"kind": "contacts", "source": f"per_spot.{SECTION}", "ratio_to": "all", "smooth": 5},
                        {"kind": "track", "track": "scAB@1Mb"},
                        {"kind": "track", "track": "scAB@1Mb", "cells": "all"},
                        {"kind": "intervals", "key": "genes"}]
                panes["genome"] = ws.open_view(title="genome", type="genome", rows=rows)["pane"]
                ws.set_locus(chrom=window[0], start=window[1], end=window[2])
            ws.arrange(layout="grid")
            page.wait_for_timeout(9000)                             # every view re-scopes to the selection
            OUT.mkdir(parents=True, exist_ok=True)
            for title, pane in panes.items():
                r = ws.capture(pane=pane, path=str(OUT / f"ed_melanoma_{title}.png"))
                print(r["path"], r["width"], "x", r["height"], flush=True)
            page.screenshot(path=str(OUT / "ed_melanoma_browser.png"))
            # where each view sits on the page (CSS px; the PNG is at device scale 2): the figure numbers the views
            # at the same place of each one's header
            rects = page.evaluate("""() => Array.from(document.querySelectorAll('section.pane:not(.hidden)')).map((el) => {
                const r = el.getBoundingClientRect(), h = el.querySelector('header.pane-head').getBoundingClientRect()
                return {x: r.left, y: r.top, w: r.width, h: r.height, head: h.height} })""")
            rects.sort(key=lambda r: (round(r["y"]), r["x"]))        # reading order: the numbers 1-5
            vp = page.viewport_size
            (OUT / "ed_melanoma_browser.json").write_text(json.dumps({"scale": 2, "viewport": vp, "views": rects}, indent=1))
            if debug:
                print(ws.describe(), flush=True)
            browser.close()
    finally:
        proc.terminate()
        proc.wait(10)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--store", type=Path, required=True,
                    help="a local hicrna_human_melanoma.chromdata.zarr (the atlas store; it is not changed)")
    ap.add_argument("--size", type=int, default=720, help="width of one view column (CSS px; device scale 2)")
    ap.add_argument("--scab", type=Path, default=None,
                    help="scab_human_melanoma_r1.parquet (stage compartments); adds the genome view")
    ap.add_argument("--window", default="chr12:0-133275309",
                    help="locus of the genome view (default: chromosome 12, whose copy number the patch keeps)")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()
    if not SPOTS.exists():
        raise SystemExit(f"{SPOTS} missing: run U-Chrom's benchmarks/spatial_variable (stage summary) first")
    t = time.time()
    c, r = args.window.split(":")
    window = (c, int(r.split("-")[0]), int(r.split("-")[1]))
    store, patch = capture_store(args.store, args.scab, window if args.scab is not None else None)
    capture(store, patch, args.size, args.debug, genome=args.scab is not None, window=window)
    print(f"{time.time() - t:.0f} s", flush=True)


if __name__ == "__main__":
    main()
