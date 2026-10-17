"""Captures of the U-Chrom web browser for Fig. 2d (Extended Data Fig. 3d until 2026-10-12): chromosome 2 in one typical
ExN and one typical Oligo cell of the imaging and of two NucDynamics reconstructions, drawn by the browser's 3D view.

The structures are those of ``make_fig2.structure_views`` (on the imaged loci, each reconstruction superposed
onto the imaged cell of its class, all in that cell's principal axes, um).  The script writes them into a small store
(``_cache/fig2_structures.chromdata.zarr``: the imaged chr2 loci as bins, one cell and one trace per structure, the
loci a structure lacks left out, so the tube runs through the ones it has), starts the browser on it, opens a 3D view
per structure (``open_view(type="3d", cells=<id>)``, Paper theme, device scale 2; in the inspector one tube radius for
all, ``--radius``, and high-quality tubes) and saves each with ``capture_view``.  The browser colours a single chromosome by position along it (violet, start; red, end) and frames
a view on the structure's bounding box looking along -z (its camera: perspective, 45 degree vertical field of view,
at 1.1 x the box's half diagonal / sin 22.5 degree from the box centre), so the scale at the box centre is known:
``captures/fig2_structures.json`` records, per capture, the centre (x, y, um) and the pixels per um there, which the
figure uses to put every capture on one scale.

    python manuscripts/nmeth-bc/figures/capture_fig2_structures.py [--size 520] [--radius 0.12] [--debug]
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))  # this repository's root (paper_paths)
from paper_paths import REPRO, UCHROM_REPO  # noqa: E402
sys.path.insert(0, str(REPRO / "benchmarks" / "fig6"))   # this repository
from bench_frontend import LAUNCH, free_port, start_server  # noqa: E402
sys.path.insert(0, str(HERE))
from make_fig2 import PANEL, STRUCT_SOURCES, structure_views  # noqa: E402

OUT = HERE / "captures"
CACHE = HERE / "_cache"
STORE = CACHE / "fig2_structures.chromdata.zarr"
FOV = 45.0                                     # the 3D view's vertical field of view (Viewport.tsx)


def cell_id(cls: str, src: str) -> str:
    return f"{cls} {src.split(' ')[0]}"


def write_store(views: dict) -> Path:
    from chromdata import ChromData

    loci = PANEL[PANEL["chrom"].astype(str) == "chr2"].sort_values("start").reset_index(drop=True)
    bins = pd.DataFrame({"chrom": "chr2", "start": loci["start"].to_numpy(np.int64),
                         "end": loci["end"].to_numpy(np.int64)})
    spots, coords = [], []
    for (cls, src), X in views.items():
        ok = np.flatnonzero(np.isfinite(X).all(1))
        cid = cell_id(cls, src)
        spots.append(pd.DataFrame({"bin_id": ok, "trace_id": cid, "cell_id": cid}))
        coords.append(X[ok])
    cells = pd.DataFrame({"cls": [c for c, _ in views], "source": [s for _, s in views]},
                         index=pd.Index([cell_id(c, s) for c, s in views], name="cell_id"))
    cd = ChromData(np.vstack(coords), pd.concat(spots, ignore_index=True), bins=bins, cells=cells,
                   uns={"xyz_unit": "um", "genome_assembly": "mm10",
                        "description": "Fig. 2d example structures (make_fig2.structure_views)"})
    if STORE.exists():
        shutil.rmtree(STORE)
    CACHE.mkdir(exist_ok=True)
    cd.write(str(STORE))
    return STORE


def set_tubes(page, radius: float) -> None:
    """The inspector's tube radius (one for every structure: the browser's auto radius follows each structure's bond
    length) and its high-quality tubes, on the active view's layer."""
    page.locator('button[title^="Show the inspector"]').first.click()
    page.wait_for_timeout(600)
    page.locator('input[type="number"][placeholder="auto"]').first.fill(f"{radius}")
    hq = page.locator('label.check:has-text("High quality") input[type="checkbox"]').first
    if not hq.is_checked():
        hq.check()
    page.wait_for_timeout(400)
    page.locator('button[title^="Hide the inspector"]').first.click()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--size", type=int, default=520, help="one view, CSS px (square; device scale 2)")
    ap.add_argument("--radius", type=float, default=0.12, help="tube radius, um (one for all structures)")
    ap.add_argument("--debug", action="store_true", help="also save the whole page")
    a = ap.parse_args()
    import requests
    from playwright.sync_api import sync_playwright

    import uchrom_browser as ub

    views = structure_views()
    store = write_store(views)
    port = free_port()
    proc = start_server(port)
    base = f"http://127.0.0.1:{port}"
    meta = {}
    try:
        requests.post(base + "/api/datasets", json={"path": str(store.resolve())}, timeout=600).raise_for_status()
        with sync_playwright() as p:
            browser = p.chromium.launch(**LAUNCH)
            page = browser.new_page(viewport={"width": a.size, "height": a.size + 120}, device_scale_factor=2)
            page.add_init_script("localStorage.setItem('uchrom.browser.theme', 'paper')")
            page.goto(base + "/")
            ws = ub.connect(base)
            for _ in range(120):
                try:
                    if ws.describe().get("views"):
                        break
                except ub.client.BrowserError:
                    pass
                page.wait_for_timeout(500)
            page.wait_for_timeout(1500)
            for side in ("left side bar", "inspector"):
                b = page.locator(f'button[title^="Hide the {side}"]')
                if b.count():
                    b.first.click()
            default = [v["pane"] for v in ws.describe()["views"]]
            OUT.mkdir(parents=True, exist_ok=True)
            for (cls, src), X in views.items():
                cid = cell_id(cls, src)
                pane = ws.open_view(type="3d", cells=cid, title=cid)["pane"]
                for pn in default:
                    ws.close_view(pane=pn)
                default = []
                ws.arrange(layout="tabs")
                page.wait_for_timeout(2500)
                set_tubes(page, a.radius)
                page.wait_for_timeout(2500)
                name = f"fig2_{cls}_{src.split(' ')[0]}".replace("-", "").lower()
                r = ws.capture(pane=pane, path=str(OUT / f"{name}.png"))
                P = X[np.isfinite(X).all(1)]
                lo, hi = P.min(0), P.max(0)
                half = np.linalg.norm(hi - lo) / 2
                dist = half / np.sin(np.radians(FOV / 2)) * 1.1
                px_per_um = (r["height"] / 2) / (dist * np.tan(np.radians(FOV / 2)))
                c = (lo + hi) / 2
                meta[f"{cls}|{src}"] = {"png": f"{name}.png", "width": r["width"], "height": r["height"],
                                         "center_um": [float(c[0]), float(c[1])], "px_per_um": float(px_per_um),
                                         "n_loci": int(len(P))}
                print(r["path"], r["width"], "x", r["height"], f"{px_per_um:.1f} px/um", flush=True)
                if a.debug:
                    page.screenshot(path=str(OUT / f"{name}_page.png"))
                ws.close_view(pane=pane)
            browser.close()
    finally:
        proc.terminate()
        proc.wait(10)
    (OUT / "fig2_structures.json").write_text(json.dumps(meta, indent=1) + "\n")


if __name__ == "__main__":
    main()
