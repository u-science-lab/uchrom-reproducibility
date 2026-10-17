"""Captures of the U-Chrom web browser for Fig. 2a: one interaction on the three kinds of data.

For each dataset the script starts ``python -m uchrom.browser --web`` (the U-Chrom checkout of ``paper_paths``),
opens the page in headless Chromium (Playwright, Paper theme, device scale 2) and drives it with the browser's own
actions (``uchrom.browser.connect()``, the actions of ``web/ui_tools.json``): it opens the views below, colours
the cells by cell type, focuses one cell type -- the selection every view then follows -- and saves each view
with ``capture_view`` (the view alone, without the interface):

* ``spatial``  spatial Hi-C-RNA, mouse brain section r1 (Guo et al. 2026): the tissue view; the authors' RNA UMAP
               (scVI); a genome view with the focused spots' pseudo-bulk contact map (per-spot maps summed, chr2)
               and their single-spot A/B scores against all spots;
* ``scrna``    single-cell Hi-C + RNA, adult mouse brain (HiRES, Liu et al. 2023): the RNA UMAP; the 3D structure
               of one focused cell (a view that keeps its own cell); a genome view with the focused cells'
               pseudo-bulk contact map and all cells' (chr2);
* ``tracing``  DNA seqFISH+, mouse cerebellum (Takei et al. 2025): the UMAP; the focused cells' nuclei in 3D;
               the median distance maps of the focused cells and of all cells, and the focused cells'
               immunofluorescence, chr15:20-25 Mb (25-kb loci; a large store's distance map uses a random subset of
               the traces, about 1,500 here).

Writes ``captures/<dataset>_<view>.png`` next to this script (and ``*_page.png``, the whole page, with --debug).

    python manuscripts/nmeth-bc/figures/capture_fig2a.py [--only spatial scrna tracing] [--debug]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))  # this repository's root (paper_paths)
from paper_paths import EXAMPLE_DATA, REPRO, UCHROM_REPO  # noqa: E402
sys.path.insert(0, str(REPRO / "benchmarks" / "fig6"))   # this repository
from bench_frontend import LAUNCH, free_port, start_server  # noqa: E402

OUT = HERE / "captures"

#: dataset -> store, the group to focus, the views (pane title, open_view arguments), the locus
SCENES = {
    "spatial": {
        "store": EXAMPLE_DATA / "spatial_hicrna" / "hicrna_mouse_brain.chromdata.zarr",
        "color_by": "cell_type", "group": "field:cell_type=L6 CT",
        "locus": ("chr2", 0, 182_000_000),
        "views": [("tissue", {"type": "spatial", "spatial": "image", "region": "mouse_brain_r1"}),
                  ("umap", {"type": "embedding", "embedding": "published_scvi_umap"}),
                  ("genome", {"type": "genome", "rows": [{"kind": "contacts", "source": "per_spot.mouse_brain_r1"},
                                                           {"kind": "track", "track": "ssAB@500kb"},
                                                           {"kind": "contacts", "source": "per_spot.mouse_brain_r1", "cells": "all"},
                                                           {"kind": "track", "track": "ssAB@500kb", "cells": "all"}]})],
    },
    "scrna": {
        "store": EXAMPLE_DATA / "hires2023" / "hires_brain.chromdata.zarr",
        "color_by": "cell_type", "group": "field:cell_type=Oli",
        "cell": ("structure", "ValaB8w3043"),     # then focus one cell: the best-determined Oli structure (rmsd_20k);
                                                  # not in the sampled overview, so the browser fetches it on focus
                                                  # and capture that view
        "locus": ("chr2", 0, 182_000_000),
        "views": [("umap", {"type": "embedding", "embedding": "rna_umap"}),
                  ("structure", {"type": "3d"}),
                  ("genome", {"type": "genome", "rows": [{"kind": "contacts", "source": "per_cell"},
                                                           {"kind": "contacts", "source": "per_cell", "cells": "all"}]})],
    },
    "tracing": {
        "store": EXAMPLE_DATA / "takei2025_cerebellum" / "takei2025_cerebellum_rep1.chromdata.zarr",
        "color_by": "cell_type", "group": "field:cell_type=Bergmann",
        "locus": ("chr15", 20_000_000, 25_000_000),    # 25-kb loci: the large-store trace budget covers 5 Mb well
        "views": [("umap", {"type": "embedding", "embedding": "umap"}),
                  ("nuclei", {"type": "3d"}),
                  ("genome", {"type": "genome", "rows": [{"kind": "distance", "stat": "median"},
                                                           {"kind": "distance", "stat": "median", "cells": "all"},
                                                           {"kind": "track", "track": "H3K27ac"},
                                                           {"kind": "track", "track": "H3K9me3"},
                                                           {"kind": "track", "track": "LaminB1"}]})],
    },
}


def scene(name: str, spec: dict, size: int, debug: bool) -> None:
    import requests
    from playwright.sync_api import sync_playwright

    import uchrom.browser as ub

    port = free_port()
    proc = start_server(port)
    base = f"http://127.0.0.1:{port}"
    try:
        requests.post(base + "/api/datasets", json={"path": str(Path(spec["store"]).resolve())}, timeout=1800).raise_for_status()
        with sync_playwright() as p:
            browser = p.chromium.launch(**LAUNCH)
            n = len(spec["views"])
            page = browser.new_page(viewport={"width": size * n, "height": int(size * 1.15)}, device_scale_factor=2)
            page.add_init_script("localStorage.setItem('uchrom.browser.theme', 'paper')")
            page.goto(base + "/")
            ws = ub.connect(base)
            for _ in range(120):                                   # the page connects and opens its default views
                try:
                    if ws.describe().get("views"):
                        break
                except ub.client.BrowserError:                    # not connected yet
                    pass
                page.wait_for_timeout(500)
            page.wait_for_timeout(2000)
            default = [v["pane"] for v in ws.describe()["views"]]
            for side in ("left side bar", "inspector"):            # the views get the whole width
                b = page.locator(f'button[title^="Hide the {side}"]')
                if b.count():
                    b.first.click()
            panes = {}
            for title, args in spec["views"]:
                panes[title] = ws.open_view(title=title, **args)["pane"]
            for v in ws.describe()["views"]:                      # then close the page's default views
                if v["pane"] not in panes.values():
                    ws.close_view(pane=v["pane"])
            ws.arrange(layout="columns")
            ws.set_locus(chrom=spec["locus"][0], start=spec["locus"][1], end=spec["locus"][2])
            ws.color_by(field=spec["color_by"])
            page.wait_for_timeout(4000)
            ws.focus(group=spec["group"])
            page.wait_for_timeout(9000)                             # every view re-scopes to the group
            OUT.mkdir(parents=True, exist_ok=True)
            later = spec.get("cell", (None, None))[0]
            for title, pane in panes.items():
                if title == later:
                    continue
                r = ws.capture(pane=pane, path=str(OUT / f"{name}_{title}.png"))
                print(r["path"], r["width"], "x", r["height"], flush=True)
            if later:                                              # one cell of the group, in the view that follows it
                ws.focus(cell=spec["cell"][1])
                page.wait_for_timeout(3000)
                # the other cells off (the inspector's "others" slider): every structure is centred on the origin
                page.locator('button[title^="Show the inspector"]').first.click()
                page.wait_for_timeout(800)
                page.locator('div.emb-row[title^="Opacity of the other cells"] input[type="range"]').first.evaluate(
                    """(el) => { const set = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set
                       set.call(el, '0'); el.dispatchEvent(new Event('input', {bubbles: true}))
                       el.dispatchEvent(new Event('change', {bubbles: true})) }""")
                page.locator('button[title^="Hide the inspector"]').first.click()
                page.wait_for_timeout(8000)
                r = ws.capture(pane=panes[later], path=str(OUT / f"{name}_{later}.png"))
                print(r["path"], r["width"], "x", r["height"], flush=True)
            if debug:
                page.screenshot(path=str(OUT / f"{name}_page.png"))
            browser.close()
    finally:
        proc.terminate()
        proc.wait(10)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", choices=sorted(SCENES))
    ap.add_argument("--size", type=int, default=640, help="width of one view (CSS px; device scale 2)")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()
    for name, spec in SCENES.items():
        if not args.only or name in args.only:
            t = time.time()
            scene(name, spec, args.size, args.debug)
            print(f"{name}: {time.time() - t:.0f} s", flush=True)


if __name__ == "__main__":
    main()
