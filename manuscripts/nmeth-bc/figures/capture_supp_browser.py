"""Supplementary Figs 2-4: the whole U-Chrom browser page on each kind of data, one figure per dataset.

For each dataset the script starts ``python -m uchrom_browser`` on a free port (the U-Chrom checkout of
``paper_paths``), opens the dataset's store of the public atlas over HTTP (the store a reader opens with the link in
the figure: ``HOSTED/?open=<store URL>``), opens the page in headless Chromium (Playwright, Paper theme) at a page size that keeps the
interface's text legible across a printed page (``--width`` x ``--height`` CSS px, device scale ``--scale``), and
drives it with the browser's own actions (``uchrom_browser.connect()``): it opens the views that suit the dataset,
closes the default ones, arranges them as a grid (three views: the first two side by side, the genome view across the
page below; four: two by two), colours the cells by cell type and focuses one cell type, which every view follows.  The left side bar stays (dataset, fields,
groups); the inspector is hidden.  The page is saved as it is shown:

* ``spatial``  spatial Hi-C-RNA, mouse brain section r1 (Guo et al. 2026), three views (no 3D coordinates, so no
               distance statistics): the tissue with its image, the RNA UMAP, a genome view of chr2 with the selected
               spots' contacts and A/B score against all spots';
* ``scrna``    single-cell Hi-C + RNA, adult mouse brain (HiRES, Liu et al. 2023), four views: the RNA UMAP, the 3D
               structure of one oligodendrocyte (a view that keeps its own cell), the chr2 contacts of the selection
               and of all cells, and 3D distance against genomic separation, the selection over all cells;
* ``tracing``  DNA seqFISH+, mouse cerebellum (Takei et al. 2025), four views: the UMAP, the nuclei in 3D,
               chr15:20-25 Mb (the median 3D distances of the selection and of all cells, the selection's
               immunofluorescence) and 3D distance against genomic separation, the selection over all cells.

Writes ``captures/supp_browser_<dataset>.png`` and, framed as a browser window (title bar with the window buttons and
the padlock of a secure address and the address of the dataset in the hosted browser; rounded border, shadow),
``captures/supp_browser_<dataset>_framed.png``, and ``captures/supp_browser_<dataset>_framed.tex``, which the supplement
inputs: the image with a link to that address over its address field only.

    python manuscripts/nmeth-bc/figures/capture_supp_browser.py [--only spatial scrna tracing] [--frame-only]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))  # this repository's root (paper_paths)
from paper_paths import REPRO, UCHROM_REPO  # noqa: E402
sys.path.insert(0, str(REPRO / "benchmarks" / "fig6"))   # this repository
from bench_frontend import LAUNCH, free_port, start_server  # noqa: E402

OUT = HERE / "captures"
ATLAS = "https://uchrom-atlas-r2.u-science.org"       # the public atlas stores (catalog.json)
HOSTED = "https://uchrom-browser.u-science.org"      # the hosted browser; ?open=<store URL> opens a store
OLI_CELL = "ValaB8w3043"           # the best-determined Oli structure (rmsd_20k), as in Fig. 2a's capture

#: dataset -> store, the field to colour by, the group to focus, the locus, the views (title, open_view arguments)
SCENES = {
    "spatial": {
        "store": f"{ATLAS}/hicrna_mouse_brain.chromdata.zarr",
        "color_by": "cell_type", "group": "field:cell_type=L6 CT", "locus": ("chr2", 0, 182_000_000),
        "views": [("tissue", {"type": "spatial", "spatial": "image", "region": "mouse_brain_r1"}),
                  ("RNA UMAP", {"type": "embedding", "embedding": "published_scvi_umap"}),
                  ("genome", {"type": "genome", "rows": [{"kind": "contacts", "source": "per_spot.mouse_brain_r1"},
                                                         {"kind": "track", "track": "ssAB@500kb"},
                                                         {"kind": "contacts", "source": "per_spot.mouse_brain_r1",
                                                          "cells": "all"},
                                                         {"kind": "track", "track": "ssAB@500kb", "cells": "all"}]})],
    },
    "scrna": {
        "store": f"{ATLAS}/hires_brain.chromdata.zarr",
        "color_by": "cell_type", "group": "field:cell_type=Oli", "locus": ("chr2", 0, 182_000_000),
        "views": [("RNA UMAP", {"type": "embedding", "embedding": "rna_umap"}),
                  ("one cell in 3D", {"type": "3d", "cells": OLI_CELL}),
                  ("genome", {"type": "genome", "rows": [{"kind": "contacts", "source": "per_cell"},
                                                         {"kind": "contacts", "source": "per_cell",
                                                          "cells": "all"}]}),
                  ("distance by separation", {"type": "stats", "stats": "scaling", "cells": "all",
                                              "compare_groups": ["field:cell_type=Oli"]})],
    },
    "tracing": {
        "store": f"{ATLAS}/takei2025_cerebellum.chromdata.zarr",
        "color_by": "cell_type", "group": "field:cell_type=Bergmann", "locus": ("chr15", 20_000_000, 25_000_000),
        "views": [("UMAP", {"type": "embedding", "embedding": "umap"}),
                  ("nuclei", {"type": "3d"}),
                  ("genome", {"type": "genome", "rows": [{"kind": "distance", "stat": "median"},
                                                         {"kind": "distance", "stat": "median", "cells": "all"},
                                                         {"kind": "track", "track": "H3K27ac"},
                                                         {"kind": "track", "track": "H3K9me3"},
                                                         {"kind": "track", "track": "LaminB1"}]}),
                  ("distance by separation", {"type": "stats", "stats": "scaling", "cells": "all",
                                              "compare_groups": ["field:cell_type=Bergmann"]})],
    },
}


def scene(name: str, spec: dict, width: int, height: int, scale: float, settle: float = 30.0) -> Path:
    import requests
    from playwright.sync_api import sync_playwright

    import uchrom_browser as ub

    port = free_port()
    proc = start_server(port)
    base = f"http://127.0.0.1:{port}"
    try:
        requests.post(base + "/api/datasets", json={"path": spec["store"]}, timeout=1800).raise_for_status()
        with sync_playwright() as p:
            browser = p.chromium.launch(**LAUNCH)
            page = browser.new_page(viewport={"width": width, "height": height}, device_scale_factor=scale)
            page.add_init_script("localStorage.setItem('uchrom.browser.theme', 'paper')")
            page.goto(base + "/")
            ws = ub.connect(base, timeout=300)                     # views read the store over HTTP
            for _ in range(120):                                   # the page connects and opens its default views
                try:
                    if ws.describe().get("views"):
                        break
                except ub.client.BrowserError:
                    pass
                page.wait_for_timeout(500)
            page.wait_for_timeout(2000)
            default = [v["pane"] for v in ws.describe()["views"]]
            b = page.locator('button[title^="Hide the inspector"]')
            if b.count():
                b.first.click()
            # the statistics last, once the other views have read what they need: a remote store answered a
            # statistic requested together with the other views' reads with a server error once
            first = [(t, v) for t, v in spec["views"] if v["type"] != "stats"]
            later = [(t, v) for t, v in spec["views"] if v["type"] == "stats"]
            panes = [ws.open_view(title=title, **args)["pane"] for title, args in first]
            for pn in default:
                if pn not in panes:
                    ws.close_view(pane=pn)
            ws.arrange(layout="grid")
            ws.set_locus(chrom=spec["locus"][0], start=spec["locus"][1], end=spec["locus"][2])
            ws.color_by(field=spec["color_by"])
            page.wait_for_timeout(15000)                           # the store is read over HTTP
            ws.focus(group=spec["group"])
            try:
                page.wait_for_load_state("networkidle", timeout=180000)
            except Exception:                                      # a view keeps polling: wait a fixed time instead
                pass
            page.wait_for_timeout(30000)                           # every view re-scopes to the group
            if later:
                panes += [ws.open_view(title=title, **args)["pane"] for title, args in later]
                ws.arrange(layout="grid")
                page.wait_for_timeout(20000)
            page.wait_for_timeout(settle * 1000)                   # a large remote store: the views keep reading
            try:
                page.wait_for_load_state("networkidle", timeout=300000)
            except Exception:
                pass
            # no view still loading (a genome row's or a map pane's "loading…", the activity bar): a fixed wait once
            # caught the group's map before it came, showing all spots' map under the group's name (2026-10-14)
            page.wait_for_function("() => !document.querySelector('.gv-row-loading, .pane-loading, .activity-bar')",
                                   timeout=1800000, polling=2000)
            page.wait_for_timeout(5000)                            # the canvases draw what arrived
            OUT.mkdir(parents=True, exist_ok=True)
            path = OUT / f"supp_browser_{name}.png"
            page.screenshot(path=str(path))
            browser.close()
        return path
    finally:
        import subprocess
        proc.terminate()
        try:
            proc.wait(30)
        except subprocess.TimeoutExpired:                  # still reading over HTTP
            proc.kill()
            proc.wait(10)


def link(name: str) -> str:
    """The dataset in the hosted browser."""
    return f"{HOSTED}/?open={SCENES[name]['store']}"


def frame(path: Path, address: str) -> Path:
    """The page as a browser window: a title bar with the three window buttons and an address field above it, a
    rounded border and a soft shadow; sizes in mm of the printed page (the capture spans 160 mm)."""
    from matplotlib import font_manager
    from PIL import Image, ImageDraw, ImageFilter, ImageFont

    page = Image.open(path).convert("RGB")
    w, h = page.size
    px = w / 160.0                                          # pixels per mm
    bar, r, pad, line = round(6.5 * px), round(1.6 * px), round(2.2 * px), max(2, round(0.18 * px))
    W, H = w + 2 * pad, h + bar + 2 * pad
    out = Image.new("RGB", (W, H), "white")
    shadow = Image.new("L", (W, H), 0)
    ImageDraw.Draw(shadow).rounded_rectangle((pad, pad + round(0.5 * px), pad + w, pad + bar + h + round(0.5 * px)),
                                             radius=r, fill=70)
    out.paste(Image.new("RGB", (W, H), "#9aa3b2"), mask=shadow.filter(ImageFilter.GaussianBlur(round(1.0 * px))))
    win = (pad, pad, pad + w - 1, pad + bar + h - 1)
    body = Image.new("RGB", (w, bar + h), "#e9eef6")       # the window: title bar and page, its corners rounded
    body.paste(page, (0, bar))
    ImageDraw.Draw(body).rectangle((0, bar - line, w - 1, bar - 1), fill="#d3d9e3")
    mask = Image.new("L", (w, bar + h), 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, w - 1, bar + h - 1), radius=r, fill=255)
    out.paste(body, (pad, pad), mask)
    d = ImageDraw.Draw(out)
    cy = pad + bar / 2
    for k, c in enumerate(("#ec6a5e", "#f4bf4f", "#61c554")):
        x = pad + (3.4 + 3.0 * k) * px
        d.ellipse((x - 0.95 * px, cy - 0.95 * px, x + 0.95 * px, cy + 0.95 * px), fill=c)
    ax0, ax1, ah = pad + 15 * px, pad + w - 15 * px, 2.1 * px
    d.rounded_rectangle((ax0, cy - ah, ax1, cy + ah), radius=round(ah), fill="white")
    # the padlock a browser shows before a secure (https) address, then the address
    lx, lw, lh, grey = ax0 + 1.8 * px, 1.7 * px, 1.35 * px, "#5b6270"
    body_top = cy - 0.15 * px
    d.rounded_rectangle((lx, body_top, lx + lw, body_top + lh), radius=round(0.25 * px), fill=grey)
    sw = max(2, round(0.28 * px))                            # the shackle: an arch above the body
    d.arc((lx + 0.3 * px, body_top - 1.25 * px, lx + lw - 0.3 * px, body_top + 0.35 * px), 180, 360, fill=grey, width=sw)
    d.line((lx + 0.3 * px + sw / 2, body_top - 0.45 * px, lx + 0.3 * px + sw / 2, body_top + 0.1 * px), fill=grey, width=sw)
    d.line((lx + lw - 0.3 * px - sw / 2, body_top - 0.45 * px, lx + lw - 0.3 * px - sw / 2, body_top + 0.1 * px),
           fill=grey, width=sw)
    font = ImageFont.truetype(font_manager.findfont(font_manager.FontProperties(family=["Arial", "DejaVu Sans"])),
                              round(2.3 * px))
    d.text((lx + lw + 1.3 * px, cy), address.removeprefix("https://"), fill=grey, font=font, anchor="lm")
    d.rounded_rectangle(win, radius=r, outline="#9fb0cc", width=line)
    dest = path.with_name(path.stem + "_framed.png")
    out.save(dest, optimize=True)
    # the LaTeX of the figure: the image, and a link over the address field only (picture units: the image width)
    rel = dest.relative_to(HERE).as_posix()
    box = dict(x=ax0 / W, y=(H - (cy + ah)) / W, w=(ax1 - ax0) / W, h=2 * ah / W)
    dest.with_suffix(".tex").write_text(
        "% written by capture_supp_browser.py: the framed capture, linked to the dataset in the hosted browser from its\n"
        "% address field only\n"
        "\\setlength{\\unitlength}{\\textwidth}%\n"
        f"\\begin{{picture}}(1,{H / W:.5f})%\n"
        f"\\put(0,0){{\\includegraphics[width=\\textwidth]{{{rel}}}}}%\n"
        # a link annotation of the field's size (xdvipdfmx: \\href would size it by what it draws, here nothing)
        f"\\put({box['x']:.5f},{box['y']:.5f}){{\\special{{pdf:ann width \\the\\dimexpr{box['w']:.5f}\\textwidth\\relax"
        f"\\space height \\the\\dimexpr{box['h']:.5f}\\textwidth\\relax\\space depth 0pt "
        f"<< /Type /Annot /Subtype /Link /Border [0 0 0] /A << /S /URI /URI ({address}) >> >>}}}}%\n"
        "\\end{picture}%\n")
    return dest


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", choices=sorted(SCENES))
    ap.add_argument("--width", type=int, default=1200, help="page width, CSS px")
    ap.add_argument("--height", type=int, default=1300, help="page height, CSS px")
    ap.add_argument("--scale", type=float, default=2.5, help="device scale factor")
    ap.add_argument("--frame-only", action="store_true", help="only frame the existing captures")
    ap.add_argument("--settle", type=float, default=30.0, help="seconds to wait before the screenshot")
    a = ap.parse_args()
    for name, spec in SCENES.items():
        if not a.only or name in a.only:
            t = time.time()
            path = (OUT / f"supp_browser_{name}.png" if a.frame_only
                    else scene(name, spec, a.width, a.height, a.scale, a.settle))
            print(frame(path, link(name)), link(name), f"{time.time() - t:.0f} s", flush=True)


if __name__ == "__main__":
    main()
