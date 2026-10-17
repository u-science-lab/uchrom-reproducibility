"""Storyboards (wireframes) of the redesigned Brief Communication figures: layout + story, not final art.

Grey boxes = panels to draw; real thumbnails where they exist (atlas figures, browser captures); orange
outline = depends on the feasibility analysis.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image

HERE = Path(__file__).resolve().parent
REPO = Path("/Users/weizexu/orca/workspaces/U-Chrom/fig6")
PAPER = Path("/Users/weizexu/Projects/u-chrom-paper/manuscripts/nmeth-bc/figures")
THUMBS = REPO / "atlas/site/figures"
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial"], "pdf.fonttype": 42})
INK, GREY, LGREY = "#2b2f36", "#8a8f98", "#e3e5e8"
TEAL, ROSE, ORANGE, INDIGO, PURPLE = "#2a9d8f", "#b5507a", "#e07a3f", "#34568b", "#8e7cc3"


class Canvas:
    def __init__(self, w, h):
        self.W, self.H = w, h
        self.fig = plt.figure(figsize=(w * MM, h * MM))
        self.A = self.fig.add_axes([0, 0, 1, 1])
        self.A.set_xlim(0, w); self.A.set_ylim(h, 0); self.A.axis("off")

    def t(self, x, y, s, size=6, color=INK, ha="left", va="top", weight="normal", **kw):
        self.A.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, fontweight=weight, linespacing=1.15, **kw)

    def box(self, x, y, w, h, fc="#f6f7f9", ec=LGREY, lw=0.6, ls="-", r=1.2, z=1):
        self.A.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec,
                                        lw=lw, ls=ls, zorder=z))

    def panel(self, x, y, w, h, letter, title, body=None, pending=False, fc="#f6f7f9"):
        self.box(x, y, w, h, fc=fc, ec=ORANGE if pending else LGREY, lw=0.9 if pending else 0.6,
                 ls=(0, (3, 1.5)) if pending else "-")
        self.t(x + 1.5, y + 1.4, letter, 8, weight="bold")
        self.t(x + 5.5, y + 1.7, title, 6.5, weight="bold")
        if body:
            self.t(x + 5.5, y + 5.6, body, 5.2, GREY)

    def img(self, path, x, y, w=None, h=None, border=True):
        im = Image.open(path).convert("RGB")
        if w is None:
            w = h * im.width / im.height
        if h is None:
            h = w * im.height / im.width
        ax = self.fig.add_axes([x / self.W, 1 - (y + h) / self.H, w / self.W, h / self.H])
        ax.imshow(np.asarray(im)); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(border); sp.set_color(LGREY); sp.set_linewidth(0.4)
        return w, h

    def arrow(self, a, b, color=INDIGO, lw=0.8):
        self.A.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>,head_length=1.6,head_width=0.9", color=color, lw=lw,
                                         mutation_scale=1, shrinkA=0, shrinkB=0, zorder=5))

    def chip(self, x, y, w, s, color=INDIGO, num=None):
        self.box(x, y, w, 7.5, fc="white", ec=color, lw=0.7, r=1.5, z=3)
        if num:
            self.t(x + 1.5, y + 3.75, num, 7, color, va="center", weight="bold", zorder=4)
        self.t(x + (5.5 if num else 1.5), y + 3.75, s, 5.0, INK, va="center", zorder=4)

    def save(self, name):
        self.fig.savefig(HERE / f"{name}.png", dpi=220)
        self.fig.savefig(HERE / f"{name}.pdf")
        return HERE / f"{name}.png"


def fig1():
    c = Canvas(180, 176)
    c.t(0, 0, "Fig. 1  U-Chrom at a glance: many assays → one data model → analysis, a linked browser and an open atlas",
        7, weight="bold")
    # row 1: problem -> model -> store
    y1, h1 = 7, 62
    c.panel(0, y1, 54, h1, "a", "Many assays, three kinds of data",
            "each with its own formats and viewers")
    groups = [("Single-cell Hi-C (+ RNA, ATAC, 3-D)", TEAL,
               "HiRES · Dip-C · GAGE-seq · Droplet / Paired Hi-C\ndscHi-C · scHiCAR · Uni-C",
               "contact map per cell", ".pairs · .scool  |  Juicebox, HiGlass"),
              ("Spatial Hi-C (+ RNA, ATAC)", ROSE, "Hi-C-RNA · spatial Hi-C · spatial ATAC-Hi-C",
               "contact map per spot, tissue position", "pairs + Seurat  |  no viewer for per-spot maps"),
              ("Chromatin tracing (+ IF, RNA)", ORANGE, "DNA seqFISH+ · ORCA · Hi-M · MERFISH",
               "x, y, z per traced locus", "FOF-CT, CSV  |  custom scripts")]
    for k, (name, col, assays, what, fmt) in enumerate(groups):
        yy = y1 + 11 + k * 17
        c.box(2, yy, 50, 15.5, fc="white", ec=col, lw=0.6)
        c.box(3.2, yy + 1.2, 9, 9, fc=col + "33", ec="none")
        c.t(7.7, yy + 5.7, "glyph", 4.5, GREY, ha="center", va="center")
        c.t(13.5, yy + 1.2, name, 5.6, col, weight="bold")
        c.t(13.5, yy + 4.4, assays, 4.8, INK)
        c.t(13.5, yy + 10.2, what, 4.8, INK)
        c.t(13.5, yy + 12.6, fmt, 4.6, GREY)
    c.arrow((54.5, y1 + 34), (58.5, y1 + 34))
    c.panel(59, y1, 79, h1, "b", "ChromData: one index for all three",
            "cells (rows) × genomic bins (columns); spots = cell › trace › bin with x, y, z;\n"
            "contact maps linked per cell, not copied; signals per spot, per bin, per cell")
    c.box(64, y1 + 14, 70, 44, fc="white", ec=LGREY)
    c.t(99, y1 + 36, "the cells × bins schematic\n(current Fig. 1a, refined)", 6, GREY, ha="center", va="center")
    c.arrow((138.5, y1 + 34), (142.5, y1 + 34))
    c.panel(143, y1, 37, h1, "c", "One store",
            "x.chromdata.zarr\nZarr + Parquet, open spec")
    for k, s in enumerate(["read in parts:\none cell, one chromosome", "from disk or\nover HTTP (S3, R2)",
                           "maps linked\nor embedded", ".cdz: one zip file"]):
        c.box(145, y1 + 17 + k * 10.5, 33, 9, fc="white", ec=LGREY)
        c.t(147, y1 + 18.3 + k * 10.5, s, 5.0, INK)
    # row 2: what you do with it
    y2, h2 = 72, 54
    c.panel(0, y2, 54, h2, "d", "Analysis in Python", "one API on every kind of data")
    ops = [("subset / stream", "cell, chromosome, trace;\nlarger than memory"),
           ("aggregate", "pseudo-bulk maps, median\ndistance maps per group"),
           ("convert", "3-D structures ↔ contacts"),
           ("call", "compartments, TADs, loops"),
           ("features", "per locus, per cell"),
           ("link", ".scool, .h5ad, SpatialData")]
    for k, (a, b) in enumerate(ops):
        xx, yy = 2 + (k % 2) * 25.5, y2 + 11 + (k // 2) * 14
        c.box(xx, yy, 24.5, 12.5, fc="white", ec=LGREY)
        c.box(xx + 1, yy + 1, 5, 5, fc=INDIGO + "22", ec="none")
        c.t(xx + 7, yy + 1.2, a, 5.4, INDIGO, weight="bold")
        c.t(xx + 1.2, yy + 6.6, b, 4.6, INK)
    c.panel(57, y2, 81, h2, "e", "One browser: views linked through cells",
            "select cells anywhere; every view follows (people and agents: the same actions)")
    caps = PAPER / "captures"
    hub = (97.5, y2 + 34)
    c.box(hub[0] - 10, hub[1] - 3.5, 20, 7, fc="white", ec=INDIGO, lw=0.9, r=3)
    c.t(hub[0], hub[1], "selected cells", 5.6, INDIGO, ha="center", va="center", weight="bold")
    spokes = [("spatial_tissue", "tissue", (60, y2 + 12)), ("scrna_umap", "embedding", (60, y2 + 33)),
              ("scrna_structure@", "3-D", (121, y2 + 12)), ("scrna_genome", "genome: maps, tracks", (121, y2 + 33))]
    for cap, lab, (x, y) in spokes:
        p = caps / f"{cap.rstrip('@')}.png"
        w, h = c.img(p, x, y, h=17)
        c.t(x + w / 2, y + h + 0.6, lab, 4.8, GREY, ha="center")
        cx = x + w if x < hub[0] else x
        c.A.plot([cx, hub[0] + (-10 if x < hub[0] else 10)], [y + h / 2, hub[1]], color=INDIGO, lw=0.5, zorder=2)
    c.panel(141, y2, 39, h2, "f", "Where you work", None)
    for k, (a, b) in enumerate([("notebook", "Python API, Jupyter"), ("your laptop", "python -m uchrom.browser"),
                                ("the web", "hosted, no install"), ("coding agents", "WebMCP actions, skills")]):
        c.box(143, y2 + 8 + k * 11, 35, 9.5, fc="white", ec=LGREY)
        c.t(145, y2 + 9.4 + k * 11, a, 5.4, INK, weight="bold")
        c.t(145, y2 + 12.8 + k * 11, b, 4.8, GREY)
    # row 3: the atlas
    y3 = 129
    c.panel(0, y3, 180, 46, "g", "An open atlas, ready to browse",
            "22 datasets · 11 studies · 315,516 cells and spots · 2.0 × 10⁸ 3-D loci · 25.6 GB · "
            "contacts per cell checked against the published numbers")
    meta = json.loads((THUMBS / "figures.json").read_text()) if (THUMBS / "figures.json").exists() else {}
    at = json.loads((PAPER / "results/atlas.json").read_text())["datasets"]
    groups = [("Spatial", "spatial Hi-C", ROSE), ("Single-cell", "single-cell Hi-C", TEAL), ("Imaging", "tracing", ORANGE)]
    x = 2.0
    for cat, lab, col in groups:
        ds = [d for d in at if d["category"] == cat]
        c.t(x, y3 + 11, f"{lab} ({len(ds)})", 5.6, col, weight="bold")
        for k, d in enumerate(ds):
            f = THUMBS / f"{d['id']}.space.webp"
            if not f.exists():
                f = THUMBS / f"{d['id']}.map.webp"
            xx, yy = x + (k % 6) * 12.2, y3 + 15 + (k // 6) * 15
            try:
                c.img(f, xx, yy, w=11.5, h=11.5)
            except Exception:
                c.box(xx, yy, 11.5, 11.5)
        x += min(len(ds), 6) * 12.2 + 3
    return c.save("fig1_storyboard")


def fig2():
    c = Canvas(180, 196)
    c.t(0, 0, "Fig. 2  From a tissue section to the genome folding of its cells — and the same interaction on every kind "
              "of data", 7, weight="bold")
    # a: workflow strip
    c.panel(0, 6, 180, 21, "a", "The workflow, end to end (browser ⇄ Python, data from the atlas)", None, fc="white")
    steps = [("1", "open the brain section\nfrom the atlas (HTTP)"), ("2", "spot compartments\non the tissue"),
             ("3", "spot composition\nfrom its RNA"), ("4", "cell-type compartments\nfrom single-cell Hi-C"),
             ("5", "expected vs observed;\nresidual on the tissue"), ("6", "replicate: section 2,\n10-µm spatial Hi-C")]
    for k, (n, s) in enumerate(steps):
        x = 3 + k * 29.5
        c.chip(x, 14, 26.5, s, color=[INDIGO, ROSE, PURPLE, TEAL, INDIGO, GREY][k], num=n)
        if k < len(steps) - 1:
            c.arrow((x + 26.7, 17.75), (x + 29.3, 17.75), lw=0.6)
    # b-d: the data
    caps = PAPER / "captures"
    c.panel(0, 30, 58, 52, "b", "Compartments of each spot", "browser: per-spot A/B on the section")
    c.img(caps / "spatial_tissue.png", 4, 42, h=38)
    c.panel(61, 30, 58, 52, "c", "What each spot is made of", "RNA → cell-type fractions (API)", pending=False)
    c.t(90, 58, "maps of L2/3, L6, Oligo\nfractions on the section", 6, GREY, ha="center", va="center")
    c.panel(122, 30, 58, 52, "d", "Expected from the cell types", "Σ fraction × cell-type A/B\n(GAGE-seq, Paired Hi-C, scHiCAR)")
    c.t(151, 62, "observed vs expected\n(r, null)", 6, GREY, ha="center", va="center")
    # e-g: the result (depends on the feasibility analysis)
    c.panel(0, 85, 58, 52, "e", "What the cell types do not explain", "browser: residual per spot on the section", pending=True)
    c.t(29, 115, "residual map\n(outcome-dependent)", 6, ORANGE, ha="center", va="center")
    c.panel(61, 85, 58, 52, "f", "Replicated", "residual vs cortical depth:\nsection 1, section 2, 10-µm cortex", pending=True)
    c.t(90, 115, "depth profile × 3 datasets\n(outcome-dependent)", 6, ORANGE, ha="center", va="center")
    c.panel(122, 85, 58, 52, "g", "One locus, three views", "maps of deep vs superficial spots;\nthe locus in HiRES 3-D structures", pending=True)
    c.t(151, 115, "example locus\n(outcome-dependent)", 6, ORANGE, ha="center", va="center")
    # h: generality
    c.panel(0, 140, 180, 55, "h", "The same interaction on the other kinds of data",
            "select a cell type; every view follows (browser captures)")
    rows = [("single-cell Hi-C + RNA (HiRES brain): Oli", ["scrna_umap", "scrna_structure", "scrna_genome"], TEAL),
            ("chromatin tracing + IF (Takei cerebellum): Bergmann glia", ["tracing_umap", "tracing_nuclei", "tracing_genome"], ORANGE)]
    for k, (lab, ims, col) in enumerate(rows):
        x0 = 4 + k * 89
        c.t(x0, 151, lab, 5.4, col, weight="bold")
        for j, n in enumerate(ims):
            c.img(caps / f"{n}.png", x0 + j * 28.5, 156, w=27.5)
    c.t(0, 193, "orange dashed: content depends on the feasibility analysis (residual present → position-dependent folding; "
                "absent → cell types explain folding in situ)", 5, ORANGE)
    return c.save("fig2_storyboard")


if __name__ == "__main__":
    print(fig1())
    print(fig2())
