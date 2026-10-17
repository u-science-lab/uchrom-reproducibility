"""Fig. 1 -- one data model and a linked browser for three kinds of 3-D genome data.

  a  one cell, measured three ways: spatial Hi-C, single-cell Hi-C and chromatin tracing (illustrations
     ``art/a_*.png``) measure one cell (the nucleus of the full Fig. 1a illustration,
     ``../../full/figures/fig1/art3/a_flat_1.png``, turned so that its territories match the plates); the cell
     opens into ChromData's layers on shared genomic bins -- geometry, contacts, molecular features and the bridges
     between them -- drawn as vectors with the full Fig. 1b code (``make_fig1_final.panel_b``, real-data shapes)
  b  one browser, linked through cells (a schematic), numbered in the order the selection flows: (1) an RNA
     embedding with a lasso, (2) a piece of cerebellar cortex drawn with anatomically fixed features only (layers;
     granule, Purkinje, Bergmann, molecular layer interneurons), the lassoed cells marked, (3) nuclei and chromatin
     traces in 3-D, (4) genome tracks and a 3-D distance map, lasso vs all; (1), (3) and (4) redrawn from
     ``../../full/figures/fig6/_cache/takei_fov0.npz``
  c  an open ecosystem: interfaces, methods, the ChromData core and the formats it reads, writes or links (dashed: can
     be linked rather than copied, explained in the legend; FOF-CT on both sides, as ChromData reads and writes it),
     with the communities' logos (drawn with the full Fig. 1c helpers,
     ``../../full/figures/fig1/make_fig1_final.py``)

The panels carry their letters only, as in the journal's published figures: the titles and subtitles of version 7
are in the legend.

The illustrations in ``art/`` and the nucleus are generated (``art/PROMPTS.md``); nothing in them is data.
Benchmarks are in the Supplementary Information.

    python figures/fig1/make_fig1.py      # -> fig1.pdf / fig1.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
from PIL import Image

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "shared"
sys.path.insert(0, str(SHARED / "ecosystem"))
import make_fig1_final as f1  # noqa: E402  (ecosystem helpers and logos)

ART = SHARED / "art"
MM = 1 / 25.4
W, H = 180.0, 150.2                                    # panel c ends with three format rows, 8 mm apart
f1.W, f1.H = W, H                                    # the full figure's helpers place in this canvas

mpl.rcParams.update({"pdf.fonttype": 42, "font.family": "sans-serif", "font.sans-serif": ["Arial"]})

ORANGE, TEAL, PURPLE, INDIGO = "#e07a3f", "#2a9d8f", "#8e7cc3", "#34568b"
GREY, LGREY, INK, LINE, FRAME = "#8a8f98", "#c9ccd1", "#2b2f36", "#3a3f47", "#d9dce1"
ROSE = "#b5507a"

fig = plt.figure(figsize=(W * MM, H * MM))
A = fig.add_axes([0, 0, 1, 1])
A.set_xlim(0, W); A.set_ylim(H, 0); A.axis("off")
TOP = fig.add_axes([0, 0, 1, 1], zorder=20)          # marks drawn over the images (same mm frame)
TOP.set_xlim(0, W); TOP.set_ylim(H, 0); TOP.axis("off"); TOP.patch.set_alpha(0)


# --------------------------------------------------------------------------- canvas helpers (mm, y down)

def T(x, y, s, size=5.5, color=INK, ha="left", va="center", weight="normal", ax=None, **kw):
    return (ax or A).text(x, y, s, fontsize=size, color=color, ha=ha, va=va, fontweight=weight, linespacing=1.05,
                          **kw)


def rect(x, y, w, h, fc, ec=LINE, lw=0.4, ls="-", z=2, r=0.0, ax=None):
    if r:
        p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, ls=ls,
                           zorder=z)
    else:
        p = Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, ls=ls, zorder=z)
    (ax or A).add_patch(p)


def line(xs, ys, color=LINE, lw=0.3, ls="-", z=3, ax=None):
    (ax or A).plot(xs, ys, color=color, lw=lw, ls=ls, zorder=z, solid_capstyle="butt")


def arrow(a, b, color=INDIGO, lw=0.9, rad=0.0, z=6, ax=None):
    (ax or A).add_patch(FancyArrowPatch(a, b, arrowstyle="-|>,head_length=1.8,head_width=1.0", color=color, lw=lw,
                                        zorder=z, shrinkA=0, shrinkB=0, mutation_scale=1,
                                        connectionstyle=f"arc3,rad={rad}"))


def panel(x, y, w, h, letter):
    """A panel's frame and letter (no title: the content starts below the letter's line, ~5 mm)."""
    rect(x, y, w, h, fc="white", ec=FRAME, lw=0.5, r=1.5, z=0)
    T(x + 1.6, y + 1.5, letter, 8, weight="bold", va="top")


def pic(img, x, y, w=None, h=None, z=3, border=None):
    """An image placed at (x, y) (top left, mm); give w or h, the aspect ratio is kept."""
    im = img if isinstance(img, Image.Image) else Image.open(img).convert("RGB")
    if w is None:
        w = h * im.width / im.height
    if h is None:
        h = w * im.height / im.width
    ax = fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], zorder=z)
    ax.imshow(np.asarray(im), interpolation="none")     # vector backends embed the image as is
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(border is not None)
        if border:
            sp.set_color(border); sp.set_linewidth(0.6)
    return w, h


def fit(art: str, x, y, w, h, z=3):
    """An illustration centred in the box (x, y, w, h); returns its placed (x, y, w, h) and the image."""
    im = Image.open(ART / f"{art}.png").convert("RGB")
    s = min(w / im.width, h / im.height)
    pw, ph = im.width * s, im.height * s
    px, py = x + (w - pw) / 2, y + (h - ph) / 2
    pic(im, px, py, pw, ph, z=z)
    return px, py, pw, ph, im


def spots_on(img: Image.Image, x, y, w, h, step=1.0, col=ROSE):
    """A grid of capture spots over the tissue of an illustration (its non-white pixels)."""
    a = np.asarray(img.convert("L"))
    for gx in np.arange(x + step / 2, x + w, step):
        for gy in np.arange(y + step / 2, y + h, step):
            px = int((gx - x) / w * a.shape[1]); py = int((gy - y) / h * a.shape[0])
            r0, r1 = max(py - 3, 0), min(py + 4, a.shape[0]); c0, c1 = max(px - 3, 0), min(px + 4, a.shape[1])
            if (a[r0:r1, c0:c1] < 238).all():
                TOP.add_patch(Circle((gx, gy), 0.25, fc=col, ec="none", alpha=0.55, zorder=8))


# --------------------------------------------------------------------------- a  one cell, measured three ways

#: the kinds of data: name, colour, illustration, what they measure
CELL_KINDS = [("Spatial Hi-C", ROSE, "a_tissue", "contacts per spot\n+ RNA, ATAC"),
              ("Single-cell Hi-C", TEAL, "a_cells", "contacts per cell\n+ RNA, ATAC"),
              ("Chromatin tracing", ORANGE, "a_nucleus", "x, y, z per locus\n+ IF, RNA")]
NUC_C, NUC_R = (506, 479), 147                                  # the nucleus in the 1024-px art of full Fig. 1a


def _nucleus_rgba():
    """The nucleus of the full Fig. 1a illustration, cut out as a round RGBA image."""
    a = np.asarray(Image.open(f1.HERE / "art3" / "a_flat_1.png").convert("RGBA")).copy()
    (cx, cy), r = NUC_C, NUC_R + 4
    sub = a[cy - r:cy + r, cx - r:cx + r].copy()
    yy, xx = np.mgrid[-r:r, -r:r]
    sub[..., 3] = np.where(xx ** 2 + yy ** 2 <= r * r, 255, 0)
    return sub



#: chromatin territories of the nucleus art (1024-px coordinates), coloured by the facet they stand for
TERRITORY = {"contacts": (463, 431), "geometry": (576, 458), "state": (502, 543)}
NUC_TURN = 60.0                                                  # degrees: geometry, contacts, state top to bottom


def _cell_image():
    """The nucleus, turned so that its geometry, contacts and state territories lie top to bottom (as the plates)."""
    im = Image.fromarray(_nucleus_rgba())
    return np.asarray(im.rotate(NUC_TURN, resample=Image.BICUBIC))


def _territory_xy(name, cx, cy, cr):
    """Where a territory lands in the placed cell (mm)."""
    tx, ty = TERRITORY[name]
    dx, dy = tx - NUC_C[0], -(ty - NUC_C[1])                     # y up
    t = np.deg2rad(NUC_TURN)
    rx, ry = dx * np.cos(t) - dy * np.sin(t), dx * np.sin(t) + dy * np.cos(t)
    k = cr / (NUC_R + 4)
    return np.array([cx + rx * k, cy - ry * k])


class _ShiftCanvas(f1.Canvas):
    """A full-figure Canvas whose source region (mm of the full Fig. 1) is drawn into a box of this figure."""

    def __init__(self, fig_, box, src):
        x, y, w, h = box
        self.fig, self.src = fig_, src
        self.ax = fig_.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], zorder=10)
        self.ax.axis("off")
        self.ax.patch.set_visible(False)
        self.fix()

    def fix(self):
        x0, x1, y0, y1 = self.src
        self.ax.set_xlim(x0, x1)
        self.ax.set_ylim(y1, y0)

    def letter(self, *a, **k):
        pass

    def title(self, *a, **k):
        return self.ax.text(0, 0, "")


def _stack(x, y, w, h):
    """The modality layers on shared genomic bins and the bridges between them: full Fig. 1b without its level
    strip and title, drawn as vectors into the box (x, y, w, h), aspect ratio kept."""
    import pandas as pd
    cell = pd.read_csv(f1.CACHE / "cell_1_0_34.csv")
    hic = np.load(f1.CACHE / "schicar_chr2.npy").astype(float)
    sx0, sx1, sy0, sy1 = 108.0, 171.0, 13.0, 98.5                  # the layers and bridges in full-figure mm
    k = min(w / (sx1 - sx0), (h - 4.0) / (sy1 - sy0))
    bw, bh = (sx1 - sx0) * k, (sy1 - sy0) * k
    bx, by = x + (w - bw) / 2, y + 4.0 + (h - 4.0 - bh) / 2
    mk = f1.Canvas(fig)                                              # the ChromData mark, then the heading
    mk.ax.set_zorder(21); mk.ax.patch.set_visible(False)
    f1.chromdata_mark(mk, bx, y - 1.3, 5.2)
    mk.fix()
    T(bx + 6.4, y + 2.0, "ChromData: layers on shared genomic bins", 6.2, INDIGO, weight="bold")
    keep = f1.level_strip, f1.chromdata_mark
    f1.level_strip = lambda *a, **kw: None
    f1.chromdata_mark = lambda *a, **kw: None
    # the plates' back-left edges (o, o + V), as f1.panel_b places them, in this figure's mm
    au, av, L_, D_, th, py0 = np.radians(21.0), np.radians(23.0), 28.5, 12.5, 1.0, 25.0
    U_ = np.array([np.cos(au), np.sin(au)]) * L_
    V_ = np.array([-np.cos(av), np.sin(av)]) * D_
    gap = (f1.A_BOTTOM - 4.2 - th - (U_[1] + V_[1]) - py0) / 2
    ox = 92.5 + 15.2 + 1.8 + np.cos(av) * D_
    to_mm = lambda q: np.array([bx + (q[0] - sx0) * k, by + (q[1] - sy0) * k])  # noqa: E731
    edges = {}
    for n_, key in enumerate(("geometry", "contacts", "state")):
        o = np.array([ox, py0 + n_ * gap])
        edges[key] = (to_mm(o), to_mm(o + V_))
    try:
        cv = _ShiftCanvas(fig, (bx, by, bw, bh), (sx0, sx1, sy0, sy1))
        f1.panel_b(cv, cell, hic)
        for t in list(cv.ax.texts):                               # the title line ("not a joint model") is not used
            if t.get_position()[0] < sx0:
                t.remove()
        cv.fix()
    finally:
        f1.level_strip, f1.chromdata_mark = keep
    return edges


def panel_a(x0=0.0, y0=0.0, w=118.0, h=87.0):
    """Three kinds of data, one cell, and the cell opening into ChromData's layers on shared genomic bins (full
    Fig. 1b, drawn from real data)."""
    from matplotlib.patches import Polygon
    panel(x0, y0, w, h, "a")
    for k, (name, col, art, what) in enumerate(CELL_KINDS):         # left: the kinds of data
        by = y0 + 8.0 + 26.0 * k
        px, py, pw, ph, im = fit(art, x0 + 2.0, by, 15.0, 11.5)
        if art == "a_tissue":
            spots_on(im, px, py, pw, ph)
        T(x0 + 2.0, by + 14.6, name, 5.8, col, weight="bold")
        T(x0 + 2.0, by + 16.9, what, 5.0, GREY, va="top")
    cx, cy, cr = x0 + 35.0, y0 + 45.0, 10.0                         # centre: the cell
    ax = fig.add_axes([(cx - cr) / W, 1 - (cy + cr) / H, 2 * cr / W, 2 * cr / H], zorder=12)
    ax.imshow(_cell_image(), interpolation="lanczos"); ax.axis("off")
    bx, b0, b1, r_ = x0 + 21.0, y0 + 11.0, y0 + 79.0, 1.2              # a brace: the three measure one cell
    t = np.linspace(0, np.pi / 2, 12)
    top = np.c_[bx - r_ + r_ * np.sin(t), b0 + r_ - r_ * np.cos(t)]
    bot = np.c_[bx - r_ + r_ * np.sin(t), b1 - r_ + r_ * np.cos(t)]
    mid = np.c_[bx + r_ - r_ * np.cos(t), cy - r_ + r_ * np.sin(t)]
    mid2 = np.c_[bx + r_ - r_ * np.cos(t), cy + r_ - r_ * np.sin(t)]
    for seg in (top, bot, mid, mid2):
        line(seg[:, 0], seg[:, 1], color=GREY, lw=0.6, z=30, ax=TOP)
    line([bx, bx], [b0 + r_, cy - r_], color=GREY, lw=0.6, z=30, ax=TOP)
    line([bx, bx], [cy + r_, b1 - r_], color=GREY, lw=0.6, z=30, ax=TOP)
    line([bx + r_, cx - cr - 0.4], [cy, cy], color=GREY, lw=0.6, z=30, ax=TOP)
    T(cx, cy + cr + 2.8, "one cell", 6.0, INK, ha="center", weight="bold", ax=TOP)
    edges = _stack(x0 + 46.0, y0 + 4.0, w - 47.0, h - 5.0)          # right: the layers
    for key, col in (("geometry", ORANGE), ("contacts", TEAL), ("state", PURPLE)):
        tx_, ty_ = _territory_xy(key, cx, cy, cr)
        e0, e1 = edges[key]                                         # the beam opens onto the plate and its label
        xe = e1[0] - 0.6
        top = e0[1] - {"geometry": 9.5, "contacts": 6.0, "state": 8.5}[key]
        A.add_patch(Polygon([(tx_, ty_ - 0.4), (xe, top), (xe, e1[1] + 1.8), (tx_, ty_ + 0.4)], closed=True,
                            fc=col, ec="none", alpha=0.13, zorder=1))
        TOP.add_patch(Circle((tx_, ty_), 0.75, fc="white", ec=col, lw=0.8, zorder=41))


# --------------------------------------------------------------------------- b  the browser (vector redraw)

REAL = np.load(SHARED / "agents_panels" / "_cache" / "takei_fov0.npz")      # shared/agents_panels/extract_data.py
CT = {"Granule": "#5b8db8", "Bergmann": "#5aa577", "Purkinje": "#c9745a", "MLI1": "#a58bc7", "MLI2+PLI": "#d4a24c",
      "Other": "#b9bdc4"}
LASSO = "Bergmann"
SHOW_CHROMS = {"chr7": "#dd6b2f", "chr15": "#3b63a8"}             # chromosomes drawn in colour in 3-D
HALO = [pe.withStroke(linewidth=1.5, foreground="white")]


def ax_mm(x, y, w, h):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], zorder=4)


def _clean(ax):
    ax.set_xticks([]), ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def _hull(pts, pad):
    """A smooth closed curve around pts (convex hull, padded, Chaikin-smoothed)."""
    from scipy.spatial import ConvexHull
    c = pts.mean(0)
    h = pts[ConvexHull(pts).vertices]
    h = h[np.argsort(np.arctan2(h[:, 1] - c[1], h[:, 0] - c[0]))]
    d = h - c
    h = c + d * (1 + pad / np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-9))
    for _ in range(4):
        q = np.roll(h, -1, 0)
        h = np.stack([0.75 * h + 0.25 * q, 0.25 * h + 0.75 * q], 1).reshape(-1, 2)
    return np.vstack([h, h[:1]])


#: label offsets per cell type, as fractions of the embedding's span
UMAP_LABEL = {"Granule": (0.0, 0.13), "MLI1": (0.0, 0.1), "MLI2+PLI": (-0.02, 0.1), "Purkinje": (0.0, -0.1),
              LASSO: (0.0, -0.17)}                                   # MLI1, Purkinje: 3 cells each, label beside


def view_umap(ax):
    u, ct = REAL["umap"], REAL["cell_type"]
    sel = ct == LASSO
    span = np.ptp(u, 0)
    hull = _hull(u[sel], pad=span.mean() * 0.07)
    ax.fill(hull[:, 0], hull[:, 1], color=CT[LASSO], alpha=0.12, lw=0, zorder=1)
    ax.plot(hull[:, 0], hull[:, 1], color=INK, lw=0.6, ls=(0, (3, 1.5)), zorder=2)
    for g in dict.fromkeys(ct):
        m = ct == g
        on = g == LASSO
        ax.scatter(u[m, 0], u[m, 1], s=13 if on else 9, color=CT.get(g, GREY), lw=0.5, edgecolors="white",
                   alpha=1 if on else 0.75, zorder=4 if on else 3)
        if g in UMAP_LABEL and m.any():
            c = u[m].mean(0) + np.array(UMAP_LABEL[g]) * span
            ax.text(*c, g, fontsize=5.0, color=CT[g], ha="center", va="center", weight="bold", zorder=6,
                    path_effects=HALO)
    right = hull[np.argmax(hull[:, 0])]                              # in the blank right of the lasso, its colour
    ax.text(right[0] + span[0] * 0.04, right[1] - span[1] * 0.02, "lasso", fontsize=5.0, color=CT[LASSO],
            ha="left", va="center", weight="bold", zorder=6, path_effects=HALO)
    pad = span * 0.1
    ax.set_xlim(u[:, 0].min() - pad[0] * 2.6, u[:, 0].max() + pad[0] * 2.2)
    ax.set_ylim(u[:, 1].min() - pad[1] * 1.6, u[:, 1].max() + pad[1])
    _clean(ax)


def _sphere(ax, c, r, col, alpha=1.0, z=1):
    """A shaded sphere: concentric discs, lighter towards the upper left."""
    base = np.array(mpl.colors.to_rgb(col))
    for k in range(9):
        f = k / 8
        cc = base + (1 - base) * (0.15 + 0.7 * f ** 1.5)
        ax.add_patch(Circle((c[0] - r * 0.28 * f, c[1] + r * 0.28 * f), r * (1 - 0.75 * f), fc=cc, ec="none",
                            alpha=alpha, zorder=z))


def view_nuclei(ax, aspect):
    """The nuclei around the selection, with their traced chromosomes (chr7 and chr15 in colour)."""
    xyz, cell, tr, chrom = REAL["xyz"], REAL["xyz_cell"], REAL["xyz_trace"], REAL["xyz_chrom"]
    ct = REAL["cell_type"]
    sel = np.flatnonzero(ct == LASSO)
    cen = np.array([np.nanmedian(xyz[cell == i, :2], 0) if (cell == i).any() else [np.nan, np.nan]
                    for i in range(len(ct))])
    sc = cen[sel]
    dd = np.linalg.norm(sc[:, None] - sc[None], axis=-1)
    ctr = sc[np.argmax((dd < 22).sum(1))]                          # the selected cell with most selected neighbours
    hy = 8.5
    hx = hy * aspect
    for i in range(len(ct)):
        m = cell == i
        if not m.any() or np.any(np.abs(cen[i] - ctr) > [hx * 1.3, hy * 1.3]):
            continue
        on = i in sel
        r = np.nanpercentile(np.linalg.norm(xyz[m, :2] - cen[i], axis=1), 90)
        _sphere(ax, cen[i], r * 1.1, CT[ct[i]] if on else "#c9ced6", alpha=0.5 if on else 0.35, z=1)
        for t in np.unique(tr[m]):
            mm = m & (tr == t)
            q = xyz[mm, :2].astype(float)
            col = SHOW_CHROMS.get(REAL["draw_chroms"][chrom[mm][0]]) if on else None
            step = np.linalg.norm(np.diff(q, axis=0), axis=1)
            if len(step):                                          # no bond across a gap > 10x the median step
                q = np.insert(q, np.flatnonzero(step > 10 * np.nanmedian(step)) + 1, np.nan, axis=0)
            if col:
                ax.plot(q[:, 0], q[:, 1], lw=0.35, color=col, alpha=0.65, zorder=4)
                ax.scatter(*xyz[mm, :2].T, s=0.8, color=col, lw=0, zorder=5)
            else:
                ax.scatter(*xyz[mm, :2].T, s=0.3, color="#9aa1ab" if on else "#c3c8cf", lw=0,
                           alpha=0.55 if on else 0.4, zorder=3 if on else 2)
    ax.set_xlim(ctr[0] - hx, ctr[0] + hx), ax.set_ylim(ctr[1] - hy, ctr[1] + hy)
    ax.axis("off")
    x0, y0 = ctr[0] - hx * 0.93, ctr[1] + hy * 0.62
    ax.add_patch(Rectangle((x0 - 0.5, y0 - 0.9), 3.0, 4.4, fc="white", ec="none", alpha=0.85, zorder=6))
    ax.plot([x0, x0 + 2], [y0, y0], color=INK, lw=0.9, solid_capstyle="butt", zorder=7)
    ax.text(x0 + 1, y0 + hy * 0.06, "2 µm", fontsize=5.0, color=INK, ha="center", va="bottom", zorder=7)
    for k, (n, col) in enumerate(list(SHOW_CHROMS.items())):
        ax.text(0.99, 0.05 + (1 - k) * 0.13, n, transform=ax.transAxes, fontsize=5.0, color=col, ha="right",
                va="bottom", weight="bold", zorder=6, path_effects=HALO)


def _triangle(ax, mat, base, depth, vmin, vmax):
    """HiGlass-style band of a symmetric matrix (up to `depth` bins apart), x in [0, n]."""
    i, j = np.triu_indices(mat.shape[0], 0)
    v = mat[i, j]
    ok = np.isfinite(v) & ((j - i) <= depth)
    polys = [[((a + b) / 2, base + (b - a) / 2), ((a + b + 1) / 2, base + (b - a) / 2 + 0.5),
              ((a + b + 2) / 2, base + (b - a) / 2), ((a + b + 1) / 2, base + (b - a) / 2 - 0.5)]
             for a, b in zip(i[ok], j[ok])]
    pc = PolyCollection(polys, array=v[ok], cmap="viridis_r", edgecolors="face", linewidths=0.05, zorder=3)
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)


def genome_window(n=30):
    """The n-bin window (1 Mb bins of chr15) with the most finite pixels in both groups' distance maps."""
    ma0, mb0 = REAL[f"dist_{LASSO}"], REAL["dist_Granule"]
    ok = np.isfinite(ma0) & np.isfinite(mb0)
    return max(range(ma0.shape[0] - n), key=lambda k: ok[k:k + n, k:k + n].sum())


def view_genome(ax, n=30, depth=10):
    """IF tracks (selection vs all cells) above the selection's median 3-D distances, on one genomic axis."""
    ma0 = REAL[f"dist_{LASSO}"]
    s0 = genome_window(n)
    x = np.arange(n) + 0.5
    rows = [("H3K27ac", ORANGE), ("H3K9me3", PURPLE), ("LaminB1", "#6b7280")]
    rh, gap = 2.6, 0.9
    base = 0.6
    top = base + depth / 2 + 1.4
    for k, (trk, col) in enumerate(rows[::-1]):
        y0 = top + k * (rh + gap)
        a, b = REAL[f"trk_{trk}_{LASSO}"][s0:s0 + n], REAL[f"trk_{trk}_all"][s0:s0 + n]
        both = np.r_[a[np.isfinite(a)], b[np.isfinite(b)]]
        lo, hi = np.nanpercentile(both, [2, 98])
        sc = lambda v: y0 + np.clip((v - lo) / max(hi - lo, 1e-9), 0, 1) * rh   # noqa: E731
        ax.fill_between(x, y0, sc(np.nan_to_num(a, nan=lo)), step="mid", color=CT[LASSO], alpha=0.55, lw=0)
        ax.step(x, sc(np.nan_to_num(b, nan=lo)), where="mid", color=INK, lw=0.45)
        ax.text(-0.4, y0 + rh / 2, trk, fontsize=5.0, color=col, ha="right", va="center", weight="bold")
    ma = ma0[s0:s0 + n, s0:s0 + n]
    vmin, vmax = np.nanpercentile(ma[np.isfinite(ma)], [3, 97])
    _triangle(ax, ma, base, depth, vmin, vmax)
    ax.text(-0.4, base + depth / 4, "3D\ndistance", fontsize=5.0, color=INK, ha="right", va="center",
            linespacing=0.95)
    ax.set_xlim(-8.5, n), ax.set_ylim(base - 0.4, top + len(rows) * (rh + gap))
    ax.axis("off")


def pane(x, y, w, h, num, title, scope=None, active=None):
    """A browser pane with a numbered one-line header; returns its content box (active: the pane the selection
    was made in, outlined)."""
    active = num == 1 if active is None else active
    rect(x, y, w, h, fc="white", ec=INDIGO if active else "#d5dae2", lw=0.7 if active else 0.4, r=0.6, z=3)
    rect(x + 0.3, y + 0.3, w - 0.6, 3.1, fc="#f5f6f8", ec="none", z=3)
    A.add_patch(Circle((x + 2.0, y + 1.85), 1.15, fc=INK, ec="none", zorder=5))
    T(x + 2.0, y + 1.9, str(num), 5.0, "white", ha="center", weight="bold", zorder=6)
    T(x + 4.0, y + 1.85, title, 5.2, INK, zorder=5)
    if scope:
        A.text(x + w - 1.2, y + 1.85, scope, fontsize=5.0, color=CT[LASSO], ha="right", va="center", zorder=5,
               bbox=dict(boxstyle="round,pad=0.2,rounding_size=0.5", fc="#e8f4ec", ec="none"))
    return x + 0.6, y + 3.6, w - 1.2, h - 4.0


def view_tissue(ax):
    """A piece of cerebellar cortex (illustration), drawing only what anatomy fixes: the molecular layer (ML) over
    the Purkinje cell layer (PCL) over the granular layer (GL); granule cells packed in the GL; Purkinje cells in a
    single row in the PCL with Bergmann glia somata beside them (more Bergmann glia than Purkinje cells); molecular
    layer interneurons (MLI1, MLI2) in the ML with no depth preference drawn, and Purkinje layer interneurons (in the
    MLI2+PLI group) also at the PCL. Not drawn: the white matter, and "Other" (a mixed group whose position is not
    fixed by its label). Colours as in the embedding; the lassoed cells marked as in (1) (ringed, the others
    faded)."""
    rng = np.random.default_rng(3)
    xs = np.linspace(0, 1, 200)
    at = lambda x, yb: yb + 0.10 * np.sin(np.pi * np.asarray(x))  # noqa: E731  (a gently bowed piece of cortex)
    for y0_, y1_, c in [(0.04, 0.42, "#e3d6e6"), (0.42, 0.52, "#eadbe2"), (0.52, 0.93, "#f5eef2")]:
        ax.fill_between(xs, at(xs, y0_), at(xs, y1_), color=c, lw=0, zorder=1)
    ax.plot(xs, at(xs, 0.93), color="#c9b3c2", lw=0.5, zorder=2)                 # pia

    def dots(x, yb, g, size):
        on = g == LASSO
        ax.scatter(x, at(x, yb), s=size * (1.6 if on else 1), color=CT[g], lw=0.5 if on else 0.3,
                   edgecolors=INK if on else "white", alpha=1 if on else 0.4, zorder=5 if on else 4)

    dots(rng.uniform(0.02, 0.98, 170), rng.uniform(0.06, 0.40, 170), "Granule", 3.2)
    px_ = np.linspace(0.07, 0.93, 6)
    dots(px_, np.full(6, 0.47), "Purkinje", 15)
    bx_ = np.concatenate([p + np.array([-0.05, -0.03, 0.04]) for p in px_])
    dots(bx_, 0.475 + rng.uniform(-0.015, 0.015, len(bx_)), "Bergmann", 5)
    n1, n2 = 13, 9                                                               # interneurons, intermixed
    dots(rng.uniform(0.04, 0.96, n1), rng.uniform(0.56, 0.90, n1), "MLI1", 5)
    dots(rng.uniform(0.04, 0.96, n2), rng.uniform(0.56, 0.90, n2), "MLI2+PLI", 5)
    pl = (px_[1:4] + px_[2:5]) / 2 + 0.012
    dots(pl, np.full(3, 0.455), "MLI2+PLI", 5)
    for g, x, yb in (("MLI1", 0.03, 0.85), ("MLI2+PLI", 0.03, 0.62), ("Purkinje", 0.40, 0.545),
                     ("Bergmann", 0.67, 0.545), ("Granule", 0.03, 0.12)):
        ax.text(x, at(x, yb), g, fontsize=4.2, color=CT[g], ha="left", va="center", weight="bold", zorder=7,
                path_effects=HALO)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0.0, 1.08)
    _clean(ax)


def panel_b(x0=121.0, y0=0.0, w=59.0, h=87.0):
    panel(x0, y0, w, h, "b")
    wx, wy, ww, wh = x0 + 2.0, y0 + 5.5, w - 4.0, h - 7.5            # a browser window
    rect(wx, wy, ww, wh, fc="white", ec="#9fb0cc", lw=0.6, r=1.0, z=1)
    rect(wx, wy, ww, 3.6, fc="#e9eef6", ec="none", z=2, r=1.0)
    for k, c in enumerate(["#ec6a5e", "#f4bf4f", "#61c554"]):
        A.add_patch(Circle((wx + 2.0 + 1.8 * k, wy + 1.8), 0.55, fc=c, ec="none", zorder=3))
    T(wx + ww - 1.5, wy + 1.8, "U-Chrom", 5.0, INDIGO, ha="right", weight="bold", zorder=3)
    px, pw = wx + 1.2, ww - 2.4
    ph = (wh - 4.6 - 2 * 1.0) / 3
    ys = [wy + 4.4 + k * (ph + 1.0) for k in range(3)]
    hw = (pw - 1.0) / 2                                              # top row: the embedding with the lasso | tissue
    c1 = pane(px, ys[0], hw, ph, 1, "Embedding · RNA UMAP", active=False)
    view_umap(ax_mm(c1[0] + 0.4, c1[1] + 0.2, c1[2] - 0.8, c1[3] - 0.4))
    c1t = pane(px + hw + 1.0, ys[0], hw, ph, 2, "Tissue · cell types", scope="lasso", active=False)
    view_tissue(ax_mm(c1t[0] + 0.3, c1t[1] + 0.2, c1t[2] - 0.6, c1t[3] - 0.4))
    c2 = pane(px, ys[1], pw, ph, 3, "3D · nuclei, traces", scope="lasso")
    view_nuclei(ax_mm(c2[0] + 0.3, c2[1] + 0.2, c2[2] - 0.6, c2[3] - 0.4), (c2[2] - 0.6) / (c2[3] - 0.4))
    s0 = genome_window()
    c3 = pane(px, ys[2], pw, ph, 4, f"Genome · chr15:{s0}–{s0 + 30} Mb", scope="lasso")
    view_genome(ax_mm(c3[0] + 0.4, c3[1] + 0.3, c3[2] - 0.8, c3[3] - 0.6))


# --------------------------------------------------------------------------- c  an open ecosystem

def _spatial_icon(cv, x, y, s):
    """Icon of uchrom.spatial in an s x s mm box at (x, y) (top-left): spots on a section, coloured by a smooth
    pattern."""
    n = 5
    for i in range(n):
        for j in range(n):
            u, v = (i + 0.5) / n, (j + 0.5) / n
            if (u - 0.5) ** 2 + (v - 0.5) ** 2 > 0.27:                 # the section's outline
                continue
            t = np.clip(0.15 + 0.85 * np.exp(-((u - 0.68) ** 2 + (v - 0.35) ** 2) / 0.08), 0, 1)
            col = mpl.colors.to_hex(np.array(mpl.colors.to_rgb(ROSE)) * t + np.array([0.93, 0.93, 0.94]) * (1 - t))
            cv.ax.add_patch(Circle((x + u * s, y + v * s), s * 0.075, fc=col, ec=f1.OUTLINE, lw=0.25, zorder=23.5))


def _group(cv, x, y, label, color, tags, lg=None, lg_h=f1.LOGO_H, plate=False):
    """As ``make_fig1_final.group`` ([logo] label over tags), but the logos sit on the tags' baseline: every item of
    a row ends at y + 7.0, the bottom of its tags."""
    blk = 7.0
    if lg:
        for key, h in ([(lg, lg_h)] if isinstance(lg, str) else lg):
            lw_, _ = f1.logo(cv, key, x, y + blk - h, h, z=22, plate=f1.LIGHTGREY if plate else None)
            x += lw_ + 1.8
    x_lab_end = x
    if label:
        t = cv.text(x + 0.2, y + 1.0, label, size=f1.F_SEC, color=color, weight="bold")
        x_lab_end = x + 0.2 + f1.measure(cv, t)
    ty = y + 2.6
    for t in tags:
        dashed, icon = False, None
        if isinstance(t, tuple):
            name, icon = t[0], (t[1] if len(t) > 1 else None)
            dashed = len(t) == 1 or (len(t) > 2 and t[2] == "linked")
        else:
            name = t
        if icon:
            ih = {"cooler": 3.0}.get(icon, 4.0)
            iw, _ = f1.logo(cv, icon, x, ty + 4.4 - ih, ih, z=22)
            x += iw + 0.9
        x = f1.tag(cv, x, ty, name, ec=color, fc=f1.PALE.get(color, "#f4f5f7"), dashed=dashed) + 1.1
    return max(x, x_lab_end + 1.0) + 3.0


def panel_c(x0=0.0, y0=90.0):
    """After the full Fig. 1c: interfaces, methods, the ChromData core, formats read / written / linked (dashed;
    the key is in the legend)."""
    panel(x0, y0, W, H - y0, "c")
    cv = f1.Canvas(fig)
    cv.ax.set_zorder(2); cv.ax.patch.set_visible(False)
    x_l, x_r = x0 + 2.0, W - 2.0
    xs0 = x_l + 21.0
    BAND = "#dfe3e8"

    def band(y, h, title, sub=None):
        cv.box(x_l, y, x_r - x_l, h, fc="white", ec=BAND, lw=0.5, r=1.6, z=2)
        cv.text(x_l + 2.0, y + h / 2 - (1.5 if sub else 0), title, size=6.5, weight="bold")
        if sub:
            cv.text(x_l + 2.0, y + h / 2 + 1.7, sub, size=5.5, color=GREY)

    y = y0 + 5.0
    h1 = 8.6                                                         # interfaces
    band(y, h1, "Interfaces", "people · agents")
    yc = y + h1 / 2
    import json
    at = json.loads((HERE.parents[1] / "benchmarks" / "atlas" / "results" / "atlas.json").read_text())["datasets"]
    n_ds, n_st = len(at), len({d["study"] for d in at})
    S, V = 4.6, 3.7                 # icon slot width S; every icon's drawn height V, centred on yc (measured)
    # no CLI: uchrom has no command of its own (a few `python -m` entry points only); the browser's agent
    # interface is WebMCP (Methods), not an MCP server
    names = ["Python API", "web browser, local or hosted", f"atlas: {n_ds} datasets, {n_st} studies",
             "WebMCP + agent skills"]
    widths = []
    for s_ in names:                                                 # the items spread evenly over the band
        t = cv.text(0, yc, s_, size=6.5, va="center")
        widths.append(f1.measure(cv, t))
        t.remove()
    gap_ = (x_r - 2.5 - xs0 - sum(S + 1.6 + w for w in widths)) / (len(names) - 1)

    def label(x, s_):
        t = cv.text(x + 1.6, yc + 0.13, s_, size=6.5, va="center")      # text centre on the icons' centre
        return x + 1.6 + f1.measure(cv, t) + gap_

    x = xs0
    snake = Image.open(f1.CACHE / "logos" / "python.png").convert("RGBA").crop((0, 0, 580, 713))  # the symbol
    sh = V * 4.6 / 3.81                                              # the file's drawn part is 3.81 of 4.6 mm
    sw = sh * 580 / 713
    pic(Image.alpha_composite(Image.new("RGBA", snake.size, "white"), snake).convert("RGB"),
        x + (S - sw) / 2, yc - sh / 2 + 0.36, w=sw, h=sh, z=5)
    x = label(x + S, names[0])
    st = V / 0.76                                                    # the browser box is 0.76 s tall
    f1.icon_browser(cv, x + (S - st) / 2, yc - st / 2, st, INDIGO)
    x = label(x + S, names[1])
    t_ = (V - 0.3) / 3                                               # the atlas: a 3 x 3 grid of dataset tiles
    for i in range(3):
        for j, c in enumerate((ROSE, TEAL, ORANGE)):
            if i == 2 and j == 0:
                continue                                             # 8 tiles
            cv.box(x + (S - V) / 2 + (t_ + 0.15) * i, yc - V / 2 + (t_ + 0.15) * j, t_, t_, fc=c, ec="white",
                   lw=0.2, r=0.2, z=22)
    x = label(x + S, names[2])
    mw, _ = f1.logo(cv, "mcp_icon", x + 0.3, yc - V / 2, V, z=22)
    label(x + S, names[3])
    y += h1 + 1.5
    h2 = 15.0                                                        # the analysis toolkit
    band(y, h2, "Analysis", "uchrom, one API")
    mods = [("recon", f1.SEQ, "3D structures\nfrom contacts"), ("im", f1.IMG, "spots → traces;\nmissing loci"),
            ("strc", f1.CORE, "TADs, loops,\ncompartments"), ("fea", f1.GREY, "per-locus and\nper-cell features"),
            ("emb", f1.OMX, "embeddings,\nclusters"), ("spatial", ROSE, "spatial tests,\nstripes removed")]
    gap = 1.6
    wd = (x_r - 1.5 - xs0 - gap * (len(mods) - 1)) / len(mods)
    for k, (name, c, txt) in enumerate(mods):
        xm = xs0 + k * (wd + gap)
        cv.box(xm, y + 1.3, wd, h2 - 2.6, fc=f1.PALE.get(c, "#f8eef2" if name == "spatial" else "#f4f5f7"), ec=c,
               lw=0.5, r=1.2, z=3)
        # each icon's drawn bottom on the baseline of its description's last line (offsets measured, mm)
        up = {"recon": 0.51, "im": 0.59, "strc": 0.34, "fea": 0.0, "emb": 0.17, "spatial": 0.85}[name]
        if name == "spatial":
            _spatial_icon(cv, xm + 0.9, y + 3.4 - up, 7.6)
        else:
            f1.mod_icon(cv, name, xm + 0.9, y + 3.4 - up, 7.6)
        cv.text(xm + 9.2, y + 4.6, "uchrom." + name, size=5.4, color=c if c != f1.GREY else "#5d636b",
                weight="bold")
        cv.text(xm + 9.2, y + 6.6, txt, size=5.2, color=INK, va="top", linespacing=1.15)
    y += h2 + 1.5
    h4 = H - 2.0 - y                                                 # formats, with ChromData in the middle
    band(y, h4, "Formats", "& standards")
    fx0 = x_l + 15.5                                                 # formats start right after the band label
    # ChromData between the formats it reads (left) and those it writes or links (right): the flow runs left to right
    cx0, cw_ = x_l + 70.5, 22.0
    cv.box(cx0, y + 1.2, cw_, h4 - 2.4, fc=f1.CORE, ec=f1.CORE, lw=0.5, r=1.4, z=4)
    ymid = y + h4 / 2                                                # the mark, the name and its file, centred
    f1.chromdata_mark(cv, cx0 + cw_ / 2 - 3.5, ymid - 9.6, 7.0)
    cv.text(cx0 + cw_ / 2, ymid + 1.6, "ChromData", size=8.0, color="white", weight="bold", ha="center", z=25)
    cv.text(cx0 + cw_ / 2, ymid + 5.8, ".chromdata.zarr", size=6.0, color="white", ha="center", z=25)
    blk = 7.0                                                        # a group: label row + tag row
    r1 = y + 1.8                                                     # three rows, evenly over the band
    pitch = (h4 - 3.6 - blk) / 2                                     # 8 mm: H is set for it
    r2, r3 = r1 + pitch, r1 + 2 * pitch

    def struct(x, yr, tag_):                                         # "structures" in the word slot, then its file
        tt = cv.text(x, yr + 2.6 + 2.2, "structures", size=6.5, color="#5d636b", weight="bold",
                     va="center_baseline")
        return f1.tag(cv, x + f1.measure(cv, tt) + 1.4, yr + 2.6, tag_, ec="#5d636b", fc="#f4f5f7")
    ya = r2 + 2.6 + 2.2                                              # the arrows run on the middle row's tag line
    xr0, xw = cx0 - 9.5, cx0 + cw_ + 8.5                            # end of the read zone, start of the write zone
    cv.arrow((xr0, ya), (cx0 - 0.2, ya), color=f1.CORE, lw=0.8, head=2.6)
    cv.text((xr0 + cx0) / 2, ya - 1.6, "read", size=5.0, color=GREY, weight="bold", ha="center")
    cv.arrow((cx0 + cw_ + 0.2, ya), (xw - 1.0, ya), color=f1.CORE, lw=0.8, head=2.6)
    cv.text((cx0 + cw_ + xw - 1.0) / 2, ya - 1.4, "write ·\nlink", size=5.0, color=GREY, weight="bold",
            ha="center", va="bottom", linespacing=1.0)
    LG = 5.0                                                         # one logo size down the read column
    x = _group(cv, fx0, r1, None, f1.IMG, ["FOF-CT"], lg="4dn", lg_h=LG)
    t = cv.text(x + 2.0, r1 + 2.6 + 2.2, "Juicer", size=7.5, color=f1.SEQ, weight="bold",   # the word as its
                va="center_baseline")                                                      # logo (none exists)
    f1.tag(cv, x + 2.0 + f1.measure(cv, t) + 1.8, r1 + 2.6, ".hic", ec=f1.SEQ, fc=f1.PALE[f1.SEQ])
    x = _group(cv, fx0, r2, None, f1.SEQ, [".pairs"], lg="open2c", lg_h=LG) - 3.0      # the logo names it
    f1.tag(cv, x + 1.1, r2 + 2.6, ".cool · .mcool · .scool", ec=f1.SEQ, fc=f1.PALE[f1.SEQ])
    x = _group(cv, fx0, r3, None, f1.IMG, ["pyHiM", "SnapFISH"], lg="pyhim",   # tracing software U-Chrom reads
               lg_h=LG / 1.08)                                       # same width as the square logos
    struct(x + 1.0, r3, ".3dg")
    zx, L, yc1 = xw, 2.2, r1 + 2.6 + 4.4 / 2                         # the store's own formats: one letter height
    for key, frac, mid in (("zarr", 0.538, 0.496), ("parquet", 0.856, 0.516)):  # word height / centre in the file
        fh = L / frac
        lw_, _ = f1.logo(cv, key, zx, yc1 - mid * fh, fh, z=22)
        zx += lw_ + 3.0
    _group(cv, zx + 3.0, r1, None, f1.IMG, ["FOF-CT"], lg="4dn", lg_h=LG)  # written as well as read (to_fofct)
    x = _group(cv, xw, r2, None, f1.OMX, [("AnnData", "anndata", "linked"), ("MuData", "mudata", "linked")],
               lg="scverse", lg_h=3.6)
    f1.tag(cv, x - 3.0 + 1.1, r2 + 2.6, "SpatialData", ec=ROSE, fc="#f9eef3", dashed=True)   # scverse too
    yt3 = r3 + 2.6                                                   # the Open2C logo, then its suffixes
    lw_, _ = f1.logo(cv, "open2c", xw, r3 + 7.0 - LG, LG, z=22)
    x = xw + lw_ + 1.8
    for name_ in (".mcool", ".scool"):
        x = f1.tag(cv, x, yt3, name_, ec=f1.SEQ, fc=f1.PALE[f1.SEQ], dashed=True) + 1.1
    struct(x + 3.0, r3, "PDB → PyMOL / ChimeraX")
    cv.fix()


def make_figure(stem: str = "fig1") -> Path:
    panel_a(); panel_b(); panel_c()
    for ax in (A, TOP):
        ax.set_xlim(0, W); ax.set_ylim(H, 0)
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"{stem}.{ext}", dpi=300)
    return HERE / f"{stem}.pdf"


if __name__ == "__main__":                                           # --out NAME writes NAME.pdf / NAME.png
    print(make_figure(sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else "fig1"))
