"""Fig. 2 -- the mouse cortex in three kinds of data:
single-cell and spatial Hi-C contacts, the same contacts reconstructed into 3D structures with U-Chrom's NucDynamics,
and chromatin tracing.

Fig. 2 (``fig2.pdf``), classes and loci:
  a  one data model, the same calls on every kind of data
  b  the class ranking of every dataset (rank of the class median of each cell's short- over long-range proximity),
     and ExN against Oligo (Cliff's delta, 95% CI; open, haploid reconstructions)
  c  where ExN and Oligo differ: the difference of the median distance of every pair of imaged loci on chr2 (corrected
     for genomic distance), imaging (upper triangle) against two reconstructions (lower)
  d  chr2 in one typical ExN and Oligo cell of each, captures of the U-Chrom browser's 3D view
     (``capture_fig2_structures.py``, run first)
  e  agreement of each reconstruction's difference map with imaging, against the upper limit the noise of the two
     maps allows (their split-half reliabilities)
  f  the same difference per locus along chr2: imaging, each reconstruction and their consensus
  g  every pair of classes: agreement with imaging, and imaging's own split-half agreement
(c, d and f, all chromosome 2, in the right column; b, e and g on the left.)

The single-cell results (the former Extended Data Fig. 4) are numbers of the Supplementary Note
(``results/cross_assay_ed_cells.json``).

Data: ``results/cross_assay_ed_{classes,loci}.json``, ``cross_assay_ed_{maps,tracks}.parquet`` (U-Chrom
``benchmarks/cross_assay/figdata.py``; provenance in the JSON).

    python figures/fig2/make_fig2.py      # -> fig2.pdf / .png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "benchmarks" / "cross_assay" / "results"
MM = 1 / 25.4
W = 180.0

mpl.rcParams.update({
    "pdf.fonttype": 42, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "font.size": 6, "axes.titlesize": 6, "axes.labelsize": 6,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "xtick.minor.size": 1, "ytick.minor.size": 1,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.0,
    "axes.labelpad": 1.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
})

INK, GREY, LGREY, BAND = "#2b2f36", "#8a8f98", "#c9ccd1", "#eceef1"
EXN, OLIGO, ROSE, INDIGO = "#b5507a", "#34568b", "#b5507a", "#34568b"
GROUP_COLOR = {"contacts": "#8a8f98", "reconstructions": "#b5507a", "imaging": "#34568b"}
GROUP_LABEL = {"contacts": "Contacts", "reconstructions": "Reconstructed", "imaging": "Imaging"}


def need(path: Path) -> Path:
    if not path.exists():
        raise SystemExit(f"{path} missing (see the docstring)")
    return path


def load(name: str) -> dict:
    d = json.loads(need(RESULTS / name).read_text())
    return d.get("results", d)


CL, LO = load("cross_assay_ed_classes.json"), load("cross_assay_ed_loci.json")
MAPS = pd.read_parquet(need(RESULTS / "cross_assay_ed_maps.parquet"))
TRACKS = pd.read_parquet(need(RESULTS / "cross_assay_ed_tracks.parquet"))
CLASSES = CL["classes"]
SHORT_LABEL = {"Dip-C (Tan), diploid": "Dip-C, diploid", "HiRES, diploid": "HiRES, diploid",
               "GAGE-seq, haploid": "GAGE-seq", "dscHi-C, haploid": "dscHi-C", "Droplet Hi-C, haploid": "Droplet Hi-C",
               "scHiCAR, haploid": "scHiCAR"}


def blocks(M, k: int) -> np.ndarray:
    """Each value as a k x k block of pixels: a viewer that scales the embedded image up (Preview smooths it) cannot
    blur one cell into the next."""
    return np.kron(np.asarray(M, float), np.ones((k, k)))


class Canvas:
    def __init__(self, h: float):
        self.H = h
        self.fig = plt.figure(figsize=(W * MM, h * MM))
        self.A = self.fig.add_axes([0, 0, 1, 1])
        self.A.set_xlim(0, W); self.A.set_ylim(h, 0); self.A.axis("off")

    def T(self, x, y, s, size=6, color=INK, ha="left", va="center", weight="normal", **kw):
        return self.A.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, fontweight=weight, **kw)

    def head(self, x, y, letter):
        """The panel letter, alone (the journal's figures carry no panel titles: they are in the legend)."""
        self.T(x, y, letter, 8, weight="bold", va="top")

    def ax(self, x, y, w, h):
        return self.fig.add_axes([x / W, 1 - (y + h) / self.H, w / W, h / self.H])

    def box(self, x, y, w, h, text, color=GREY, size=5.0, fill="white"):
        self.A.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.0", lw=0.6,
                                        edgecolor=color, facecolor=fill))
        self.T(x + w / 2, y + h / 2, text, size, ha="center", va="center", linespacing=1.15)

    def arrow(self, x0, y0, x1, y1, color=GREY):
        self.A.annotate("", xy=(x1, y1), xytext=(x0, y0),
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color=color, shrinkA=0, shrinkB=0, mutation_scale=5))


# =========================================================================== Extended Data Fig. 3


CODE_COLOR = {"mod": "#7a828e", "fn": INDIGO, "str": "#1f8a7d", "p": INK, "arg": "#b0573a"}
ART = HERE.parents[1] / "shared" / "art"


def code_chip(C, x, y, w, h, lines, size=5.0):
    """Code block as in Fig. 2a: rounded light-grey box, tokens coloured by kind, the font shrunk to fit.  Tokens are
    placed by character count times the font's advance, measured once on a long string (monospaced, so exact; the
    extent of a single token is its ink and would drop trailing punctuation)."""
    C.A.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.9", lw=0, fc="#eef0f3"))
    r = C.fig.canvas.get_renderer()
    mm_per_px = W / C.fig.bbox.width
    probe = C.A.text(0, 0, "0" * 50, fontsize=10, family="DejaVu Sans Mono")
    adv10 = probe.get_window_extent(renderer=r).width * mm_per_px / 50      # mm per character at 10 pt
    probe.remove()
    n = max(sum(len(tok) for tok, _ in line) for line in lines)
    sz = min(size, (w - 2.0) / (n * adv10 / 10))
    adv = adv10 * sz / 10
    lh = sz * 0.42
    y0 = y + h / 2 - (len(lines) - 1) * lh / 2
    for k, line in enumerate(lines):
        cx = x + 1.0
        for tok, kind in line:
            C.A.text(cx, y0 + k * lh, tok, fontsize=sz, family="DejaVu Sans Mono", color=CODE_COLOR[kind], va="center")
            cx += len(tok) * adv                        # (regular weight only: the bold face advances less here)

A_STEPS = [
    ("contacts", "8 datasets", "a2_contacts", GREY,
     [[("uc.datasets.", "mod"), ("load", "fn"), ("(", "p")], [("  \"gageseq_…\"", "str"), (")", "p")]]),
    ("3D from contacts", "6 datasets", "x3_recon", ROSE,
     [[("uc.recon.sc.", "mod")], [("reconstruct_nucdyn", "fn"), ("(…)", "p")]]),
    ("imaged traces", "2 datasets", "x3_trace", INDIGO,
     [[("ChromData.", "mod"), ("read", "fn"), ("(", "p")], [("  \"liu2025_mop…\"", "str"), (")", "p")]]),
    ("short vs long range", "b", "x3_separation", INK,
     [[("uc.fea.", "mod")], [("separation_profile", "fn"), ("(", "p"), ("cd", "arg"), (")", "p")]]),
    ("where classes differ", "c-g", "x3_classes", INK,
     [[("uc.fea.", "mod"), ("distance_map", "fn"), ("(", "p")], [("  cd", "arg"), (", \"chr2\", cells=…)", "p")]]),
]


#: the separation icon redrawn as vectors (``art/x3_separation.png``, pixels, y down): the bead centres along the chain,
#: and the two pairs whose distance it compares, adjacent beads (short) and the chain's far ends (long), as dashed arcs
SEP_BEADS = [(23.0, 156.3), (21.7, 120.0), (41.8, 89.2), (60.2, 58.8), (110.3, 34.2), (135.2, 57.9), (118.1, 85.0),
             (86.9, 95.9), (76.2, 126.6), (93.0, 157.2), (125.8, 168.3), (156.4, 146.8), (188.0, 142.5),
             (209.7, 165.0), (230.6, 183.7)]
SEP_SHORT, SEP_LONG = (3, 4), (4, 14)


def draw_separation(ax, box):
    """The chain of beads with a dashed teal arc over two adjacent beads and a dashed blue arc between far-apart ones;
    ``box`` = (left, top, right, bottom) of the icon in its pixels."""
    from matplotlib.patches import Circle, FancyArrowPatch
    P = np.array(SEP_BEADS)
    ax.plot(P[:, 0], P[:, 1], color=INK, lw=0.9, zorder=1, solid_capstyle="round")
    for (i, j), col, rad in ((SEP_SHORT, "#1f8a7d", -0.8), (SEP_LONG, INDIGO, -0.55)):   # ends just off the beads
        ax.add_patch(FancyArrowPatch(P[i], P[j], arrowstyle="-", connectionstyle=f"arc3,rad={rad}", color=col,
                                     lw=1.0, ls=(0, (2.2, 1.3)), shrinkA=3.2, shrinkB=3.2, zorder=4))
    for p in P:
        ax.add_patch(Circle(p, 11.0, fc="#d4d6da", ec=INK, lw=0.7, zorder=3))
    ax.set_xlim(box[0], box[2]); ax.set_ylim(box[3], box[1]); ax.set_aspect("equal")


def panel_a(C, x0, y0, w, h):
    from PIL import Image

    C.head(x0, y0, "a")
    gap, big = 2.6, 7.0
    bw = (w - 1.0 - (len(A_STEPS) - 2) * gap - big) / len(A_STEPS)
    top, bh = y0 + 3.5, h - 4.5
    x = x0 + 0.5
    for i, (title, note, icon, col, code) in enumerate(A_STEPS):
        C.A.add_patch(FancyBboxPatch((x, top), bw, bh, boxstyle="round,pad=0,rounding_size=1.2", lw=0.7, ec=col,
                                     fc="white"))
        C.T(x + 1.2, top + 2.2, title, 5.6, col if col != GREY else INK, weight="bold")
        C.T(x + bw - 1.0, top + 2.2, note, 4.4, GREY, ha="right")
        im = Image.open(need(ART / f"{icon}.png")).convert("RGB")
        ink = np.asarray(im).min(2) < 235                                # the icon without its white margin
        r_, c_ = np.flatnonzero(ink.any(1)), np.flatnonzero(ink.any(0))
        crop = (max(c_[0] - 3, 0), max(r_[0] - 3, 0), min(c_[-1] + 4, im.width), min(r_[-1] + 4, im.height))
        im = im.crop(crop)
        area = bh - 13.0                                                 # between the title and the code
        ih = area
        iw = min(ih * im.width / im.height, bw - 4.0)
        ih = iw * im.height / im.width
        ax = C.ax(x + (bw - iw) / 2, top + 4.0 + (area - ih) / 2, iw, ih)
        if icon == "x3_separation":
            draw_separation(ax, crop)
        else:
            ax.imshow(np.asarray(im), interpolation="lanczos")
        ax.axis("off")
        if icon == "a2_contacts":
            ax.text(1.0, -0.05, "short", fontsize=4.4, color="#1f8a7d", ha="right", va="top", transform=ax.transAxes)
            ax.text(0.64, 0.82, "long", fontsize=4.4, color=INDIGO, ha="left", va="center", transform=ax.transAxes)
        if icon == "x3_separation":                                     # the teal and the blue arc
            ax.text(0.20, 0.97, "short", fontsize=4.4, color="#1f8a7d", ha="right", va="center", transform=ax.transAxes)
            ax.text(1.03, 0.55, "long", fontsize=4.4, color=INDIGO, ha="left", va="center", transform=ax.transAxes)
        if code is not None:
            code_chip(C, x + 0.8, top + bh - 8.2, bw - 1.6, 7.0, code)
        else:
            C.T(x + bw / 2, top + bh - 4.7, "a locus signature\nscored on each cell", 4.6, GREY, ha="center",
                linespacing=1.2)
        nx = x + bw + (big if i == 2 else gap)
        if i == 2:          # the three kinds of data, each one ChromData, into the same calls
            C.A.annotate("", (nx - 0.6, top + bh / 2), (x + bw + 0.6, top + bh / 2),
                         arrowprops=dict(arrowstyle="-|>,head_length=0.8,head_width=0.45", lw=1.3, color=GREY,
                                         shrinkA=0, shrinkB=0))
        elif 3 <= i < len(A_STEPS) - 1:
            C.A.annotate("", (nx - 0.3, top + bh / 2), (x + bw + 0.3, top + bh / 2),
                         arrowprops=dict(arrowstyle="-|>,head_length=0.5,head_width=0.3", lw=0.9, color=GREY,
                                         shrinkA=0, shrinkB=0))
        x = nx

def panel_b(C, x0, y0, w, h):
    C.head(x0, y0, "b")
    y0, h = y0 - 4.0, h + 4.0  # no title line: the content moves up
    rows = CL["rows"]
    lab_w, cw, gap_g = 38.0, 7.0, 1.6                                # room for the group labels and bars
    n = len(rows)
    rh = (h - 15.5 - 2 * gap_g) / n                                  # room below for the axis title
    cmap = plt.get_cmap("RdPu_r")
    y = y0 + 10.0
    prev = None
    for r in rows:
        if prev is not None and r["group"] != prev:
            y += gap_g
        prev = r["group"]
        label = r["label"] + (" (pooled)" if r.get("pooled_only") else "")
        r["_t"] = C.T(x0 + lab_w - 1.0, y + rh / 2, label, 4.8, ha="right")
        r["_y0"] = y
        for k, c in enumerate(CLASSES):
            rk = r["ranks"].get(c)
            xx = x0 + lab_w + k * cw
            if rk is None:
                C.A.add_patch(Rectangle((xx, y), cw - 0.4, rh - 0.4, color=BAND, lw=0))
                continue
            C.A.add_patch(Rectangle((xx, y), cw - 0.4, rh - 0.4, color=cmap((rk - 1) / 6 * 0.85), lw=0))
            C.T(xx + (cw - 0.4) / 2, y + (rh - 0.4) / 2, str(rk), 4.6, "white" if rk <= 3 else INK, ha="center")
        r["_y"] = y + rh / 2
        y += rh
    rend, inv = C.fig.canvas.get_renderer(), C.A.transData.inverted()
    left = lambda t: inv.transform(t.get_window_extent(rend))[:, 0]  # noqa: E731  (a text's x extent, mm)
    xl = min(left(r["_t"]).min() for r in rows) - 1.2                # one bar per group, all on one line left of
    for g in dict.fromkeys(r["group"] for r in rows):                # the longest label; the group's name centred
        rr = [r for r in rows if r["group"] == g]                     # on its bar
        yt, yb_ = rr[0]["_y0"] + 0.15, rr[-1]["_y0"] + rh - 0.55
        C.A.plot([xl, xl], [yt, yb_], color=GROUP_COLOR[g], lw=0.9, solid_capstyle="butt")
        t = C.T(x0 + 0.5, (yt + yb_) / 2, GROUP_LABEL[g], 4.8, GROUP_COLOR[g], weight="bold")
        assert left(t).max() < xl - 0.8, f"{g}: label reaches the bar"
    for k, c in enumerate(CLASSES):
        C.T(x0 + lab_w + k * cw + (cw - 0.4) / 2, y0 + 8.2, c, 5.0, ha="center")
    ax0 = x0 + lab_w + len(CLASSES) * cw + 3.0
    aw = x0 + w - ax0 - 1.0
    ax = C.ax(ax0, y0 + 10.0, aw, y - (y0 + 10.0))
    ax.set_ylim(y, y0 + 10.0)
    ax.set_xlim(-0.05, 1.05)
    ax.axvline(0, color=LGREY, lw=0.5)
    for r in rows:
        e = r["cliff"] or {}
        if "cliffs_delta" in e:
            col = GROUP_COLOR[r["group"]]
            ax.plot(e["ci95"], [r["_y"]] * 2, color=col, lw=0.8)
            ax.scatter(e["cliffs_delta"], r["_y"], s=7, color=col, lw=0, zorder=3)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xticks([0, 0.5, 1]); ax.tick_params(labelsize=4.8)
    ax.set_xlabel("Cliff's delta, ExN vs Oligo", fontsize=4.9, labelpad=1.5)
    # the rank key and Kendall's W are in the legend


def _matrix(src: str, col: str = "resid"):
    t = MAPS[MAPS["source"] == src]
    loci = np.unique(np.r_[t["i"].to_numpy(), t["j"].to_numpy()])
    pos = {l: k for k, l in enumerate(loci)}
    M = np.full((len(loci), len(loci)), np.nan)
    M[[pos[i] for i in t["i"]], [pos[j] for j in t["j"]]] = t[col].to_numpy()
    return M, loci


def _on(M, loci, target):
    """M re-indexed to the loci ``target`` (NaN where absent)."""
    pos = {l: k for k, l in enumerate(loci)}
    idx = np.array([pos.get(l, -1) for l in target])
    out = np.full((len(target), len(target)), np.nan)
    ok = idx >= 0
    out[np.ix_(ok, ok)] = M[np.ix_(idx[ok], idx[ok])]
    return out



PANEL = pd.read_parquet(need(RESULTS / "cross_assay_ed_panel.parquet"))


def _panel_xyz(t: pd.DataFrame, src: str) -> np.ndarray:
    """A chromosome 2 structure on the imaged loci (NaN where absent): imaging by its own loci, a reconstruction by
    the particle holding each locus's midpoint ((k - 1, k] Mb)."""
    loci = PANEL[PANEL["chrom"].astype(str) == "chr2"].sort_values("start")
    mid = ((loci["start"] + loci["end"]) // 2).to_numpy()
    X = np.full((len(loci), 3), np.nan)
    t = t.sort_values("start")
    xyz = t[["x", "y", "z"]].to_numpy(float)
    if src == "imaging":
        pos = {int(v): k for k, v in enumerate(loci["start"])}
        for s0, v in zip(t["start"], xyz):
            if int(s0) in pos:
                X[pos[int(s0)]] = v
    else:
        st, en = t["start"].to_numpy(), t["end"].to_numpy()
        k = np.searchsorted(en, mid, side="left")
        ok = (k < len(en)) & (st[np.minimum(k, len(en) - 1)] < mid)
        X[ok] = xyz[k[ok]]
    return X


def _superpose(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Q moved onto P over the loci both have: scaled to the same radius of gyration (a least-squares scale shrinks a
    poorly matching structure), then the least-squares rotation (a mirror allowed: contacts do not fix chirality)."""
    m = np.isfinite(P).all(1) & np.isfinite(Q).all(1)
    p0, q0 = P[m].mean(0), Q[m].mean(0)
    A, B = P[m] - p0, Q[m] - q0
    B = B * np.sqrt((A ** 2).sum() / (B ** 2).sum())
    U, _, Vt = np.linalg.svd(B.T @ A)
    sc = np.sqrt((A ** 2).sum() / ((Q[m] - q0) ** 2).sum())
    return (Q - q0) * sc @ (U @ Vt) + p0


STRUCT_SOURCES = ["imaging", "Dip-C (reconstructed)", "scHiCAR (reconstructed)"]


def structure_views() -> dict:
    """The example structures of panel c in the frame they are drawn in: {(cls, source): (n_loci, 3)} on the imaged
    chr2 loci (NaN where absent), each reconstruction superposed onto the imaged cell of its class, all in that
    cell's principal axes (x, y on the page, z towards the viewer; centred on the imaged cell), in um."""
    S = pd.read_parquet(need(RESULTS / "cross_assay_ed_structures.parquet"))
    out = {}
    for cls in ("ExN", "Oligo"):
        img = _panel_xyz(S[(S["source"] == "imaging") & (S["cls"] == cls)], "imaging")
        ok = np.isfinite(img).all(1)
        c0 = img[ok].mean(0)
        _, _, vt = np.linalg.svd(img[ok] - c0, full_matrices=False)
        vt[2] *= np.sign(np.linalg.det(vt))                              # a rotation (right-handed frame)
        for src in STRUCT_SOURCES:
            X = img if src == "imaging" else _superpose(img, _panel_xyz(S[(S["source"] == src) & (S["cls"] == cls)], src))
            out[cls, src] = (X - c0) @ vt.T
    return out


CAPTURES = HERE / "captures"


def _capture_crop(key: str):
    """A browser capture of a structure (``capture_fig2_structures.py``) cropped to what is drawn, and where it lies in
    um (left, right, bottom, top): placed by its camera (centre and pixels per um at the box centre), so every capture
    shares one scale."""
    m = json.loads((CAPTURES / "fig2_structures.json").read_text())[key]
    img = plt.imread(CAPTURES / m["png"])
    H, W = img.shape[:2]
    ink = np.flatnonzero((img[..., :3] < 0.97).any(-1).any(1)), np.flatnonzero((img[..., :3] < 0.97).any(-1).any(0))
    i0, i1 = max(ink[0][0] - 4, 0), min(ink[0][-1] + 5, H)
    j0, j1 = max(ink[1][0] - 4, 0), min(ink[1][-1] + 5, W)
    (cx, cy), s = m["center_um"], m["px_per_um"]
    return img[i0:i1, j0:j1], (cx + (j0 - W / 2) / s, cx + (j1 - W / 2) / s, cy + (H / 2 - i1) / s, cy + (H / 2 - i0) / s)


def _capture(ax, key: str, lim: tuple):
    """The capture on the axes in um (``lim``: half the width and height shown)."""
    img, ext = _capture_crop(key)
    ax.imshow(img, extent=ext, interpolation="antialiased")
    ax.set_xlim(-lim[0], lim[0]); ax.set_ylim(-lim[1], lim[1]); ax.set_aspect("equal")
    ax.axis("off")


def rainbow(t):
    """The browser's colour along a chromosome (frontend/src/colors.ts ``rainbow``): violet at the start, red at the end."""
    return mpl.colors.hsv_to_rgb(np.c_[(1 - np.clip(t, 0, 1)) * 0.75, np.full(len(t), 0.85), np.ones(len(t))])


def panel_c(C, x0, y0, w, h):
    C.head(x0, y0, "c")
    y0, h = y0 - 4.0, h + 4.0  # no title line: the content moves up
    up, loci = _matrix("imaging (Liu 2025)")
    up = up / np.nanquantile(np.abs(up), 0.97)           # each source on its own scale (97th percentile = 1)
    s = (w - 6.0) / 2
    n = len(loci)
    iu = np.triu_indices(n, 1)
    for k, (src, lab) in enumerate((("Dip-C (Tan), diploid", "Dip-C (reconstructed)"),
                                    ("scHiCAR, haploid", "scHiCAR (reconstructed)"))):
        lo = _on(*_matrix(src), loci)
        lo = lo / np.nanquantile(np.abs(lo), 0.97)
        M = np.full((n, n), np.nan)
        M[iu] = up[iu]
        M[(iu[1], iu[0])] = lo[iu]
        xx = x0 + 0.5 + k * (s + 5.0)
        ax = C.ax(xx, y0 + 9.5, s, s)
        ax.imshow(blocks(M, 6), cmap="RdBu_r", vmin=-1, vmax=1, interpolation="none",
                  extent=(-0.5, M.shape[1] - 0.5, M.shape[0] - 0.5, -0.5))
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_linewidth(0.4); sp.set_color(LGREY)
        C.T(xx + s, y0 + 8.0, "imaging (DNA-MERFISH)", 4.5, INDIGO, ha="right")
        C.T(xx, y0 + 9.5 + s + 1.8, lab, 4.5, ROSE)
    cy = y0 + 9.5 + s + 4.6
    cax = C.ax(x0 + 3.0, cy, 14.0, 1.1)
    mpl.colorbar.ColorbarBase(cax, cmap=plt.get_cmap("RdBu_r"), norm=mpl.colors.Normalize(-1, 1),
                              orientation="horizontal", ticks=[-1, 0, 1])
    cax.tick_params(length=1.5, pad=1, labelsize=4.2)
    C.T(x0 + 18.6, cy + 0.55, "log distance, ExN \u2212 Oligo", 4.5, INK)     # beside the bar; the rest: legend


D_SOURCES = ["imaging", "Dip-C (reconstructed)", "scHiCAR (reconstructed)"]


def panel_d(C, x0, y0, w, h):
    """The three kinds of structure, one typical ExN and one Oligo cell each (chromosome 2), on one scale: each
    capture in a window around what it draws, the scale the largest at which the widest fits its column and both rows
    fit the panel; spare height shared by the two rows."""
    C.head(x0, y0, "d")
    y0, h = y0 - 2.0, h + 2.0  # no title line: the content moves up
    sy = y0 + 5.0
    cols = [("imaging", "imaging (DNA-MERFISH)", INDIGO), ("Dip-C (reconstructed)", "Dip-C (reconstructed)", ROSE),
            ("scHiCAR (reconstructed)", "scHiCAR (reconstructed)", ROSE)]
    lab_w = 3.5                                                        # the class names, rotated
    tw = (w - lab_w) / 3
    aw = tw - 0.8
    for j, (src, lab, col) in enumerate(cols):
        C.T(x0 + lab_w + j * tw + tw / 2, sy + 1.0, lab, 4.3, col, ha="center")
    if not (RESULTS / "cross_assay_ed_structures.parquet").exists() or not (CAPTURES / "fig2_structures.json").exists():
        C.T(x0 + w / 2, sy + 10.0, "(structures: benchmarks/cross_assay/structures.py, capture_fig2_structures.py)", 4.4,
            GREY, ha="center")
        return
    ext = {(cls, src): np.array(_capture_crop(f"{cls}|{src}")[1]) for cls in ("ExN", "Oligo") for src in D_SOURCES}
    yb = y0 + h - 1.0                                                  # the bottom line: 1 um, colours
    y, gap = sy + 2.3, 1.0
    avail = yb - 1.0 - y - gap
    rows_um = [max(e[3] - e[2] for (c, _), e in ext.items() if c == cls) for cls in ("ExN", "Oligo")]
    mm = min(min(aw / (e[1] - e[0]) for e in ext.values()), avail / sum(rows_um))   # mm per um, one for all six
    spare = (avail - sum(rows_um) * mm) / 2
    for cls, r in zip(("ExN", "Oligo"), rows_um):
        ah = r * mm + spare
        for j, (src, lab, _) in enumerate(cols):
            e = ext[cls, src]
            cx, cy = (e[0] + e[1]) / 2, (e[2] + e[3]) / 2
            ax = C.ax(x0 + lab_w + j * tw + 0.4, y, aw, ah)
            img, ex = _capture_crop(f"{cls}|{src}")
            ax.imshow(img, extent=ex, interpolation="antialiased")
            ax.set_xlim(cx - aw / mm / 2, cx + aw / mm / 2); ax.set_ylim(cy - ah / mm / 2, cy + ah / mm / 2)
            ax.set_aspect("equal"); ax.axis("off")
        C.T(x0 + 1.2, y + ah / 2, cls, 4.8, EXN if cls == "ExN" else OLIGO, weight="bold", rotation=90, ha="center")
        y += ah + gap
    C.A.plot([x0 + lab_w + 1.0, x0 + lab_w + 1.0 + mm], [yb, yb], color=INK, lw=0.8, solid_capstyle="butt")
    C.T(x0 + lab_w + 1.8 + mm, yb, "1 \u00b5m", 4.0, INK)
    bx, by = x0 + w - 8.0, yb - 0.45                                   # colour along chromosome 2 (the browser's)
    cax = C.ax(bx - 10.0, by, 10.0, 0.9)
    cax.imshow(rainbow(np.linspace(0, 1, 256))[None], aspect="auto", extent=(0, 1, 0, 1))
    cax.set_xticks([]); cax.set_yticks([])
    for sp in cax.spines.values():
        sp.set_visible(False)
    C.T(bx - 10.6, by + 0.45, "chr2 start", 4.0, GREY, ha="right", va="center")
    C.T(bx + 0.6, by + 0.45, "end", 4.0, GREY, va="center")


def panel_e(C, x0, y0, w, h):
    C.head(x0, y0, "e")
    y0, h = y0 - 4.0, h + 4.0  # no title line: the content moves up
    rows = LO["agreement"]
    ax = C.ax(x0 + 23.0, y0 + 7.5, w - 25.0, h - 16.0)
    sb = lambda r: 2 * r / (1 + r)  # noqa: E731
    for k, r in enumerate(rows):
        L, cr, cl = r["loci_liu"], r["ceiling_recon"], r["ceiling_liu"]
        ceil = np.sqrt(sb(cl["resid"]) * sb(cr["resid"])) if cr.get("resid") else np.nan
        ax.plot([0, ceil], [k, k], color=BAND, lw=4.0, solid_capstyle="butt", zorder=1)
        ax.scatter(L["resid"], k, s=14, color=ROSE, lw=0, zorder=3)
        ax.scatter(L["pair"], k, s=12, facecolor="white", edgecolor=ROSE, lw=0.6, zorder=3)
        if r.get("contacts"):                            # the same cells' contacts, without a 3D model
            ax.scatter(r["contacts"]["loci_liu"]["resid"], k, s=13, marker="x", color=INK, lw=0.6, zorder=4)
        ax.text(1.0, k, f"{100 * r['disattenuated_resid']:.0f}%", fontsize=4.6, va="center", ha="right", color=INK)
    ax.set_yticks(range(len(rows)), [r["label"] for r in rows], fontsize=4.8)
    ax.set_ylim(len(rows) - 0.5, -0.6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xticks([0, 0.25, 0.5, 0.75])
    ax.set_xlabel("Spearman with imaging over pairs of loci", fontsize=5.0, labelpad=1.5)
    ax.tick_params(axis="x", labelsize=4.8, length=1.5, pad=1)
    ax.text(1.0, -0.7, "of the upper limit", fontsize=4.3, ha="right", va="bottom", color=GREY)
    hs = [mpl.lines.Line2D([], [], ls="", marker="o", ms=3, color=ROSE, mew=0),
          mpl.lines.Line2D([], [], ls="", marker="o", ms=3, mfc="white", mec=ROSE, mew=0.6),
          mpl.lines.Line2D([], [], ls="", marker="x", ms=3, color=INK, mew=0.6),
          Rectangle((0, 0), 1, 1, color=BAND)]
    ax.legend(hs, ["all pairs, corrected for genomic distance", "pair-specific part",
                   "contacts, no 3D model", "upper limit (noise)"],
              loc="lower center", bbox_to_anchor=(x0 + w / 2, y0 + h + 2.6), bbox_transform=C.A.transData, ncol=4,
              frameon=False, fontsize=4.5,
              handletextpad=0.3, columnspacing=0.8)



def panel_f(C, x0, y0, w, h):
    C.head(x0, y0, "f")
    y0, h = y0 - 5.0, h + 5.0  # no title line: the content moves up (and the key line: in the legend)
    t = TRACKS[TRACKS["chrom"] == "chr2"].sort_values("mid")
    srcs = ["imaging (Liu 2025)", "consensus of reconstructions"] + list(SHORT_LABEL)
    labels = ["imaging", "consensus"] + [r["label"] for r in LO["agreement"]]
    lab_w, rho_w = 21.0, 7.5
    top = y0 + 9.5
    th = (h - 9.5 - 6.0) / len(srcs)
    rho = LO["locus_effect_vs_imaging"]
    for k, (s, lab) in enumerate(zip(srcs, labels)):
        ax = C.ax(x0 + lab_w, top + k * th, w - lab_w - rho_w, th * 0.86)
        v = t[s].to_numpy()
        v = (v - np.nanmean(v)) / np.nanstd(v)
        x = t["mid"].to_numpy() / 1e6
        col = INDIGO if k == 0 else ROSE if k == 1 else GREY
        ax.fill_between(x, 0, v, where=v > 0, color=col, lw=0, alpha=0.9 if k < 2 else 0.6)
        ax.fill_between(x, 0, v, where=v < 0, color=col, lw=0, alpha=0.35 if k < 2 else 0.25)
        ax.axhline(0, color=LGREY, lw=0.3)
        ax.set_xlim(x.min(), x.max()); ax.set_ylim(-3.2, 3.2)
        ax.set_yticks([]); ax.spines["left"].set_visible(False)
        if k < len(srcs) - 1:
            ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
        else:
            ax.set_xlabel("chromosome 2 (Mb)", fontsize=5.0, labelpad=1.0)
            ax.tick_params(labelsize=4.6, length=1.5, pad=1)
        C.T(x0 + lab_w - 1.0, top + k * th + th * 0.43, lab, 4.5, col if k < 2 else INK, ha="right",
            weight="bold" if k < 2 else "normal")
        if k > 0:
            key = "consensus of reconstructions" if k == 1 else s
            C.T(x0 + w - 0.5, top + k * th + th * 0.43, f"ρ {rho[key]:.2f}", 4.4, GREY, ha="right")
    C.T(x0 + w - 0.5, top - 1.2, "ρ, all loci", 4.2, GREY, ha="right")


def panel_g(C, x0, y0, w, h):
    C.head(x0, y0, "g")
    y0, h = y0 - 4.0, h + 4.0  # no title line: the content moves up
    P = pd.DataFrame(LO["class_pairs"])
    order = ["ExN-Oligo", "InN-Oligo", "ExN-Micro", "InN-Micro", "ExN-Astro", "InN-Astro", "Astro-Oligo", "Astro-Micro",
             "Oligo-Micro", "OPC-Micro", "OPC-Oligo", "InN-ExN", "ExN-OPC", "InN-OPC", "OPC-Astro"]
    P = P.reindex(order)
    rows = [c for c in P.columns if c != "imaging split-half"] + ["imaging split-half"]
    P = P[rows].T
    ax = C.ax(x0 + 24.0, y0 + 7.0, w - 33.0, h - 17.0)
    A = P.to_numpy(float)
    im = ax.imshow(blocks(A, 12), cmap="RdPu", vmin=0, vmax=0.75, aspect="auto", interpolation="none",
                   extent=(-0.5, A.shape[1] - 0.5, A.shape[0] - 0.5, -0.5))
    ax.set_yticks(range(len(rows)), rows, fontsize=4.6)
    ax.set_xticks(range(len(order)), [o.replace("-", "–") for o in order], fontsize=4.4, rotation=45, ha="right",
                  rotation_mode="anchor")
    ax.tick_params(length=0, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.axhline(len(rows) - 1.5, color="white", lw=1.2)
    cx = x0 + w - 7.5
    cax = C.ax(cx, y0 + 7.0, 1.1, h - 17.0)
    plt.colorbar(im, cax=cax, orientation="vertical", ticks=[0, 0.25, 0.5, 0.75])
    cax.tick_params(length=1.5, pad=1, labelsize=4.4)
    cax.set_ylabel("Spearman with imaging", fontsize=4.5, color=GREY, labelpad=2)



def make_ed3() -> Path:
    yb = 33.5                                                       # panel letters only, no titles (version of
    ye = yb + 56.0 + 1.5                                            # 2026-10-14): e under b
    top = ye + 24.0 + 3.0                                           # f and g side by side, tops aligned; e's key
    hfg = 42.0                                                      # between
    C = Canvas(top + hfg)
    panel_a(C, 0.5, 0.0, W - 1.0, 31.5)
    panel_b(C, 0.5, yb, 102.0, 56.0)
    panel_c(C, 106.0, yb, W - 106.5, 46.0)             # c, d and f, all chromosome 2, in the right column
    panel_d(C, 106.0, yb + 47.0, W - 106.5, ye + 24.0 - yb - 47.0)  # its bottom with e's
    panel_e(C, 0.5, ye, 102.0, 24.0)
    panel_f(C, 106.0, top, W - 106.5, hfg)
    panel_g(C, 0.5, top, 102.0, hfg)
    print(f"Fig. 2: {top + hfg:.1f} mm tall (f and g from {top:.1f} mm)")
    for ext in ("pdf", "png"):
        C.fig.savefig(HERE / f"fig2.{ext}", dpi=300)
    return HERE / "fig2.pdf"

if __name__ == "__main__":
    print(make_ed3())
