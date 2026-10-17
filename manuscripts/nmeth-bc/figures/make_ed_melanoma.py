"""Brief Communication, Figure 2 -- the 3D genome across a tumour section: copy number, read from the contacts,
explains much of the spatial variation of the contacts that the RNA leaves (against random covariates as smooth in
space; the mouse brain is Extended Data Fig. 3).  Spatial
Hi-C-RNA of a human melanoma (Guo et al. 2026, section r1; r2 replicates) from the atlas, analysed with U-Chrom.
Every panel has its own letter; the explanations and statistics are in the legend (main.tex).

  a  the workflow, one call per step, each with what it gives on this section
  b  the U-Chrom web browser on the section, a copy-number patch selected, five linked views (numbered 1-5;
     ``captures/ed_melanoma_browser.png``, ``capture_ed_melanoma.py``)
  c  the eight contact features of each spot on the tissue (destriped), beside b
  d  how much of each feature's spatial pattern the spot's RNA and its copy number explain, both sections
  e  what the RNA leaves of the short / long ratio = the part copy number explains + the rest; the patch outlined
(2026-10-14: panel letters only; step 5 of a is the spatial test; f, which repeated d with the null beside it, is
gone: the null comparison is a number of the legend.)

The patch is an illustration of a selection, not a finding: against random patches of its size and shape it stands
out in r1 but not in r2 (``benchmarks/spatial_variable`` stage ``patch_null``; Supplementary Note).

Data: ``results/spatial_variable_{melanoma,brain}{.json,_spots.parquet}`` (``benchmarks/spatial_variable``,
stages through patch + summary; provenance in the JSON).

    python manuscripts/nmeth-bc/figures/make_ed_melanoma.py      # -> ed_melanoma.pdf / ed_melanoma.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle
from PIL import Image
from scipy.ndimage import gaussian_filter

HERE = Path(__file__).resolve().parent
CAPTURES = HERE / "captures"
RESULTS = HERE / "results"
MM = 1 / 25.4
W, H = 180.0, 183.5                                   # panel letters only (2026-10-14; 188 with titles)
SEC, SEC2 = "human_melanoma_r1", "human_melanoma_r2"

mpl.rcParams.update({
    "pdf.fonttype": 42, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "font.size": 6, "axes.titlesize": 6, "axes.labelsize": 6,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "xtick.minor.size": 1, "ytick.minor.size": 1,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.0,
    "axes.labelpad": 1.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    # glyph advances as the font gives them (hinting snaps them to whole pixels in the PNG: the code tokens of
    # panel a, placed one after another, would drift)
    "text.hinting": "none", "text.hinting_factor": 1,
})

INK, GREY = "#2b2f36", "#8a8f98"
INDIGO, ROSE = "#34568b", "#b5507a"
PATCH = "#9b3ef0"                 # the browser's colour of a saved selection
RNA_C, CNV_C, LEFT_C = "#8fa8c8", "#d98cb0", "#eceef1"
GROUPS = {"tumour": ["Melanoma_Invasive", "Melanoma_Melanocytic", "Melanoma_NeuralCrest", "Melanoma_Transitory"],
          "T / NK": ["Activated CD8", "Deletional tolerance CD8", "Dysfunc CD8", "Precursor exhaustion", "Tcm", "Tem",
                     "Terminal exhaustion", "Tfh", "Treg", "NK"],
          "B / plasma": ["B cells memory", "B cells naive", "Plasmablasts"],
          "myeloid": ["MDSC", "Macrophage M1", "Macrophage M2", "Monocytes", "cDC", "pDC"],
          "stroma": ["CAFs MSC iCAF-like", "CAFs myCAF-like", "PVL differentiated", "PVL immature"],
          "endothelial": ["Endothelial"]}
GROUP_COLOR = {"tumour": ROSE, "T / NK": INDIGO, "B / plasma": "#5aa0d8", "myeloid": "#e07a3f", "stroma": "#7a9a3a",
               "endothelial": "#d9a441"}
#: the features of panels c and d (U-Chrom names) and their labels
FEATURES = [("short_long_ratio", "short / long ratio"), ("frac_short", "short range"),
            ("frac_long", "long range"), ("frac_2M_12M", "2-12 Mb"),
            ("ps_slope_100k_1M", "decay 0.1-1 Mb"), ("ps_slope_1M_10M", "decay 1-10 Mb"),
            ("cis_fraction", "cis fraction"), ("window:boundary_strength_oe", "TAD insulation")]

fig = plt.figure(figsize=(W * MM, H * MM))
A = fig.add_axes([0, 0, 1, 1])
A.set_xlim(0, W); A.set_ylim(H, 0); A.axis("off")


def T(x, y, s, size=6, color=INK, ha="left", va="center", weight="normal", **kw):
    return A.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, fontweight=weight, **kw)


def head(x, y, letter):
    """The panel letter, alone (the journal's figures carry no panel titles; what they said is in the legend)."""
    T(x, y, letter, 8, weight="bold", va="top")


def ax_mm(x, y, w, h):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def need(path: Path) -> Path:
    if not path.exists():
        raise SystemExit(f"{path} missing (see the docstring)")
    return path


S = json.loads(need(RESULTS / "spatial_variable_melanoma.json").read_text())
S = S.get("results", S)
B = json.loads(need(RESULTS / "spatial_variable_brain.json").read_text())
B = B.get("results", B)
SPOTS = pd.read_parquet(need(RESULTS / "spatial_variable_melanoma_spots.parquet"))
R1 = SPOTS[SPOTS.section == SEC]


def fmt_p(p: float) -> str:
    if p <= 0:
        return "P < 1e-300"
    e = int(np.floor(np.log10(p)))
    return f"P = {p / 10 ** e:.1f}e{e}" if e < -2 else f"P = {p:.2g}"


def grid(t: pd.DataFrame, col) -> np.ndarray:
    x, y = t["x"].to_numpy(int), t["y"].to_numpy(int)
    img = np.full((y.max() + 1, x.max() + 1), np.nan)
    img[y, x] = t[col].to_numpy(float) if isinstance(col, str) else np.asarray(col, float)
    return img


def smoothed(img: np.ndarray, sigma: float = 2.0) -> np.ndarray:
    ok = np.isfinite(img)
    num = gaussian_filter(np.where(ok, img, 0.0), sigma, mode="constant")
    den = gaussian_filter(ok.astype(float), sigma, mode="constant")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(ok, num / den, np.nan)


def show(ax, img, cmap="viridis", lim=None):
    v = img[np.isfinite(img)]
    lo, hi = lim if lim is not None else np.quantile(v, [0.02, 0.98])
    ax.imshow(img, cmap=cmap, vmin=lo, vmax=hi, interpolation="none")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    return lo, hi


def patch_outline(ax, t=R1):
    """The patch's outline, closed along the edges of the section."""
    m = np.pad(np.nan_to_num(grid(t, t["in_patch"].astype(float))), 1)
    sm = gaussian_filter(m, 1.0)
    ny, nx = m.shape
    ax.contour(np.arange(nx) - 1, np.arange(ny) - 1, sm, levels=[0.5], colors=PATCH, linewidths=0.8)
    ax.set_xlim(-0.5, nx - 2.5); ax.set_ylim(ny - 2.5, -0.5)


def split(sec: str, feat: str, res: dict = S):
    """Shares of a feature's spatial pattern: explained by the RNA, more by copy number, unexplained (stage cnv for
    the contact features, stage local for the pooled window scores: the same smooth variance relative to depth)."""
    if feat.startswith("window:"):
        v = res["local_per_spot"][f"{sec}/{feat.split(':')[1]}"]
        rna, both = v["depth_rna30.smooth_left"], v["depth_rna30_cnv.smooth_left"]
    else:
        left = res["cnv"][sec]["smooth_left"]
        rna, both = left["rna30"][f"q30.{feat}"], left["rna30+cnv"][f"q30.{feat}"]
    return 1 - rna, rna - both, both


def column(feat: str) -> str:
    """The spots table's column of a feature (destriped)."""
    return f"window.{feat.split(':')[1]}" if feat.startswith("window:") else f"destriped.q30.{feat}"


# --------------------------------------------------------------------------- a

ART = HERE / "art"
#: code tokens: (text, kind); kinds coloured as in an editor
CODE_COLOR = {"mod": "#7a828e", "fn": INDIGO, "str": "#1f8a7d", "p": INK, "arg": "#b0573a"}


def code_chip(x, y, w, h, lines, size=5.2):
    """A code block: rounded light-grey box, one or two lines of tokens coloured by kind, the font shrunk to fit.
    DejaVu Sans Mono advances every glyph (regular and bold) by 0.602 em, so the tokens are placed by their length."""
    A.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.9", lw=0, fc="#eef0f3"))
    em_mm = 25.4 / 72 * 0.602                               # one glyph, per point of font size
    n = max(sum(len(tok) for tok, _ in line) for line in lines)
    sz = min(size, (w - 2.4) / (n * em_mm))
    lh = sz * 0.42                                          # line height, mm
    y0 = y + h / 2 - (len(lines) - 1) * lh / 2
    for k, line in enumerate(lines):
        cx = x + 1.2
        for tok, kind in line:
            A.text(cx, y0 + k * lh, tok, fontsize=sz, family="DejaVu Sans Mono", color=CODE_COLOR[kind], va="center",
                   fontweight="bold" if kind == "fn" else "normal")
            cx += len(tok) * sz * em_mm


STEPS = [("open the section", "a2_atlas",
          [[("uc.datasets.", "mod"), ("load", "fn"), ("(", "p")], [("  \"hicrna_human_melanoma\"", "str"), (")", "p")]]),
         ("contacts by distance", "a2_contacts",
          [[("cd.cellm[", "p")], [("  \"contacts.distance\"", "str"), ("]", "p")]]),
         ("features of each spot", "a2_features",
          [[("uc.fea.", "mod"), ("contact_features", "fn"), ("(", "p"), ("cd", "arg"), (")", "p")]]),
         ("barcode stripes out", "a2_destripe",
          [[("uc.spatial.", "mod"), ("destripe", "fn"), ("(", "p"), ("f", "arg"), (", ...)", "p")]]),
         ("spatially variable?", "a2_spatial",                     # a pattern, tested against a null
          [[("uc.spatial.", "mod"), ("spatially_variable", "fn"), ("(", "p")],
           [("  ", "p"), ("cd", "arg"), (", ", "p"), ("f", "arg"), (", covariates=", "p"), ("X", "arg"), (")", "p")]])]


def panel_a(y0):
    head(0.5, y0, "a")
    gap, x = 5.0, 0.8
    w = (W - 2 * x - 4 * gap) / 5
    top, bh = y0 + 3.5, 27.5
    for i, (title, icon, code) in enumerate(STEPS):
        col = INDIGO if i < 4 else ROSE
        A.add_patch(FancyBboxPatch((x, top), w, bh, boxstyle="round,pad=0,rounding_size=1.2", lw=0.7, ec=col,
                                   fc="white"))
        T(x + 1.4, top + 2.3, f"{i + 1}  {title}", 6.0, col, weight="bold")
        im = Image.open(need(ART / f"{icon}.png")).convert("RGB")
        ih = 12.5
        iw = min(ih * im.width / im.height, w - 4.0)
        ih = iw * im.height / im.width
        ax = ax_mm(x + (w - iw) / 2, top + 4.6 + (12.5 - ih) / 2, iw, ih)
        ax.imshow(np.asarray(im), interpolation="lanczos")
        ax.axis("off")
        if icon == "a2_contacts":                          # the two ranges of the icon's triangle (on its axes: on top)
            ax.text(-0.03, 0.12, "short", fontsize=4.8, color="#1f8a7d", ha="right", va="center", transform=ax.transAxes)
            ax.text(0.64, 0.80, "long", fontsize=4.8, color=INDIGO, ha="left", va="center", transform=ax.transAxes)
        code_chip(x + 1.0, top + bh - 9.0, w - 2.0, 7.6, code)
        if i < len(STEPS) - 1:
            A.annotate("", (x + w + gap - 0.5, top + bh / 2), (x + w + 0.5, top + bh / 2),
                       arrowprops=dict(arrowstyle="-|>,head_length=0.8,head_width=0.45", lw=1.3, color=GREY,
                                       shrinkA=0, shrinkB=0))
        x += w + gap
    return top + bh + 1.0


# --------------------------------------------------------------------------- b

def panel_b(y0, w=118.0):
    head(0.5, y0, "b")
    im = Image.open(need(CAPTURES / "ed_melanoma_browser.png")).convert("RGB")
    x0 = 1.0
    h = (w - 0.6) * im.height / im.width
    top = y0 + 3.5
    A.add_patch(FancyBboxPatch((x0, top), w, h + 4.2, boxstyle="round,pad=0,rounding_size=1.5", lw=0.6, ec=GREY,
                               fc="#f4f5f7"))
    for k, c in enumerate(("#ec6a5e", "#f4bf4f", "#61c554")):
        A.add_patch(mpl.patches.Circle((x0 + 2.4 + k * 2.3, top + 2.1), 0.7, color=c, lw=0))
    T(x0 + w / 2, top + 2.1, "U-Chrom web browser  ·  hicrna_human_melanoma, section r1", 5.0, GREY, ha="center")
    ax = ax_mm(x0 + 0.3, top + 4.2, w - 0.6, h)
    ax.imshow(np.asarray(im), interpolation="none")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    # the five views, numbered for the legend: each number at the left end of its view's header (where the view's
    # type icon is), from the views' places on the page (captures/ed_melanoma_browser.json)
    lay = json.loads(need(CAPTURES / "ed_melanoma_browser.json").read_text())
    k = lay["scale"]
    for num, v in enumerate(lay["views"], start=1):
        cx, cy = (v["x"] + 13) * k, (v["y"] + v["head"] / 2) * k
        ax.plot([cx], [cy], "o", ms=6.5, color=INK, mec="white", mew=0.6, zorder=5)
        ax.text(cx, cy, str(num), fontsize=5.0, color="white", ha="center", va="center", fontweight="bold", zorder=6)
    return top + h + 4.2


# --------------------------------------------------------------------------- c, d

def panel_c(x0, y0, w, h):
    """The eight features on the tissue, two columns by four rows beside the browser (as large as its height
    allows)."""
    head(x0, y0, "c")
    ncol, gap, lab, key = 2, 1.5, 4.6, 5.0                          # key: the colour bar below
    s_ = min((w - (ncol - 1) * gap) / ncol, (h - 3.5 - 4 * lab - key) / 4)
    xs = x0 + (w - ncol * s_ - (ncol - 1) * gap) / 2
    for k, (f, label) in enumerate(FEATURES):
        r, c = divmod(k, ncol)
        x, y = xs + c * (s_ + gap), y0 + 3.5 + r * (s_ + lab)
        show(ax_mm(x, y, s_, s_), grid(R1, column(f)))
        T(x + s_ / 2, y + s_ + 1.8, label, 4.6, INK, ha="center")
    yk = y0 + 3.5 + 4 * (s_ + lab) + 0.6                              # one bar for all: each map is scaled to its
    bw_ = s_ * 0.9                                                    # own 2nd-98th percentile (show)
    cax = ax_mm(xs + s_ + gap / 2 - bw_ / 2, yk, bw_, 1.2)
    mpl.colorbar.ColorbarBase(cax, cmap=plt.get_cmap("viridis"), orientation="horizontal", ticks=[])
    cax.set_frame_on(False)
    T(xs + s_ + gap / 2 - bw_ / 2 - 1.0, yk + 0.6, "low", 4.6, INK, ha="right")
    T(xs + s_ + gap / 2 + bw_ / 2 + 1.0, yk + 0.6, "high", 4.6, INK)
    T(xs + s_ + gap / 2, yk + 3.4, "each feature on its own scale", 4.4, GREY, ha="center")


def panel_d(x0, y0, w, h):
    """Per feature, sections r1 (left bar) and r2: the share of its spatial pattern the spot's RNA explains, what its
    copy number adds and what is left.  (The comparison with random covariates as smooth as the copy-number
    components, stage cnv_null, is a number of the legend: drawn on the bars it was not readable.)"""
    head(x0, y0, "d")
    ax = ax_mm(x0 + 9.0, y0 + 5.0, w - 9.5, h - 15.0)
    bw_ = 0.36
    for k, (f, _) in enumerate(FEATURES):
        for j, sec in enumerate((SEC, SEC2)):
            a_, b_, c_ = split(sec, f)
            xx = k + (j - 0.5) * 0.40
            for lo, v, col in ((0, a_, RNA_C), (a_, b_, CNV_C), (a_ + b_, c_, LEFT_C)):
                ax.bar(xx, v, bottom=lo, width=bw_, color=col, lw=0)
    ax.set_xlim(-0.6, len(FEATURES) - 0.4)
    ax.set_ylim(0, 1); ax.set_yticks([0, 0.5, 1], ["0", "50", "100 %"])
    ax.tick_params(axis="y", labelsize=4.6, length=1.5, pad=1)
    ax.set_ylabel("share of the pattern", fontsize=4.8, labelpad=1.0)
    ax.set_xticks(range(len(FEATURES)), [lab for _, lab in FEATURES], rotation=35, ha="right", fontsize=4.6,
                  rotation_mode="anchor")
    ax.tick_params(axis="x", length=0, pad=1.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for xx, sec in ((-0.2, "r1"), (0.2, "r2")):                     # the sections, over the first pair of bars
        ax.text(xx, 1.02, sec, fontsize=4.2, color=GREY, ha="center", va="bottom")
    hs = [Rectangle((0, 0), 1, 1, color=c) for c in (RNA_C, CNV_C, LEFT_C)]
    ax.legend(hs, ["RNA", "+ copy number", "unexplained"], loc="lower left",
              bbox_to_anchor=(x0 + 9.0, y0 + 2.6), bbox_transform=A.transData, ncol=4, frameon=False, fontsize=4.6,
              handlelength=0.9, handletextpad=0.3, columnspacing=0.9, borderaxespad=0)


# --------------------------------------------------------------------------- e, f, g

def panel_e(x0, y0, w_all):
    """What the RNA leaves of the short / long ratio = the part copy number explains + the rest; the patch outlined."""
    head(x0, y0, "e")
    col = "q30.short_long_ratio"
    after_rna = smoothed(grid(R1, f"rna30.{col}"))
    after_both = smoothed(grid(R1, f"rna30_cnv.{col}"))
    maps = [(after_rna, "unexplained by RNA"), (after_rna - after_both, "copy number"), (after_both, "the rest")]
    lim = np.nanquantile(np.abs(after_rna), 0.98)
    gap = 3.0
    w = (w_all - 1.0 - 2 * gap) / 3
    x = x0 + 1.0
    for i, (img, title) in enumerate(maps):
        ax = ax_mm(x, y0 + 5.5, w, w)
        show(ax, img, "RdBu_r", lim=(-lim, lim))
        patch_outline(ax)
        T(x + w / 2, y0 + 4.0, title, 5.0, INK, ha="center")
        if i < 2:
            T(x + w + gap / 2, y0 + 5.5 + w / 2, "=" if i == 0 else "+", 7, GREY, ha="center")
        x += w + gap
    cax = ax_mm(x0 + 1.0 + w + gap + 3.0, y0 + 5.5 + w + 1.8, w - 6.0, 1.3)
    mpl.colorbar.ColorbarBase(cax, cmap=plt.get_cmap("RdBu_r"), norm=mpl.colors.Normalize(-lim, lim),
                              orientation="horizontal", ticks=[-lim, 0, lim], format="%.2f")
    cax.tick_params(length=1.5, pad=1, labelsize=4.6)
    T(x0 + 1.0 + w + gap + 3.0 + w - 6.0 + 1.5, y0 + 5.5 + w + 2.45, "short / long ratio left after RNA", 4.6, INK)
    return y0 + 5.5 + w + 6.0


#: the copy-number share against random covariates with the same spatial autocorrelation (U-Chrom
#: benchmarks/spatial_variable stage cnv_null, Moran spectral randomisation, the "bands" null)
NULL = json.loads(need(RESULTS / "spatial_variable_cnv_null_melanoma.json").read_text())["results"]["nulls"]["bands"]


def null_of(sec: str, feat: str) -> dict:
    key = feat if feat.startswith("window:") else f"q30.{feat}"
    return NULL["sections"][sec]["features"][key]


def make_figure() -> Path:
    y = panel_a(0.0)
    yb = y + 1.5
    wb = 118.0
    b_end = panel_b(yb, wb)
    xr, wr = wb + 4.5, W - wb - 5.5
    panel_c(xr, yb, wr, b_end - yb)                       # the features beside the browser
    yr = b_end + 3.0                                      # below: d, what explains each pattern, and e, one of them
    xe = 93.0
    y = panel_e(xe, yr, W - xe - 0.5)
    panel_d(0.5, yr, xe - 3.0, y - yr)
    print(f"figure height used: {y:.1f} of {H} mm")
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"ed_melanoma.{ext}", dpi=300)
    return HERE / "ed_melanoma.pdf"


if __name__ == "__main__":
    print(make_figure())
