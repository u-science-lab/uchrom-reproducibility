"""Extended Data Fig. 3 (mouse brain) -- the same analysis as Fig. 2 on the mouse brain sections
of Guo et al. 2026 (spatial Hi-C-RNA, sections r1 and r2): the contact features vary across the section as its cells
do.  Every panel has its own letter; the explanations and statistics are in the legend below.

  a  the authors' spatial domains, section r1 (the only section they published them for), named from marker genes
  b  the eight contact features of each spot on the tissue (destriped), section r1
  c  the short / long ratio by domain, section r1
  d  how much of each feature's spatial pattern the spot's RNA explains, both sections
  e  the yardstick: what the same RNA adjustment leaves of the features and of held-out genes and gene modules
  f  the short / long ratio after depth = the part the RNA explains + the rest, section r1
  g  each domain's departure from the section mean of the ratio, after depth and after depth + RNA, section r1

Legend (draft).  **The 3D genome across a mouse brain section follows its cells.**  Spatial Hi-C-RNA of a coronal
mouse brain section (Guo et al. 2026; sections r1 and r2, 6,050 and 5,766 spots after quality control), opened from
the U-Chrom atlas and analysed as in Fig. 2.
**a**, The authors' spatial domains (CellCharter on the RNA, published for section r1), named from their marker
genes (cortex L2-4 Cux2 / Rorb / Satb2; L5-6 Bcl11b / Tshz2; CA1 Bcl11b / Wfs1 / Fibcd1; CA3 Grik4; DG Prox1;
neuropil -- cortical layer 1 and the hippocampal strata -- Ddn / Slc1a2; striatum Pde10a / Gpr88 / Penk; thalamus
Tcf7l2 / Slc17a6 / Ntng1; habenula Pou4f1 / Tac2 / Gpr151; white matter Plp1 / Mbp; choroid plexus Ttr).
**b**, Eight contact features of each spot, destriped; TAD insulation against the spot's own expected.
**c**, The short / long ratio by domain (boxes: quartiles; whiskers: 5th-95th percentiles).  Cortex, CA and DG lie
above white matter and choroid plexus (medians 3.22 vs 2.42; two-sided Mann-Whitney P < 1e-300), as Guo et al.
report for the cell types of these regions (their cell-type labels are not public).
**d**, The share of each feature's spatial pattern (the variance of a smooth surface fitted to the spots) that the
spot's RNA (30 principal components) explains; sections r1 (upper bar) and r2.  RNA explains 80-90 % (TAD
insulation 72 and 59 %, the 0.1-1 Mb decay slope 46 and 42 %).  Copy-number covariates (Fig. 2;
the chromosomes' shares of the contacts, U-Chrom ``uchrom.fea.chromosome_dosage``) add at most 1.2 %.
**e**, The yardstick: genes detected in >= 10 % of the spots are split at random in two halves; the RNA components
are computed from one half, and the other half (held-out genes, grey) is adjusted by them as the contact features
(rose) are.  Modules (open circles) average held-out genes with similar spatial patterns.  x: the share of a
value's variance (after depth) in its smooth spatial pattern; y: the share of that pattern left after RNA.  The
weaker a value's pattern, the more of it the RNA leaves.  The features keep 7-15 % (the 0.1-1 Mb decay slope, the
weakest pattern, 43 and 55 %), held-out genes with a pattern (x > 20 %) 25 and 28 % (median); the module with the
strongest pattern in each section (112 and 35 genes, x = 49 and 47 %) keeps 11 and 10 %, as much as the features.
**f**, The short / long ratio after depth, split into the part the spot's RNA explains and the rest (section r1;
smoothed for display, Gaussian, 2 spots; one colour scale, ratio units).
**g**, Each domain's mean short / long ratio minus the section's mean (ratio units), after depth (grey) and after
depth + RNA (rose), section r1.  Of every domain departing by more than 0.25, the RNA takes at least 70 % (white matter
-0.39 -> -0.06, choroid plexus -0.71 -> 0.00, cortex L2-4 +0.38 -> +0.06; CA1 keeps the most, +0.39 -> +0.12).
Cell-cycle scores explain nothing more (the ratio keeps 11.5 % and 13.8 % of its pattern after RNA, with or without
them).

Data: ``results/spatial_variable_brain{.json,_spots.parquet,_calibration.parquet}`` (U-Chrom
``benchmarks/spatial_variable --tissue brain``, stages through summary; provenance in the JSON).

    python figures/extended_data_fig3/make_ed_brain.py      # -> ed_brain.pdf / ed_brain.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.ndimage import gaussian_filter

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "benchmarks" / "spatial_variable" / "results"
MM = 1 / 25.4
W, H = 180.0, 145.5   # panel letters only (2026-10-14; 152 with titles)
SEC, SEC2 = "mouse_brain_r1", "mouse_brain_r2"

mpl.rcParams.update({
    "pdf.fonttype": 42, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "font.size": 6, "axes.titlesize": 6, "axes.labelsize": 6,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "xtick.minor.size": 1, "ytick.minor.size": 1,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.0,
    "axes.labelpad": 1.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
})

INK, GREY, LGREY = "#2b2f36", "#8a8f98", "#c9ccd1"
INDIGO, ROSE = "#34568b", "#b5507a"
RNA_C, LEFT_C = "#8fa8c8", "#eceef1"
#: the authors' domains (names of the benchmark's TISSUES["brain"]["domains"]) in the order of panels c and g
DOMAINS = [("cortex L2-4", "#b5507a"), ("cortex L5-6", "#d98cb0"), ("CA1", "#e07a3f"), ("CA3", "#f2ab6d"),
           ("DG", "#d9a441"), ("neuropil", "#cdbbe0"), ("striatum", "#8c6bb1"), ("thalamus 1", "#6f9a3a"),
           ("thalamus 2", "#a9c76d"), ("habenula", "#3f6b2f"), ("white matter", "#34568b"),
           ("choroid plexus", "#5aa0d8")]
#: the features of panels b and d (U-Chrom names) and their labels, as in Fig. 2
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


S = json.loads(need(RESULTS / "spatial_variable_brain.json").read_text())
S = S.get("results", S)
SPOTS = pd.read_parquet(need(RESULTS / "spatial_variable_brain_spots.parquet"))
CAL = pd.read_parquet(need(RESULTS / "spatial_variable_brain_calibration.parquet"))
R1 = SPOTS[SPOTS.section == SEC]
if R1["domain"].isna().all():
    raise SystemExit("no domains in the spots table: rerun the benchmark's summary stage (U-Chrom, --tissue brain)")


def fmt_p(p: float) -> str:
    if p <= 0:
        return "P < 1e-300"
    e = int(np.floor(np.log10(p)))
    return f"P = {p / 10 ** e:.1f}e{e}" if e < -2 else f"P = {p:.2g}"


def grid(t: pd.DataFrame, col) -> np.ndarray:
    """Spots on the array grid, cropped to the section."""
    x, y = t["x"].to_numpy(int), t["y"].to_numpy(int)
    x, y = x - x.min(), y - y.min()
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


def explained_by_rna(sec: str, feat: str) -> float:
    """The share of a feature's spatial pattern the spot's RNA explains (stage cnv for the contact features, stage
    local for the pooled window scores: the smooth variance after depth + RNA over that after depth)."""
    if feat.startswith("window:"):
        return 1 - S["local_per_spot"][f"{sec}/{feat.split(':')[1]}"]["depth_rna30.smooth_left"]
    return 1 - S["cnv"][sec]["smooth_left"]["rna30"][f"q30.{feat}"]


def column(feat: str) -> str:
    """The spots table's column of a feature (destriped)."""
    return f"window.{feat.split(':')[1]}" if feat.startswith("window:") else f"destriped.q30.{feat}"


def depth_removed(t: pd.DataFrame, col: str) -> np.ndarray:
    """The destriped feature minus its least-squares fit on log depth (+ intercept): the start of the RNA step."""
    v = t[f"destriped.{col}"].to_numpy(float)
    X = np.c_[np.ones(len(t)), t["log_contacts"].to_numpy(float)]
    ok = np.isfinite(v)
    beta = np.linalg.lstsq(X[ok], v[ok], rcond=None)[0]
    return v - X @ beta


# --------------------------------------------------------------------------- a, b, c

def panel_a(x0, y0, s):
    head(x0, y0, "a")
    y0 -= 2.0  # no title line: the content moves up
    code = {name: i for i, (name, _) in enumerate(DOMAINS)}
    img = grid(R1, [code.get(d, np.nan) for d in R1["domain"]])
    cmap = mpl.colors.ListedColormap([c for _, c in DOMAINS])
    ax = ax_mm(x0 + 0.5, y0 + 5.5, s, s)
    ax.imshow(np.ma.masked_invalid(img), cmap=cmap, vmin=-0.5, vmax=len(DOMAINS) - 0.5, interpolation="none")
    ax.axis("off")
    for k, (lab, c) in enumerate(DOMAINS):
        cx, cy = x0 + 1.0 + (k % 2) * (s / 2 + 1.5), y0 + 5.5 + s + 2.4 + (k // 2) * 2.6
        A.add_patch(Rectangle((cx, cy - 0.85), 1.7, 1.7, color=c, lw=0))
        T(cx + 2.4, cy, lab, 4.8)
    return y0 + 5.5 + s + 2.4 + 6 * 2.6


def panel_b(x0, y0, w):
    head(x0, y0, "b")
    y0 -= 2.0  # no title line: the content moves up
    n, gap = 4, 1.0
    s = (w - (n - 1) * gap) / n
    for k, (f, lab) in enumerate(FEATURES):
        r, c = divmod(k, n)
        x, y = x0 + c * (s + gap), y0 + 5.5 + r * (s + 4.6)
        show(ax_mm(x, y, s, s), grid(R1, column(f)))
        T(x + s / 2, y + s + 1.8, lab, 4.6, INK, ha="center")
    yk = y0 + 5.5 + 2 * (s + 4.6) + 0.4                              # one bar for all: each map is scaled to its
    xm, bw_ = x0 + w / 2, s * 1.2                                    # own 2nd-98th percentile (show)
    cax = ax_mm(xm - bw_ / 2, yk, bw_, 1.2)
    mpl.colorbar.ColorbarBase(cax, cmap=plt.get_cmap("viridis"), orientation="horizontal", ticks=[])
    cax.set_frame_on(False)
    T(xm - bw_ / 2 - 1.0, yk + 0.6, "low", 4.6, INK, ha="right")
    T(xm + bw_ / 2 + 1.0, yk + 0.6, "high", 4.6, INK)
    T(xm + bw_ / 2 + 7.0, yk + 0.6, "each feature on its own scale", 4.4, GREY)
    return yk + 1.2


def panel_c(x0, y0, w, h):
    head(x0, y0, "c")
    y0, h = y0 - 2.0, h + 2.0  # no title line: the content moves up
    ax = ax_mm(x0 + 16.5, y0 + 6.5, w - 22.0, h - 13.5)
    v = R1["destriped.q30.short_long_ratio"]
    for i, (name, col) in enumerate(DOMAINS):
        d = v[R1.domain == name].dropna().to_numpy()
        q5, q25, q50, q75, q95 = np.quantile(d, [0.05, 0.25, 0.5, 0.75, 0.95])
        ax.add_patch(Rectangle((q25, i - 0.32), q75 - q25, 0.64, color=col, lw=0))
        ax.plot([q50, q50], [i - 0.32, i + 0.32], color="white", lw=0.8)
        ax.plot([q5, q25], [i, i], color=col, lw=0.6)
        ax.plot([q75, q95], [i, i], color=col, lw=0.6)
    ax.set_ylim(len(DOMAINS) - 0.5, -0.5)
    ax.set_yticks(range(len(DOMAINS)), [n for n, _ in DOMAINS], fontsize=4.8)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("short / long ratio")
    lo, hi = ax.get_xlim()
    xr = hi + 0.03 * (hi - lo)
    for a_, b_, s in ((0, 4, "high"), (10, 11, "low")):        # the contrast of Guo et al.'s statement
        ax.plot([xr, xr], [a_ - 0.3, b_ + 0.3], color=INK, lw=0.5, clip_on=False)
        ax.text(xr + 0.03 * (hi - lo), (a_ + b_) / 2, s, fontsize=4.6, rotation=90, ha="left", va="center",
                color=INK, clip_on=False)
    ax.set_xlim(lo, hi)
    p = S[f"domain_high_vs_low/q30/{SEC}"]["mannwhitney_p_two_sided"]
    T(x0 + 16.5 + (w - 22.0) / 2, y0 + h - 0.6, f"high vs low: {fmt_p(p)}", 4.8, INK, ha="center")


# --------------------------------------------------------------------------- d, e

def panel_d(x0, y0, w, h):
    head(x0, y0, "d")
    y0, h = y0 - 2.0, h + 2.0  # no title line: the content moves up
    ax = ax_mm(x0 + 20.0, y0 + 6.0, w - 21.0, h - 19.0)
    for k, (f, lab) in enumerate(FEATURES):
        for j, sec in enumerate((SEC, SEC2)):
            a_ = explained_by_rna(sec, f)
            yy = k * 1.0 + j * 0.42
            ax.barh(yy, a_, left=0, height=0.38, color=RNA_C, lw=0)
            ax.barh(yy, 1 - a_, left=a_, height=0.38, color=LEFT_C, lw=0)
    ax.set_yticks([k + 0.21 for k in range(len(FEATURES))], [lab for _, lab in FEATURES], fontsize=4.8)
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1); ax.set_xticks([0, 0.5, 1], ["0", "50", "100 %"])
    ax.tick_params(axis="x", labelsize=4.8)
    ax.set_xlabel("share of the spatial pattern", fontsize=5.0)
    hs = [Rectangle((0, 0), 1, 1, color=c) for c in (RNA_C, LEFT_C)]
    ax.legend(hs, ["explained by RNA", "unexplained"], loc="upper center", ncol=2, bbox_to_anchor=(0.4, -0.25),
              frameon=False, fontsize=4.8, handlelength=0.9, handletextpad=0.3, columnspacing=0.9)


def panel_e(x0, y0, w, h):
    head(x0, y0, "e")
    y0, h = y0 - 2.0, h + 2.0  # no title line: the content moves up
    names = {f"q30.{f}": lab for f, lab in FEATURES}
    gap = 4.0
    aw = (w - 10.0 - gap) / 2
    for j, sec in enumerate((SEC, SEC2)):
        ax = ax_mm(x0 + 9.0 + j * (aw + gap), y0 + 7.5, aw, h - 16.5)
        t = CAL[CAL.section == sec]
        g = t[t.kind == "gene"]
        m = t[t.kind == "module"]
        f = t[(t.kind == "hic") & t.feature.isin(names)]
        ax.axvline(0.2, color=GREY, lw=0.5, ls=(0, (2, 2)))
        ax.scatter(g.smooth_share_depth, g.smooth_left.clip(upper=1.05), s=1.0, color=LGREY, lw=0, rasterized=True)
        ax.scatter(m.smooth_share_depth, m.smooth_left.clip(upper=1.05), s=9, facecolor="none", edgecolor=INK, lw=0.5)
        ax.scatter(f.smooth_share_depth, f.smooth_left, s=9, color=ROSE, lw=0, zorder=3)
        for r in f.itertuples():
            if r.feature in ("q30.short_long_ratio", "q30.ps_slope_100k_1M"):
                up = r.feature == "q30.short_long_ratio"
                ax.annotate(names[r.feature], (r.smooth_share_depth, r.smooth_left), xytext=(0, 7) if up else (3, 2),
                            textcoords="offset points", fontsize=4.4, color=ROSE, ha="center" if up else "left")
        ax.set_xlim(0, 0.8); ax.set_ylim(0, 1.08)
        ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8], ["0", "20", "40", "60", "80"])
        ax.set_yticks([0, 0.5, 1], ["0", "50", "100"])
        ax.set_title(f"section {'r1' if j == 0 else 'r2'}", fontsize=5.5, pad=1.5)
        ax.set_xlabel("variance in the spatial pattern (%)", fontsize=5.2)
        if j == 0:
            ax.set_ylabel("pattern left after RNA (%)", fontsize=5.2)
            hs = [mpl.lines.Line2D([], [], ls="", marker="o", ms=2.6, color=ROSE, mew=0),
                  mpl.lines.Line2D([], [], ls="", marker="o", ms=1.6, color=GREY, mew=0),
                  mpl.lines.Line2D([], [], ls="", marker="o", ms=2.8, mfc="none", mec=INK, mew=0.5)]
            ax.legend(hs, ["contact features", "held-out genes", "gene modules"], loc="upper right", frameon=False,
                      fontsize=4.8, handletextpad=0.2, borderaxespad=0.2)
        else:
            ax.set_yticklabels([])


# --------------------------------------------------------------------------- f, g

def panel_f(x0, y0):
    head(x0, y0, "f")
    y0 -= 3.0  # no title line: the content moves up
    col = "q30.short_long_ratio"
    after_depth = smoothed(grid(R1, depth_removed(R1, col)))
    after_rna = smoothed(grid(R1, f"rna30.{col}"))
    maps = [(after_depth, "after depth"), (after_depth - after_rna, "explained by RNA"), (after_rna, "unexplained")]
    lim = np.nanquantile(np.abs(after_depth), 0.98)
    w, x = 25.0, x0 + 1.0
    for i, (img, title) in enumerate(maps):
        ax = ax_mm(x, y0 + 8.5, w, w)
        show(ax, img, "RdBu_r", lim=(-lim, lim))
        T(x + w / 2, y0 + 7.0, title, 5.0, INK, ha="center")
        if i < 2:
            T(x + w + 1.0, y0 + 8.5 + w / 2, "=" if i == 0 else "+", 7, GREY, ha="center")
        x += w + 2.0
    cax = ax_mm(x0 + 1.0 + 3.0, y0 + 8.5 + w + 1.8, w - 6.0, 1.3)
    mpl.colorbar.ColorbarBase(cax, cmap=plt.get_cmap("RdBu_r"), norm=mpl.colors.Normalize(-lim, lim),
                              orientation="horizontal", ticks=[-lim, 0, lim], format="%.2f")
    cax.tick_params(length=1.5, pad=1, labelsize=4.6)
    T(x0 + 1.0 + w - 3.0 + 1.5, y0 + 8.5 + w + 2.45, "short / long ratio after depth", 4.6, INK)
    return y0 + 8.5 + w + 6.0


def panel_g(x0, y0, w, h):
    head(x0, y0, "g")
    y0, h = y0 - 2.0, h + 2.0  # no title line: the content moves up
    col = "q30.short_long_ratio"
    t = R1.assign(after_depth=depth_removed(R1, col), after_rna=R1[f"rna30.{col}"])
    t = t[t[["after_depth", "after_rna"]].notna().all(axis=1)]
    t[["after_depth", "after_rna"]] -= t[["after_depth", "after_rna"]].mean()     # the section mean = 0
    means = t.groupby("domain")[["after_depth", "after_rna"]].mean()
    ax = ax_mm(x0 + 10.0, y0 + 7.0, w - 10.5, h - 19.0)
    bw = 0.38
    for i, (name, _) in enumerate(DOMAINS):
        ax.bar(i - bw / 2, means.loc[name, "after_depth"], width=bw * 0.94, color=GREY, lw=0)
        ax.bar(i + bw / 2, means.loc[name, "after_rna"], width=bw * 0.94, color=ROSE, lw=0)
    ax.axhline(0, color=INK, lw=0.5)
    ax.set_xlim(-0.6, len(DOMAINS) - 0.4)
    ax.set_xticks(range(len(DOMAINS)), [n for n, _ in DOMAINS], fontsize=4.8, rotation=40, ha="right",
                  rotation_mode="anchor")
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("short / long ratio,\ndomain \u2212 section mean", fontsize=5.2)   # (no taller than its axes)
    hs = [Rectangle((0, 0), 1, 1, color=c) for c in (GREY, ROSE)]
    ax.legend(hs, ["after depth", "after depth + RNA"], loc="lower left", frameon=False, fontsize=4.8,
              handlelength=0.9, handletextpad=0.3)


def make_figure() -> Path:
    sa = 33.0
    a_end = panel_a(0.5, 0.0, sa)
    xb, wb = 40.0, 82.0
    b_end = panel_b(xb, 0.0, wb)
    xc = xb + wb + 3.5
    panel_c(xc, 0.0, W - xc, a_end)
    y = max(a_end, b_end) + 3.0
    hd = 50.0
    panel_d(0.5, y, 66.0, hd)
    panel_e(70.0, y, W - 70.0, hd)
    y += hd + 1.0
    f_end = panel_f(0.5, y)
    panel_g(88.0, y, W - 88.0, f_end - y)
    print(f"figure height used: {f_end:.1f} of {H} mm")
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"ed_brain.{ext}", dpi=300)
    return HERE / "ed_brain.pdf"


if __name__ == "__main__":
    print(make_figure())
