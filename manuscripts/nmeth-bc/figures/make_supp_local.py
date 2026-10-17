"""Brief Communication, Supplementary Figure (local level) -- single boundaries and loops across a tissue: scored as
log ratios of their windows, they seem to vary in space the same way in both sections; scored against each spot's
own expected (U-Chrom ``uchrom.fea.window_expected``), they do not.

  a  each structure's smooth share of its per-spot score above that of shuffled spots (the excess), section r1
     against r2, after depth + RNA + copy number: log ratios (left) and observed - expected shares (right);
     boundaries coloured by chromosome X (haploid in these male tissues), loops by their span
  b  the agreement between the sections (Spearman of the excess), per tissue, kind and score
  c  pooled over all structures, a spot's insulation and loop strength against its expected do vary in space:
     smooth share after depth + RNA + copy number vs shuffled spots, both sections

Data: ``results/spatial_variable_{brain,melanoma}_local.parquet`` and ``.json`` (U-Chrom
``benchmarks/spatial_variable``, stage local + summary).

    python manuscripts/nmeth-bc/figures/make_supp_local.py      # -> supp_local.pdf / supp_local.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
MM = 1 / 25.4
W, H = 180.0, 108.0
COL = "depth_rna30_cnv.smooth_share"

mpl.rcParams.update({
    "pdf.fonttype": 42, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "font.size": 6, "axes.titlesize": 6, "axes.labelsize": 6,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelpad": 1.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
})
GREY, LGREY, INK, ROSE, TEAL, INDIGO = "#8a8f98", "#c9ccd1", "#2b2f36", "#b5507a", "#2a9d8f", "#34568b"
TISSUES = [("brain", "mouse brain"), ("melanoma", "human melanoma")]
VERSIONS = [("logratio", "log ratio of the windows"), ("oe", "observed - expected share")]
PLURAL = {"boundary": "boundaries", "loop": "loops"}

fig = plt.figure(figsize=(W * MM, H * MM))
A = fig.add_axes([0, 0, 1, 1])
A.set_xlim(0, W); A.set_ylim(H, 0); A.axis("off")


def T(x, y, s, size=6, color=INK, ha="left", va="center", weight="normal"):
    return A.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, fontweight=weight)


def head(x, y, letter, s, sub=None):
    T(x, y, letter, 8, weight="bold", va="top")
    T(x + 4.5, y + 0.4, s, 7, weight="bold", va="top")
    if sub:
        T(x + 4.5, y + 3.9, sub, 5.5, GREY, va="top")


def ax_mm(x, y, w, h):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def local(tissue: str) -> pd.DataFrame:
    return pd.read_parquet(RESULTS / f"spatial_variable_{tissue}_local.parquet")


def excess(t: pd.DataFrame, kind: str, version: str) -> pd.DataFrame:
    """Per structure in both sections: smooth share minus that of shuffled spots, with its chromosome and span."""
    d = t[(t.kind == kind) & (t.version == version)]
    secs = sorted(d.section.unique())
    a, b = (d[d.section == s].set_index("feature") for s in secs)
    common = a.index.intersection(b.index)
    return pd.DataFrame({"r1": (a[COL] - a[COL + "_null"]).loc[common], "r2": (b[COL] - b[COL + "_null"]).loc[common],
                         "chrom": a.loc[common, "chrom"], "span": a.loc[common, "span"]})


def panel_a(y0):
    head(0.5, y0, "a", "Single structures: section r1 against r2",
         "excess smooth share of each structure's per-spot score (over shuffled spots), after depth, RNA and "
         "copy number; human melanoma")
    t = local("melanoma")
    w, h, x = 38.0, 34.0, 10.0
    for kind in ("boundary", "loop"):
        for version, label in VERSIONS:
            e = excess(t, kind, version)
            ax = ax_mm(x, y0 + 12.0, w, h)
            if kind == "boundary":
                x_chr = e.chrom.isin(["chrX"])
                ax.scatter(e.r1[~x_chr], e.r2[~x_chr], s=1.2, color=LGREY, lw=0, rasterized=True)
                ax.scatter(e.r1[x_chr], e.r2[x_chr], s=2.0, color=ROSE, lw=0, rasterized=True)
            else:
                span = np.log10(e.span)
                lo, hi = np.quantile(span, [0.05, 0.95])
                ax.scatter(e.r1, e.r2, s=1.5, c=span, cmap="viridis_r", vmin=lo, vmax=hi, lw=0, rasterized=True)
            r = stats.spearmanr(e.r1, e.r2)[0]
            lim = np.nanquantile(np.r_[e.r1, e.r2], [0.002, 0.998])
            ax.set_xlim(*lim); ax.set_ylim(*lim)
            ax.axhline(0, color=GREY, lw=0.4); ax.axvline(0, color=GREY, lw=0.4)
            ax.set_title(f"{PLURAL[kind]}, {label}\nSpearman {r:.2f}", fontsize=5.5, pad=2)
            ax.set_xlabel("excess, r1")
            if x < 20:
                ax.set_ylabel("excess, r2")
            x += w + 5.0
    T(10.0, y0 + 12.0 + h + 8.5, "boundaries: chromosome X in rose (haploid: half the contacts); loops: colour = span "
      "(yellow: short, 5-95 % of the spans)", 5.0, GREY)
    return y0 + 12.0 + h + 11.0


def panel_b(y0):
    head(0.5, y0, "b", "Agreement between the sections",
         "Spearman correlation of the structures' excess, r1 vs r2, after depth + RNA + copy number")
    ax = ax_mm(12.0, y0 + 11.0, 70.0, 30.0)
    rows, x = [], 0
    for tissue, tl in TISSUES:
        t = local(tissue)
        for kind in ("boundary", "loop"):
            for i, (version, label) in enumerate(VERSIONS):
                e = excess(t, kind, version)
                r = stats.spearmanr(e.r1, e.r2)[0]
                ax.bar(x + 0.42 * i, r, width=0.4, color=[LGREY, TEAL][i], lw=0)
            rows.append((x + 0.21, f"{tl.split()[-1]}\n{PLURAL[kind]}"))
            x += 1.2
    ax.set_xticks([p for p, _ in rows], [s for _, s in rows], linespacing=1.0)
    ax.tick_params(axis="x", length=0)
    ax.axhline(0, color=INK, lw=0.5)
    ax.set_ylabel("Spearman, r1 vs r2")
    ax.legend([mpl.patches.Patch(color=LGREY), mpl.patches.Patch(color=TEAL)], [l for _, l in VERSIONS],
              loc="upper left", frameon=False, fontsize=5.0)
    return ax


def panel_c(x0, y0):
    head(x0, y0, "c", "Pooled over all structures",
         "a spot's insulation and loop strength against its expected: smooth share vs shuffled spots")
    ax = ax_mm(x0 + 12.0, y0 + 11.0, 70.0, 30.0)
    x, ticks = 0, []
    for tissue, tl in TISSUES:
        d = json.loads((RESULTS / f"spatial_variable_{tissue}.json").read_text())
        d = d.get("results", d)["local_per_spot"]
        for feat, lab in (("boundary_strength_oe", "insulation"), ("loop_strength_oe", "loops")):
            for j, sec in enumerate(sorted({k.split("/")[0] for k in d})):
                v = d[f"{sec}/{feat}"]
                ax.bar(x + 0.42 * j, 100 * v[COL], width=0.4, color=[ROSE, "#d9a0b8"][j], lw=0)
                ax.plot([x + 0.42 * j - 0.2, x + 0.42 * j + 0.2], [100 * v[COL + "_null"]] * 2, color=INK, lw=0.8)
            ticks.append((x + 0.21, f"{tl.split()[-1]}\n{lab}"))
            x += 1.2
    ax.set_xticks([p for p, _ in ticks], [s for _, s in ticks], linespacing=1.0)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("smooth share (%)")
    ax.legend([mpl.patches.Patch(color=ROSE), mpl.patches.Patch(color="#d9a0b8"), mpl.lines.Line2D([], [], color=INK)],
              ["r1", "r2", "shuffled spots"], loc="upper right", frameon=False, fontsize=5.0)
    return ax


def make_figure() -> Path:
    y = panel_a(0.0)
    panel_b(y + 1.0)
    panel_c(92.0, y + 1.0)
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"supp_local.{ext}", dpi=300)
    return HERE / "supp_local.pdf"


if __name__ == "__main__":
    print(make_figure())
