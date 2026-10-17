"""Brief Communication, Supplementary Fig. 2 -- per-cell contact maps linked, not copied, and analysed with the RNA.

The data and the drawing of the full manuscript's Supplementary Fig. S1 (``../../full/figures/fig2/make_figS_fig2.py``)
c-e, without panel titles: what the linked maps are used for.  Its access-by-format panels a, b and the cost of
the linked maps (f-h) are Supplementary Fig. 1 (``make_supp_access.py``):

  a  RNA UMAP: three Leiden clusters computed in ChromData        (full: c)
  b  pseudo-bulk contact maps of two clusters from the linked .scool (full: d)
  c  log2 ratio of the two maps, at the size of the maps of b      (full: e)

Grid (mm, 180 wide, one row): 30-mm axes at x 12 (the UMAP), 54 and 88 (the maps of b, colour bar after them) and
139 (c, its colour bar after it); letters at x 0, 44 and 124.  The numbers come from the full manuscript's results
(``../../full/figures/fig2/results``), as there.

    python manuscripts/nmeth-bc/figures/make_supp_linked.py      # -> supp_linked.pdf / supp_linked.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea

HERE = Path(__file__).resolve().parent
FULL = HERE.parents[1] / "full" / "figures" / "fig2"
sys.path.insert(0, str(FULL))
import make_figS_fig2 as S  # noqa: E402

F = S.F
W, H = 180.0, 46.0

# the grid (mm)
SIDE = 30.0                             # every axes SIDE wide
UMAP_W = SIDE
COL = [12.0, 54.0, 88.0, 139.0]         # axes x: the UMAP, the two maps of b, the ratio map
LX = [0.0, 44.0, 124.0]                 # letters: a, b, c
ROWS = [(1.0, 8.0, 30.0)]               # (letter y, axes top, axes height)


def ax_at(fig, x, y, w, h):
    return F.ax_at(fig, x, y, w, h)


def letter(fig, x, row, s):
    F.label(fig, x, ROWS[row][0], s)


def _mb_ticks(mx, w, n):
    mx.set_xticks([-0.5, n - 0.5], [f"{w['start'] // 10**6}", f"{w['end'] // 10**6} Mb"])
    mx.set_yticks([])
    mx.tick_params(axis="x", length=1.5, pad=1)
    mx.get_xticklabels()[0].set_ha("left"); mx.get_xticklabels()[1].set_ha("right")


def _frame(mx, col):
    for sp in mx.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.6); sp.set_color(col)


def _cbar(fig, x, top, h, norm, cmap):
    cax = ax_at(fig, x, top + 2, 1.6, h - 4)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=S.TICK, width=0.4, length=1.5, pad=1)
    return cb


def umap(fig, fx):
    letter(fig, LX[0], 0, "a")
    _, top, h = ROWS[0]
    ax = ax_at(fig, COL[0], top, UMAP_W, h)
    groups, um = fx["groups"], fx["umap"]
    xs, ys, cl = np.asarray(um["x"]), np.asarray(um["y"]), np.asarray(um["cluster"])
    ax.scatter(xs, ys, s=0.4, color=S.OTHER, lw=0, rasterized=True)
    for r in S.ROLES:
        m = cl == groups[r]["cluster"]
        ax.scatter(xs[m], ys[m], s=0.6, color=S.ROLE_COL[r], lw=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel("RNA UMAP 1", labelpad=1); ax.set_ylabel("RNA UMAP 2", labelpad=1)
    span = ys.max() - ys.min()
    for r in S.ROLES:                    # astrocytes and neurons above their clusters; oligodendrocytes bottom left,
        m = cl == groups[r]["cluster"]   # where the UMAP leaves room (below its cluster it would run into the neurons)
        if r == "oligodendrocyte":
            ax.text(0.0, 0.22, S.NAMES[r], transform=ax.transAxes, fontsize=6, color=S.ROLE_COL[r], ha="left",
                    va="top")
        else:
            ax.text(np.median(xs[m]), np.percentile(ys[m], 99) + 0.03 * span, S.NAMES[r], fontsize=6,
                    color=S.ROLE_COL[r], ha="center", va="bottom")


def maps(fig, fx):
    letter(fig, LX[1], 0, "b")
    _, top, h = ROWS[0]
    groups = fx["groups"]
    w = groups["neuronal"]["window"]
    cmap = LinearSegmentedColormap.from_list("wt", ["#ffffff", "#b9d9f0", F.INDIGO, "#101d33"])
    pair = ["neuronal", "oligodendrocyte"]
    mats = [np.asarray(groups[r]["window"]["matrix"], float) / groups[r]["total_cis_contacts"] * 1e6 for r in pair]
    pos = np.concatenate([m[m > 0] for m in mats])
    norm = LogNorm(vmin=np.percentile(pos, 5), vmax=np.percentile(pos, 99.5))
    for k, (r, mat) in enumerate(zip(pair, mats)):
        mx = ax_at(fig, COL[1 + k], top, SIDE, h)
        mx.imshow(mat, cmap=cmap, norm=norm, interpolation="nearest")
        _mb_ticks(mx, w, mat.shape[0])
        _frame(mx, S.ROLE_COL[r])
        mx.set_title(S.NAMES[r], fontsize=6.5, color=S.ROLE_COL[r], pad=2)
    cb = _cbar(fig, COL[2] + SIDE + 2.5, top, h, norm, cmap)
    lo, hi = np.ceil(np.log10(norm.vmin)), np.floor(np.log10(norm.vmax))
    S.ticks(cb.ax, "y", [10.0 ** k for k in range(int(lo), int(hi) + 1)], "plain")
    cb.set_label("contacts per million cis", fontsize=6, color=F.GREY, labelpad=2)


def ratio(fig, fx):
    letter(fig, LX[2], 0, "c")
    _, top, h = ROWS[0]
    rr = fx["window_log2_ratio"]
    w = fx["groups"]["neuronal"]["window"]
    num, den = rr["numerator"], rr["denominator"]
    mat = np.asarray(rr["matrix"], float)
    lim = float(np.nanpercentile(np.abs(mat[np.isfinite(mat)]), 98))
    norm = mpl.colors.TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
    cmap = LinearSegmentedColormap.from_list("div", [S.ROLE_COL[den], "#ffffff", S.ROLE_COL[num]])
    mx = ax_at(fig, COL[3], top, SIDE, h)
    mx.imshow(mat, cmap=cmap, norm=norm, interpolation="nearest")
    _mb_ticks(mx, w, mat.shape[0])
    _frame(mx, F.GREY)
    # the title names the two groups in their colours, which are the two ends of the colour scale
    parts = [TextArea(S.NAMES[num], textprops=dict(color=S.ROLE_COL[num], fontsize=6.5)),
             TextArea(" / ", textprops=dict(color=F.INK, fontsize=6.5)),
             TextArea(S.NAMES[den], textprops=dict(color=S.ROLE_COL[den], fontsize=6.5))]
    mx.add_artist(AnchoredOffsetbox(loc="lower center", child=HPacker(children=parts, align="baseline", pad=0, sep=0),
                                    pad=0, borderpad=0, frameon=False, bbox_to_anchor=(0.5, 1.0),
                                    bbox_transform=mx.transAxes))
    cx = COL[3] + SIDE + 2.5
    cb = _cbar(fig, cx, top, h, norm, cmap)
    t = np.floor(lim * 2) / 2
    cb.set_ticks([-t, 0, t], labels=[f"{-t:g}".replace("-", "−"), "0", f"+{t:g}"])
    fig.text((cx + 0.8) / W, 1 - (top + 1.2) / H, r"$\log_2$", fontsize=6, color=F.GREY, ha="center", va="bottom")


def main() -> None:
    S.W, S.H = W, H
    F.W, F.H = W, H
    S.style()
    fx = F._load("f_crossmodal")
    fig = plt.figure(figsize=(W * F.MM, H * F.MM))
    umap(fig, fx)
    maps(fig, fx)
    ratio(fig, fx)
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"supp_linked.{ext}", dpi=300)
    print(HERE / "supp_linked.pdf", f"{W:.0f} x {H:.0f} mm")


if __name__ == "__main__":
    main()
