"""Brief Communication, Supplementary Fig. 1 -- reading the store, and the contact maps linked to it.

Row 1: the access-by-format panels a, b of the full manuscript's Supplementary Fig. S1
(``../../full/figures/fig2/make_figS_fig2.py``, ``by_format``), drawn by that code from the same results
(``../../full/figures/fig2/results/c_access.json``), without panel titles.  Rows 1-2, c-e: the cost of a linked
per-cell contact map (``e_query.json``, ``e_contacts.json``; the full manuscript's S1 f-h), moved here from
Supplementary Fig. 2, which keeps what the linked maps are used for:

  a  open or load the file and get one trace, by format and size
  b  the same for one chromosome (chromdata.zarr with all per-spot tracks, solid, or coordinates only, dashed)
  c  time to attach the scHiCAR per-cell maps: link vs copy
  d  pseudo-bulk query time, first vs cached: each cell type over the genome
  e  the same for random groups on chr1
  f  on-disk size, linked vs copied (decimal MB)

Grid (mm, 180 x 100): three 60-mm columns and two rows, every axes 42 x 35 mm, 13 mm right of its letter; a and b
name their series in a key at the top left; no titles in the panels (their subjects are in the legend).

    python manuscripts/nmeth-bc/figures/make_supp_access.py      # -> supp_access.pdf / supp_access.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
FULL = HERE.parents[1] / "full" / "figures" / "fig2"
sys.path.insert(0, str(FULL))
import make_figS_fig2 as S  # noqa: E402

F = S.F
W, H = 180.0, 100.0
# a 3 x 2 grid: columns 60 mm wide, every axes AX_W x AX_H (by_format's), the letter 13 mm left of its axes
AX_W, AX_H = 42.0, S.ROW1_H
COLS = (0.0, 60.0, 120.0)
ROWS = (1.0, 51.0)                       # letter y; the axes from ROW + S.AX_DY


def _ax(fig, x0, y0):
    return F.ax_at(fig, x0 + 13, y0 + S.AX_DY, AX_W, AX_H)


def link_copy(fig, x0, y0, letter, imp):
    F.label(fig, x0, y0, letter)
    ax = _ax(fig, x0, y0)
    t = [imp["link"]["seconds"], imp["copy"]["seconds"]]
    med = np.array([v["median"] for v in t])
    err = np.array([[v["median"] - v["min"] for v in t], [v["max"] - v["median"] for v in t]])
    ax.bar([0, 1], med, width=0.55, color=[F.TEAL, F.LGREY], edgecolor=[F.TEAL, F.GREY], lw=0.5, yerr=err,
           error_kw=dict(lw=0.6, capsize=1.5, capthick=0.6, ecolor=F.INK))
    for i, v in enumerate(med):
        ax.text(i, v * 1.35, f"{v:.2f} s" if v < 1 else f"{v:.1f} s", ha="center", va="bottom", fontsize=6)
    ax.set_yscale("log")
    ax.set_ylim(0.1, 100)
    S.ticks(ax, "y", [0.1, 1, 10, 100], "time")
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0, 1], ["link", "copy"])
    ax.tick_params(axis="x", length=0, pad=1.5)
    ax.set_ylabel("time")


def query(fig, x0, y0, letter, xs, first, cached, xlabel, key=False):
    """First and repeated (cached) query time against the cells summed."""
    F.label(fig, x0, y0, letter)
    ax = _ax(fig, x0, y0)
    ax.plot(xs, first, color=F.TEAL, marker="o", ms=2.0, lw=1.1)
    ax.plot(xs, cached, color=F.TEAL, ls="--", marker="o", ms=2.0, mfc="white", lw=1.0)
    F.loglog(ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("time")
    ax.set_xlim(xs.min() / 1.6, xs.max() * 1.6)
    S.ticks(ax, "x", [v for v in (1e2, 1e3) if xs.min() / 1.6 < v < xs.max() * 1.6], "count")
    ax.set_ylim(1e-5, 30)
    S.ticks(ax, "y", [1e-4, 1e-2, 1e0], "time")
    if key:
        hs = [Line2D([], [], color=F.TEAL, marker="o", ms=2.4, lw=1.1),
              Line2D([], [], color=F.TEAL, ls="--", marker="o", ms=2.4, mfc="white", lw=1.0)]
        ax.legend(hs, ["first query", "cached"], loc="center left", frameon=False, fontsize=6,
                  handlelength=2.0, handletextpad=0.4, labelspacing=0.3, borderaxespad=0.2)


def queries(fig, e, q):
    per = q["identity"]["per_type"]
    n_chr = len(q["identity"]["scope"].split("(")[1].split(")")[0].split(","))
    n = np.array([p["n_cells"] for p in per]); o = np.argsort(n)
    first = np.array([p["first_query_s_total"] for p in per])
    cached = np.array([p["cached_query_s_median"] for p in per]) * n_chr
    query(fig, COLS[0], ROWS[1], "d", n[o], first[o], cached[o], "cells of the type", key=True)
    c = np.array([p["n_cells"] for p in e])
    query(fig, COLS[1], ROWS[1], "e", c, np.array([p["cold_s"]["median"] for p in e]),
          np.array([p["cached_s"]["median"] for p in e]), "cells in the group")


def sizes(fig, x0, y0, letter, imp):
    F.label(fig, x0, y0, letter)
    ax = _ax(fig, x0, y0)
    mb = 1e-6                            # decimal megabytes, as every size on disk in the paper
    scool = imp["link"]["scool_bytes"] * mb
    store = {"linked": imp["link"]["store_bytes"] * mb, "copied": imp["copy"]["store_bytes"] * mb}
    x = np.array([0, 1, 2.5, 3.5])
    own = [store["linked"], store["copied"], store["linked"], store["copied"]]
    ext = [0, 0, scool, scool]
    ax.bar(x, own, width=0.7, color=[F.TEAL, F.LGREY] * 2, edgecolor=[F.TEAL, F.GREY] * 2, lw=0.5)
    ax.bar(x, ext, bottom=own, width=0.7, color="white", edgecolor=F.TEAL, hatch="////", lw=0.5)
    for xi, o, ex in zip(x, own, ext):
        ax.text(xi, o + ex + 6, f"{o + ex:.0f}", ha="center", va="bottom", fontsize=6, color=F.INK)
    ax.set_xticks(x, ["link", "copy", "link", "copy"])
    ax.tick_params(axis="x", length=0, pad=1.5)
    ax.set_ylabel("size (MB)")
    ax.set_xlim(-0.6, 4.1)
    ax.tick_params(axis="x", labelsize=5.5)
    top_mb = max(o + e for o, e in zip(own, ext)) * 1.18
    ax.set_ylim(0, top_mb)
    ax.set_yticks(np.arange(0, top_mb, 100))
    for xc, s in ((0.5, "store"), (3.0, "total")):
        ax.text(xc, -0.13, s, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6, color=F.GREY)
    ax.legend([Patch(fc="white", ec=F.TEAL, hatch="////", lw=0.5)], [".scool"], frameon=False,
              loc="upper left", handlelength=1.3, handletextpad=0.4, borderaxespad=0.1, fontsize=6)


def main() -> None:
    S.W, S.H = W, H
    F.W, F.H = W, H
    S.AX_DY = 3.0                                       # no title above the axes
    S.FORMAT_AX_W, S.FORMAT_KEY_INSIDE = AX_W, True     # a, b as wide as the others, their series in a key
    S.head = lambda fig, x0, y0, letter, title=None, sub=None: F.label(fig, x0, y0, letter)
    S.style()
    acc = F.access()
    q = F._load("e_query")
    e = F._load("e_contacts")["aggregation"]["points"]
    fig = plt.figure(figsize=(W * F.MM, H * F.MM))
    S.by_format(fig, COLS[0], ROWS[0], "a", "trace", None, None, acc)
    S.by_format(fig, COLS[1], ROWS[0], "b", "chrom", None, None, acc)
    link_copy(fig, COLS[2], ROWS[0], "c", q["import"])
    queries(fig, e, q)
    sizes(fig, COLS[2], ROWS[1], "f", q["import"])
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"supp_access.{ext}", dpi=300)
    print(HERE / "supp_access.pdf", f"{W:.0f} x {H:.0f} mm")


if __name__ == "__main__":
    main()
