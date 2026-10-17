"""Brief Communication, Supplementary Fig. 5 -- the browser's response times and the atlas.  (The browser on each kind
of data, Supplementary Figs 2-4, are whole-page captures: ``capture_supp_browser.py``; ``panel_browser`` below drew
the former composite of single views and is kept for reference.)

Supplementary Fig. 5 (``supp_atlas``):
  a  response times, locally (full Fig. 6c: ``../../full/figures/fig6/results/frontend.json``, ``server.json``) and
     through the hosted browser reading the atlas stores over HTTP (``results/hosted.json``)
  b  the atlas: 22 datasets by cells or spots and contacts or traced loci per cell (``results/atlas.json``)
  c  reprocessed vs published contacts per cell (``results/atlas.json``)

What a panel shows is said in its legend (``captions/supp_atlas.tex``), not on it.

    python manuscripts/nmeth-bc/figures/make_supp_browser_atlas.py      # -> supp_atlas.pdf / .png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from PIL import Image

HERE = Path(__file__).resolve().parent
CAPTURES = HERE / "captures"
RESULTS = HERE / "results"
FIG6 = HERE.parents[1] / "full" / "figures" / "fig6" / "results"
MM = 1 / 25.4
W, H = 180.0, 213.0                                   # H: set per figure by new_canvas

mpl.rcParams.update({
    "pdf.fonttype": 42, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "font.size": 6, "axes.titlesize": 6, "axes.labelsize": 6,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "xtick.major.size": 2, "ytick.major.size": 2, "xtick.minor.size": 1, "ytick.minor.size": 1,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.0,
    "axes.labelpad": 1.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
})

ORANGE, TEAL, PURPLE, INDIGO = "#e07a3f", "#2a9d8f", "#8e7cc3", "#34568b"
GREY, LGREY, INK = "#8a8f98", "#c9ccd1", "#2b2f36"
ROSE = "#b5507a"
#: atlas category -> colour of its kind of data (as in Fig. 1a)
KIND = {"Single-cell": TEAL, "Spatial": ROSE, "Imaging": ORANGE}
KIND_LABEL = {"Single-cell": "single-cell Hi-C", "Spatial": "spatial Hi-C", "Imaging": "chromatin tracing"}

fig = A = None


def new_canvas(h: float):
    """A figure W x h mm with an mm canvas (y down) behind the panels."""
    global fig, A, H
    H = h
    fig = plt.figure(figsize=(W * MM, H * MM))
    A = fig.add_axes([0, 0, 1, 1])
    A.set_xlim(0, W); A.set_ylim(H, 0); A.axis("off")


def T(x, y, s, size=6, color=INK, ha="left", va="center", weight="normal", **kw):
    return A.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, fontweight=weight, **kw)


def head(x, y, letter, s=None, sub=None):
    """The panel letter only: no panel titles (the legend says what each panel shows)."""
    T(x, y, letter, 8, weight="bold", va="top")


def ax_mm(x, y, w, h):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def load(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"{path} missing (see the docstring)")
    return json.loads(path.read_text())


# --------------------------------------------------------------------------- a

#: kind of data, name, dataset, the selection; then the three views (capture, caption)
SCENES = [
    ("Spatial", "Spatial Hi-C-RNA", "mouse brain", "L6 CT, 630 of 6,275 spots",
     [("spatial_tissue", "tissue"), ("spatial_umap", "RNA UMAP"), ("spatial_genome", "chr2 contacts, A/B")]),
    ("Single-cell", "Single-cell Hi-C + RNA", "HiRES, mouse brain", "Oli, 20 of 399 cells",
     [("scrna_umap", "RNA UMAP"), ("scrna_structure@1.6", "3D structure"), ("scrna_genome", "chr2 contacts")]),
    ("Imaging", "Chromatin tracing + IF", "DNA seqFISH+, cerebellum", "Bergmann glia, 192 of 1,799 cells",
     [("tracing_umap", "UMAP"), ("tracing_nuclei", "nuclei in 3D"), ("tracing_genome", "chr15 distances, IF")]),
]
A_X0, A_S, A_GAP = 41.0, 44.0, 0.8                     # the view grid of Supplementary Fig. 2


def capture(name: str) -> Image.Image:
    """A capture; ``name@z`` zooms in by z around the drawn content (the aspect ratio is kept)."""
    name, _, zoom = name.partition("@")
    im = Image.open(CAPTURES / f"{name}.png").convert("RGB")
    if zoom:
        a = np.asarray(im).astype(int)
        ink = np.argwhere((np.abs(a - a[0, 0]).sum(2) > 30))
        (r0, c0), (r1, c1) = ink.min(0), ink.max(0)
        cy, cx = (r0 + r1) / 2, (c0 + c1) / 2
        h, w = im.height / float(zoom), im.width / float(zoom)
        top = int(min(max(cy - h / 2, 0), im.height - h))
        left = int(min(max(cx - w / 2, 0), im.width - w))
        im = im.crop((left, top, left + int(w), top + int(h)))
    return im


def browser_height(y0: float = 0.0) -> float:
    """The height of Supplementary Fig. 2: its rows follow the captures' aspect ratios."""
    y = y0 + 1.0
    for *_, views in SCENES:
        im = capture(views[-1][0])
        y += A_S * im.height / im.width + 4.6
    return y


def panel_browser(y0=0.0):
    x0, s, gap = A_X0, A_S, A_GAP
    y = y0 + 1.0
    for kind, name, ds, focus, views in SCENES:
        col = KIND[kind]
        T(0.8, y + 3.0, name, 6.5, col, weight="bold")
        T(0.8, y + 6.4, ds, 5.5, INK)
        T(0.8, y + 9.6, focus, 5.5, GREY)
        hh = 0.0
        for c, (cap, caption) in enumerate(views):
            im = capture(cap)
            hh = s * im.height / im.width
            ax = ax_mm(x0 + c * (s + gap), y, s, hh)
            ax.imshow(np.asarray(im), interpolation="none")      # vector backends embed the capture as is
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(True); sp.set_color(LGREY); sp.set_linewidth(0.5)
            T(x0 + c * (s + gap) + s / 2, y + hh + 2.0, caption, 5.5, GREY, ha="center")
        y += hh + 4.6
    return y


# --------------------------------------------------------------------------- b

def panel_timing(x0, y0, letter="a"):
    head(x0, y0, letter, "Interactive, locally and over HTTP")
    fe, sv = load(FIG6 / "frontend.json"), load(FIG6 / "server.json")
    ho = load(RESULTS / "hosted.json")
    cat = {d["id"]: d["category"] for d in load(RESULTS / "atlas.json")["datasets"]}
    by = {}
    for ds in sv["datasets"]:
        by.setdefault(ds["kind"], []).append(ds)

    def inter(kind, name):
        out = []
        for ds in by.get(kind, []):
            r = {r["interaction"]: r for r in ds["interactions"]}.get(name)
            if r:
                out.append(r["median_ms"] / 1e3)
        return out

    imaging = lambda name: inter("imaging", name) + inter("imaging-replicated", name)   # noqa: E731
    first = [np.median([l["first_frame_ms"] for l in d["loads"]]) / 1e3 for d in fe["datasets"]]
    local = [  # label, values (real first, then tiled copies), colour
        ("focus a cell", imaging("focus a cell"), INDIGO),
        ("pan, 2 IF tracks", imaging("pan locus (2 tracks, 20 Mb)"), INDIGO),
        ("pseudo-bulk, 500 cells", inter("hic", "pseudo-bulk, 500 cells, chr1"), TEAL),
        ("lasso → all views", imaging("lasso → re-scope views"), INDIGO),
        ("distance map, all cells", imaging("distance map, chr15, all cells"), INDIGO),
        ("open page → first frame", first, INDIGO),
    ]
    passes = {d["id"]: d["passes"][-1] for d in ho["datasets"] if d["passes"] and "error" not in d["passes"][-1]}
    hosted = [("open a store", "open"), ("first view", "first_view"), ("merged map, chr2", "bulk_map"),
              ("one cell's map, chr2", "one_cell"), ("all cells, pseudo-bulk", "all_cells")]
    ax = ax_mm(x0 + 30.0, y0 + 7.0, 46.0, 36.0)
    n_l, n_h = len(local), len(hosted)
    ys_l = np.arange(n_l)[::-1] + n_h + 1.4
    ys_h = np.arange(n_h)[::-1].astype(float)
    for y, (lab, vals, col) in zip(ys_l, local):
        vals = np.array(vals)
        ax.scatter(vals[:1], [y], s=9, color=col, lw=0, zorder=4)
        if len(vals) > 1:                                          # tiled copies of the same cells, larger
            ax.scatter(vals[1:], [y] * (len(vals) - 1), s=9, facecolor="white", edgecolor=col, lw=0.6, zorder=4)
            ax.plot([vals.min(), vals.max()], [y, y], color=col, lw=0.5, zorder=3)
    for y, (lab, key) in zip(ys_h, hosted):
        vals = [(p[key], KIND[cat[i]]) for i, p in passes.items() if p.get(key) is not None]
        jit = np.random.default_rng(1).uniform(-0.22, 0.22, len(vals))
        ax.scatter([v for v, _ in vals], y + jit, s=5, c=[c for _, c in vals], lw=0, alpha=0.9, zorder=4)
        med = np.median([v for v, _ in vals])
        ax.plot([med, med], [y - 0.38, y + 0.38], color=INK, lw=0.9, zorder=5)
    ax.set_yticks(list(ys_l) + list(ys_h), [r[0] for r in local] + [r[0] for r in hosted])
    ax.tick_params(axis="y", length=0)
    ax.set_xscale("log")
    ax.set_xlim(3e-3, 8)
    ax.set_ylim(-0.7, n_l + n_h + 1.3)
    ax.set_xticks([0.01, 0.1, 1], ["10 ms", "100 ms", "1 s"])
    ax.axvline(0.1, color=LGREY, lw=0.5, ls=(0, (3, 1.5)), zorder=1)
    ax.axvline(1.0, color=LGREY, lw=0.5, ls=(0, (3, 1.5)), zorder=1)
    ax.axhline(n_h + 0.25, color=LGREY, lw=0.5)
    ax.set_xlabel("response time (median of repeats)")
    # the two groups of rows: headings in the column of the row labels, above their rows
    tr = mpl.transforms.blended_transform_factory(ax.transAxes, ax.transData)
    for yy, s in [(ys_l[0] + 0.95, "Local, one laptop"), (ys_h[0] + 0.95, "Hosted, over HTTP")]:
        ax.text(-0.03, yy, s, transform=tr, fontsize=6, color=INK, fontweight="bold", va="center", ha="right")
    hs = [Line2D([], [], ls="", marker="o", ms=2.6, color=INDIGO),
          Line2D([], [], ls="", marker="o", ms=2.6, mfc="white", mec=INDIGO, mew=0.6)]
    ax.legend(hs, ["real cells", "tiled copies, to 11.5 M spots"], loc="upper right", bbox_to_anchor=(1.0, 0.93),
              frameon=False, fontsize=5.0, handletextpad=0.2, borderaxespad=0, labelspacing=0.2)
    return ax


# --------------------------------------------------------------------------- c, d

STUDY_LABEL = {"guo2026": "Hi-C-RNA", "chen2026": "spatial Hi-C", "wang2026": "ATAC-Hi-C",
               "stevens2017": "Stevens 2017", "liu2023": "HiRES", "zhou2024": "GAGE-seq",
               "chang2024": "Droplet / Paired Hi-C", "wu2025": "dscHi-C", "wei2026": "scHiCAR",
               "gao2025": "Uni-C", "takei2025": "DNA seqFISH+"}
#: label offsets (points) for the study labels, at the study's first plotted dataset
LABEL_AT = {"guo2026": [("hicrna_human_melanoma", (2, 5), "left")],
            "chen2026": [("chen2026_cortex_adult_10um", (4, -2), "left")],
            "wang2026": [("atac_hic_mouse_brain", (-4, -2), "right")],
            "stevens2017": [("stevens2017_mesc", (3, -5), "left")],
            "liu2023": [("hires_brain", (-4, 3), "right")],
            "zhou2024": [("gageseq_human_bm_cd34", (0, 6), "center")],
            "chang2024": [("droplet_hic_mouse_cortex", (-4, 2), "right"), ("paired_hic_mouse_cortex", (4, 0), "left")],
            "wu2025": [("dschic_aging_cortex", (4, -5), "left")],
            "wei2026": [("schicar_mouse_cortex", (3, -5), "left")],
            "gao2025": [("unic_human", (4, 1), "left")],
            "takei2025": [("takei2025_cerebellum", (-4, 0), "right")]}
#: datasets labelled by their own name (two of one study)
DATASET_LABEL = {"droplet_hic_mouse_cortex": "Droplet Hi-C", "paired_hic_mouse_cortex": "Paired Hi-C"}


def panel_atlas(x0, y0, letter="b"):
    at = load(RESULTS / "atlas.json")["datasets"]
    head(x0, y0, letter, "An open atlas")
    ax = ax_mm(x0 + 9.0, y0 + 7.0, 42.0, 31.5)
    shown, skipped = {}, []
    for d in at:
        cpu = (d.get("contacts_per_unit") or {}).get("median")
        y = cpu if cpu else (d.get("loci_per_cell") or {}).get("median")
        if not y:
            skipped.append(d["id"])
            continue
        shown[d["id"]] = (d["n_cells"], y)
        human = d["organism"].startswith("Homo")
        ax.scatter(d["n_cells"], y, s=5 * np.sqrt(d["size_mb"] / 100), marker="^" if human else "o",
                   color=KIND[d["category"]], alpha=0.85, lw=0.4, edgecolor="white", zorder=4)
    for st, places in LABEL_AT.items():
        for ds, off, ha in places:
            if ds in shown:
                ax.annotate(DATASET_LABEL.get(ds, STUDY_LABEL[st]), shown[ds], xytext=off, textcoords="offset points",
                            fontsize=5.0, color=INK, ha=ha, va="center")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(4, 2e5); ax.set_ylim(800, 3e7)
    ax.set_xlabel("cells or spots")
    ax.set_ylabel("contacts per cell or spot\n(traced loci for tracing)", linespacing=1.0)
    # columns: the two Hi-C kinds | tracing | the two species (open symbols together); area (store size): legend text
    hs = [Line2D([], [], ls="", marker="o", ms=3, color=KIND[k]) for k in KIND] + [Line2D([], [], ls="")]
    hs += [Line2D([], [], ls="", marker="^", ms=3, mfc="white", mec=INK, mew=0.5),
           Line2D([], [], ls="", marker="o", ms=3, mfc="white", mec=INK, mew=0.5)]
    ax.legend(hs, [KIND_LABEL[k] for k in KIND] + ["", "human", "mouse"], loc="upper center",
              bbox_to_anchor=(0.42, -0.27), frameon=False, fontsize=5.0, handletextpad=0.1, borderaxespad=0,
              labelspacing=0.2, columnspacing=0.8, ncol=3)
    return ax


PUB = [("hires_brain", "HiRES brain"), ("dschic_aging_cortex", "dscHi-C"), ("paired_hic_mouse_cortex", "Paired Hi-C"),
       ("hires_embryo", "HiRES embryo"), ("droplet_hic_mouse_cortex", "Droplet Hi-C")]


def panel_published(x0, y0, letter="c"):
    at = {d["id"]: d for d in load(RESULTS / "atlas.json")["datasets"]}
    head(x0, y0, letter, "Contacts per cell")
    ax = ax_mm(x0 + 19.5, y0 + 10.0, 17.0, 26.0)
    rows = [(lab, at[i]["published_counts"]) for i, lab in PUB if at.get(i, {}).get("published_counts")]
    y = np.arange(len(rows))[::-1]
    near = [100 * p["within_0.1pct"] for _, p in rows]
    exact = [100 * p["exact"] for _, p in rows]
    ax.barh(y, near, height=0.62, color="#cfe7e3", edgecolor=TEAL, lw=0.4, label="within 0.1%")
    ax.barh(y, exact, height=0.62, color=TEAL, lw=0, label="identical")
    ax.set_yticks(y, [f"{lab}\n{p['n']:,} cells" for lab, p in rows], linespacing=1.0)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 100); ax.set_xticks([0, 50, 100])
    ax.set_xlabel("cells (%)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.3, -0.33), frameon=False, fontsize=5.0, ncol=2,
              handlelength=1.0, handletextpad=0.3, columnspacing=0.6, borderaxespad=0)
    return ax


def make_figures() -> list:
    new_canvas(53.0)                                   # Supplementary Fig. 5
    panel_timing(0.5, 0.0)
    panel_atlas(82.0, 0.0)
    panel_published(141.0, 0.0)
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"supp_atlas.{ext}", dpi=300)
    plt.close(fig)
    return [HERE / "supp_atlas.pdf"]


if __name__ == "__main__":
    for p in make_figures():
        print(p)
