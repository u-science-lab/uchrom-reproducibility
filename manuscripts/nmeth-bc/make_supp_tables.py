"""Brief Communication, Supplementary Table (the atlas catalogue) -> ``supp/table_atlas.tex`` (rows of a longtable).

Reads ``figures/results/atlas.json`` (U-Chrom ``benchmarks/atlas/validate.py``: counts, sizes, the checks against the
published numbers) and the study records of the U-Chrom atlas (``apps/atlas/datasets.json`` in the U-Chrom checkout,
``paper_paths.UCHROM_REPO``: citation, journal, kind of data, source); each study is cited from ../shared/refs.bib.
The spatial checks come from ``atlas.json``'s ``spatial_checks`` (the results of U-Chrom
``benchmarks/spatial_hicrna_validation.py`` and ``benchmarks/spatial_hic_chen2026_reprocessed.py``); an embedding's
agreement is the best of the store's embeddings of that modality (PCA, UMAP, t-SNE or LSI).

    UCHROM_REPO=/path/to/U-Chrom python manuscripts/nmeth-bc/make_supp_tables.py
"""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))  # this repository's root (paper_paths)
from paper_paths import UCHROM_REPO  # noqa: E402

KIND = {"Spatial": "spatial Hi-C", "Single-cell": "single-cell Hi-C", "Imaging": "chromatin tracing"}
# The atlas study ids that differ from their keys in ../shared/refs.bib
BIBKEY = {"liu2023": "hires", "zhou2024": "gageseq", "chang2024": "droplethic", "wu2025": "dschic",
          "wei2026": "schicar", "gao2025": "unic"}
# Embeddings of the stores that are the authors' own but carry no published_ / paper_ prefix in cellm
AUTHORS_EMBEDDINGS = {"takei2025_cerebellum": {"umap"}}  # the published RNA UMAP (apps/atlas/recipes/README.md)


def tex(s: str) -> str:
    s = re.sub(r"(?<![\w-])([23])-D(?![\w-])", r"\1D", str(s))     # house style: 3D, 2D
    s = re.sub(r"\b(E\d+(?:\.\d+)?)-(E\d)", r"\1–\2", s)                # stage ranges: en dash (E7.0–E9.5)
    return (s.replace("\\", "\\textbackslash{}").replace("&", "\\&").replace("%", "\\%").replace("_", "\\_")
            .replace("#", "\\#").replace("–", "--"))


def rng(xs, fmt="{:.2f}") -> str:
    """A value, or a range when the values differ."""
    lo, hi = fmt.format(min(xs)), fmt.format(max(xs))
    return lo if lo == hi else f"{lo}--{hi}"


#: atlas.json spatial_checks: the section keys of the Guo et al. validation, per store
GUO_SECTIONS = {"hicrna_mouse_brain": ["mouse_brain_r1", "mouse_brain_r2"],
                "hicrna_mouse_embryo": ["mouse_embryo_E11.5", "mouse_embryo_E13.5_r1", "mouse_embryo_E13.5_r2"],
                "hicrna_human_melanoma": ["human_melanoma_r1", "human_melanoma_r2"]}
#: stores whose checks are only described in the Supplementary Note (from the recipes' README)
NOTE_ONLY = {"atac_hic_mouse_brain", "stevens2017_mesc"}


def spatial(d: dict, sc: dict) -> list[str]:
    """The checks of a spatial store that need no published number: RNA against contacts on the same spots and
    replicate sections (Guo et al.); per-spot totals on the grid and the bulk compartments against the authors' A/B
    (Chen et al.)."""
    out = []
    g = sc.get("spatial_hicrna_validation.json", {})
    if d["id"] in GUO_SECTIONS:
        same = g.get("same_spots", {})
        rho = [same[k]["spearman_rna_umi_vs_hic_contacts"] for k in GUO_SECTIONS[d["id"]] if k in same]
        if rho:
            out.append(f"RNA UMIs vs contacts per spot, Spearman {rng(rho)}")
        rep = [v["median_abs_r"] for k, v in g.get("replicates", {}).items()
               if k.split(" vs ")[0] in GUO_SECTIONS[d["id"]]]
        if rep:
            out.append(f"replicate sections' compartments, median $|r|$ {rng(rep, '{:.3f}')}")
    c = sc.get("spatial_hic_chen2026_reprocessed.json", {})
    secs = [v for k, v in c.items() if d["id"] == "chen2026_" + k.split("/")[0]]
    if secs:
        shuffled = math.ceil(100 * max(abs(v["moran_shuffled"]) for v in secs)) / 100
        out.append(f"per-spot contacts on the grid, Moran's $I$ {rng([v['moran_log_total_contacts'] for v in secs])} "
                   f"(shuffled spots, $|I| \\le {shuffled:.2f}$); bulk compartments vs the authors' A/B, median $|r|$ "
                   f"{rng([v['bulk_E1_vs_ssAB_median_abs_r'] for v in secs])}")
    return out


def check(d: dict, sc: dict | None = None) -> str:
    """The checks of one dataset, in words."""
    out = []
    pc = d.get("published_counts")
    if pc:
        ex, w = 100 * pc["exact"], 100 * pc["within_0.1pct"]
        out.append(f"contacts per cell vs published: {ex:.0f}\\,\\% identical" + ("" if ex >= 99.95 else
                                                                                  f", {w:.1f}\\,\\% within 0.1\\,\\%"))
    cb = d.get("cells_vs_bulk")
    if cb:
        scc = cb["scc"]
        own = {k: v[k] for k, v in scc.items() if k in v}
        other = [x for k, v in scc.items() if k in v for kk, x in v.items() if kk != k]
        none = [x for k, v in scc.items() if k not in v for x in v.values()]       # samples without a library
        s = "merged cells vs bulk (SCC): own library " + (rng(list(own.values())) if len(own) > 1 else
                                                          f"{next(iter(own.values())):.2f} ({next(iter(own))})")
        if other:
            s += f", other libraries {rng(other)}"
        if none:
            s += f"; samples without a library of their own, {rng(none)} against it"
        out.append(s)
    ea = d.get("embedding_agreement")
    if ea:
        ag = ea["agreement"]
        parts = []
        own = AUTHORS_EMBEDDINGS.get(d["id"], set())
        for name, keys in (("RNA", [k for k in ag if k.startswith("rna_")]), ("Hi-C", [k for k in ag if k.startswith("hic_")]),
                           ("published", [k for k in ag if k.startswith(("published_", "paper_")) or k in own])):
            if keys:
                parts.append(f"{name} {max(ag[k] for k in keys):.2f}")
        if not parts:
            parts = [f"{max(ag.values()):.2f}"]
        lab = "cell types" if ea["label"] == "cell_type" else "samples"
        out.append(f"embedding vs the authors' {lab} (kNN agreement): {', '.join(parts)} (chance {ea['chance']:.2f})")
    out += spatial(d, sc or {})
    if not out and d["id"] in NOTE_ONLY:
        out.append("Supplementary Note~\\ref{snote:atlas}")
    return "; ".join(out) if out else "--"


def main() -> Path:
    res = json.loads((HERE / "figures" / "results" / "atlas.json").read_text())
    # the Hi-C embeddings recomputed with the corrected scHiCluster (2026-10-15, atlas_hic_embeddings.json; the
    # RNA and published ones are unchanged): their agreements replace those of atlas.json (2026-10-05)
    hic = HERE / "figures" / "results" / "atlas_hic_embeddings.json"
    if hic.exists():
        new = json.loads(hic.read_text())["results"]
        for d in res["datasets"]:
            if d.get("id") in new and d.get("embedding_agreement"):
                d["embedding_agreement"]["agreement"].update(new[d["id"]]["agreement_new"])
    atlas = json.loads((UCHROM_REPO / "apps" / "atlas" / "datasets.json").read_text())
    studies = {s["id"]: s for s in atlas["studies"]}
    rows = []
    for d in res["datasets"]:
        st = studies[d["study"]]
        n = f"{d['n_cells']:,} " + ("spots" if d["unit"] == "spot" else "cells")
        if d["n_spots"]:
            n += f" ({d['n_spots']:,} spots)"
        rows.append(" & ".join([
            tex(d["name"]) + f" \\newline\\path{{{d['id']}}}",
            f"{tex(st['citation'])}~\\citep{{{BIBKEY.get(st['id'], st['id'])}}}, \\emph{{{tex(st['journal'])}}} \\newline {tex(st['source'])}",
            tex(KIND.get(d["category"], d["category"])) + " \\newline " + tex(", ".join(d["modalities"])),
            f"\\emph{{{tex(d['organism'])}}} ({tex(d['assembly'])})",
            n,
            f"{d['size_mb'] / 1e3:.2f}",
            str(d["contact_maps"]),
            check(d, res.get("spatial_checks"))]) + " \\\\")
    out = HERE / "supp" / "table_atlas.tex"
    out.parent.mkdir(exist_ok=True)
    out.write_text(f"% Generated by make_supp_tables.py from figures/results/atlas.json ({res.get('generated', '')}) and "
                   "U-Chrom apps/atlas/datasets.json; do not edit.\n" + "\n".join(rows) + "\n")
    return out


if __name__ == "__main__":
    print(main())
