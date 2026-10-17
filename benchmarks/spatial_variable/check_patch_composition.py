"""Is the melanoma patch's higher short / long ratio a matter of cell types?  (Fig. 2; Online Methods, Statistics.)

The patch (top fifth of the smoothed second copy-number component, ``benchmarks/spatial_variable`` stage
``patch``) holds more spots the authors label as melanoma than the rest of the section (r1).  This compares the
destriped short / long ratio of the patch and the rest (median difference, two-sided Mann-Whitney U):

* as measured; after depth and the spot's RNA (the residuals of the regression on log10 contacts and the 30 RNA
  principal components, ``rna30.*``); after depth, RNA and copy number (``rna30_cnv.*``);
* within the spots labelled as one of the four melanoma states, and within the other spots (r1; the authors'
  cell types are published for r1 only).

Reads ``results/spatial_variable_melanoma_spots.parquet`` and writes
``results/spatial_variable_melanoma_patch_composition.json``.

    python benchmarks/spatial_variable/check_patch_composition.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import scipy
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
SPOTS = HERE / "results" / "spatial_variable_melanoma_spots.parquet"
OUT = HERE / "results" / "spatial_variable_melanoma_patch_composition.json"
RATIO = "q30.short_long_ratio"
STEPS = {"measured": f"destriped.{RATIO}", "after_depth_rna": f"rna30.{RATIO}",
         "after_depth_rna_cnv": f"rna30_cnv.{RATIO}"}


def compare(a: pd.Series, b: pd.Series) -> dict:
    a, b = a.dropna(), b.dropna()
    if len(a) == 0 or len(b) == 0:
        return {"n_patch": len(a), "n_rest": len(b)}
    return {"n_patch": len(a), "n_rest": len(b), "median_patch": round(float(a.median()), 4),
            "median_rest": round(float(b.median()), 4), "median_difference": round(float(a.median() - b.median()), 4),
            "mannwhitney_p": float(mannwhitneyu(a, b, alternative="two-sided").pvalue)}


def main() -> Path:
    spots = pd.read_parquet(SPOTS)
    out = {"generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "script": "benchmarks/spatial_variable/check_patch_composition.py",
           "input": {"file": SPOTS.name, "md5": hashlib.md5(SPOTS.read_bytes()).hexdigest()},
           "versions": {"pandas": pd.__version__, "scipy": scipy.__version__}, "sections": {}}
    for section, t in spots.groupby("section", observed=True):
        patch = t["in_patch"].astype(bool)
        labelled = t["cell_type"].notna() & (t["cell_type"].astype(str) != "nan")
        melanoma = t["cell_type"].astype(str).str.startswith("Melanoma")
        res = {"n_spots": int(len(t)), "n_patch": int(patch.sum()),
               "steps": {k: compare(t.loc[patch, c], t.loc[~patch, c]) for k, c in STEPS.items()}}
        if labelled.any():
            res["melanoma_share"] = {"patch": round(float(melanoma[patch].mean()), 3),
                                     "rest": round(float(melanoma[~patch].mean()), 3)}
            col = STEPS["measured"]
            res["within_cell_types"] = {
                "melanoma": compare(t.loc[patch & melanoma, col], t.loc[~patch & melanoma, col]),
                "other": compare(t.loc[patch & ~melanoma & labelled, col], t.loc[~patch & ~melanoma & labelled, col])}
        out["sections"][str(section)] = res
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    return OUT


if __name__ == "__main__":
    print(main())
