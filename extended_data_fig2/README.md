# Extended Data Fig. 2: spatially variable 3D genome features

Also used by Extended Data Fig. 3 (the mouse brain), Supplementary Table 3 (the spatial statistics row) and
Supplementary Note 4.

| Script | What it does | Writes |
|---|---|---|
| `spatial_variable.py` | per-spot 3D genome features of the melanoma and mouse-brain spatial Hi-C, their spatial statistics, the nulls (`patch_null`, `cnv_null`) and the local structure scores | `results/*_<tissue>.json`, `spatial_variable_<tissue>*.parquet` / `.npz` |
| `bulk_calls.py` | boundaries, loops and compartments of each tissue's pooled map: the structures the local scores are measured on | the run's output folder |
| `validate_stats.py` | SPARK-X and Moran's I / Geary's C against `SPARK::sparkx` and `esda` on real spot RNA counts | `results/validate_stats.json` |
| `check_patch_composition.py` | whether the melanoma patch's higher short / long ratio is a matter of cell types | `results/spatial_variable_melanoma_patch_composition.json` |

Heavy runs: `slurm/` (`REPO` = a U-Chrom copy, `REPRO` = a copy of this repository).  `results/MANIFEST.md` lists
the committed results.  Tests: `pytest extended_data_fig2/tests`.
