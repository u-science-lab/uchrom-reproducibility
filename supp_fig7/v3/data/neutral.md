# Version 3: the stores of T1, P1 and P2 in community formats

T1, P1 and P2 read data that the shared data root holds only as ChromData stores (`takei`, `melanoma`).  An agent without U-Chrom
then has to read the store's Zarr + Parquet layout by hand, an obstacle the questions do not test.  Version 3 gives
the same data in community formats as well, in two folders under `$SCRATCH/uchrom-nd/agentbench/data_v3/`.  The
prompts name both formats, and every condition gets both: the store and the export.  `neutral_build.py`
(this folder) writes the exports with chromdata's own exporters, in the job `../sherlock/neutral_build.sbatch`
(Sherlock jobs 47856015 and, for the Takei folder with the exact-parse check, 47858981; 2026-10-14;
`venv-uchrom-v3`, chromdata of b6d8873).  It checks every export against its store value by value, and
recomputes the question's statistic from the files.  The references do not change.

| data key | folder: files | questions | content |
| --- | --- | --- | --- |
| `takei_fofct` | `takei2025_fov0_fofct/`: `takei2025_fov0_fofct_core.csv`, `takei2025_fov0_fofct_cells.csv` | T1 (with `takei`), P2 (with `takei`, `takei_rna`) | 4DN FOF-CT core table + cell table of the `takei` store |
| `melanoma_r1` | `hicrna_melanoma_r1/`: `human_melanoma_r1.spots_1Mb.scool`, `human_melanoma_r1_spots.csv` | P1 (with `melanoma`) | per-spot 1-Mb maps of section `human_melanoma_r1` of the `melanoma` store + a spot table |

## `takei2025_fov0_fofct/`

- `takei2025_fov0_fofct_core.csv`: `ChromData.to_fofct(spot_tracks=False)` of the shared store `takei2025_fov0`
  (Takei et al. 2025, cerebellum, field of view 0).  FOF-CT header lines (`##`/`#`), then
  `##columns=(Spot_ID,Trace_ID,X,Y,Z,Chrom,Chrom_Start,Chrom_End,Cell_ID, ...)` and one row per spot; X/Y/Z in µm.
  The store's spot tracks are not in the table; no question uses them.
- `takei2025_fov0_fofct_cells.csv`: the FOF-CT cell table of the same call, with every column of the store's `cells`
  (`Cell_ID` first, `cell_type` among them).
- Checks (`_build_info/neutral_build_takei.json`): the two tables read with pandas alone have one row per spot and
  per cell.  X, Y, Z, Trace_ID, Cell_ID, Chrom, Chrom_Start, Chrom_End and cell_type equal the store's values exactly.
  The coordinates are float64 in the store.  The table holds their shortest round-trip decimals, and pandas reads
  them back exactly with `float_precision="round_trip"`.  Its default fast parser can be one unit in the last place
  off (at most 5.7e-14 µm here).  P2's two correlations recomputed from the files (pandas / numpy) equal
  `results/bulk.json` to 1e-12.
- T1 from the files (`neutral_t1_check.py`, job 47857826): T1's input (chr15, 1-Mb bins) built from the core table
  with pandas' default parser, then the authors' ArcFISH (a21d889) and U-Chrom on it.  Both give the reference: 8
  domains, boundaries at 7, 11, 25, 43, 59, 66 and 72 Mb.  The binned coordinates differ from the reference's input
  by at most 5.7e-14 µm.
- Size: 39,240,100 bytes (core, md5 `a427cedf725acf1454c218b9d56d9abe`) + 47,797 bytes (cells, md5
  `fdb08d602967355853da986d4c51fe8d`); 383,674 spots, 47 cells; columns
  `Spot_ID,Trace_ID,X,Y,Z,Chrom,Chrom_Start,Chrom_End,Cell_ID,name` (core), `Cell_ID`, `leiden`, `cell_type`, `fov`,
  `Cent_ROI_x/y/z`, `nuc_volume_um3`, `doublet`, `batch` and 48 `if.*` immunofluorescence columns (cells).

## `hicrna_melanoma_r1/`

- `human_melanoma_r1.spots_1Mb.scool`: `chromdata.embedded.export_scool(store, "per_spot.human_melanoma_r1",
  cells=<the section's spots>, trans=True)` of the atlas store `hicrna_human_melanoma` (Guo et al. 2026, GEO
  GSE311199).  It holds one cooler per spot, 1-Mb bins, hg38, with the contacts within and between chromosomes.
  P1's reference makes the same export.
- `human_melanoma_r1_spots.csv`: `cell_id` (the cooler's name under `/cells/`), `array_row`, `array_col` of the
  section's spots (the store's `cells`).
- Checks (`_build_info/neutral_build_melanoma.json`): one cooler per spot of the section, 1-Mb bins.  P1's three
  per-spot tables read from the scool with cooler alone equal those read from the store's embedded map (chromdata's
  reader, `truth/bulk.py::_melanoma_legs_fast`): the legs per chromosome, the counts at bin offset 1 and those at
  offsets 19-53.  P1's observed R² recomputed from the two files equals `results/bulk.json` (0.343917).
- Size: 1,102,779,654 bytes (scool, md5 `65208a714b600d995bd9d390700b208f`) + 256,578 bytes (spots, md5
  `2227e756a12b3b784b5fcc615c6977e5`); 6,268 spots.  The export takes 328 s, and reading every cooler with cooler
  270 s (4 CPUs).

Registration: no new source.  Both folders are derived from shared data (`takei2025_fov0`, the atlas store
`hicrna_human_melanoma`), which are registered.

## Pilot

In the first pilot (`runs/pilot_v3c`, gpt-6-luna, 1 repeat, store only), the condition without U-Chrom failed on two
of the three questions because of the store's layout, not the method.  On T1 the agent found no coordinates in the
Zarr store ("that table has no genomic coordinates or x/y/z columns") and gave no answer.  On P2 it misread the
store and reported Spearman -0.110 and partial 0.073 (reference 0.827 / 0.687).  It solved P1.  With both formats
(`runs/pilot_v3d_neutral`, job 47859138, graded by 47859267, 2026-10-14), all six sessions (none, skills x T1, P1,
P2) scored 1.0:

- **none, T1**: read the FOF-CT core table and used ArcFISH's code from GitHub;
- **none, P2**: read the FOF-CT core table;
- **none, P1**: read the `.scool`;
- **skills, T1 and P2**: used U-Chrom on the store and compared it with the FOF-CT;
- **skills, P1**: read the `.scool` with cooler, without U-Chrom.
