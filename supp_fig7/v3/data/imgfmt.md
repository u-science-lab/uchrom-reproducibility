# Family `imgfmt`: data of the version-3 questions T6, T7, F1-F4 (+ anchors A1, A6, H1, S3)

Built by `v3/data/imgfmt_build.py` inside `v3/sherlock/imgfmt_setup.sbatch` (2026-10-13), from registered sources only,
into `$SCRATCH/uchrom-nd/agentbench/data_v3/<folder>/`.  F1 and F4 read public atlas stores over HTTP (no folder).
The anchors use the shared data keys (`takei`, `schicar`, `atlas_hires`, `melanoma`).

## Data keys (`run_codex.py` `DATA` / `URLS`)

| Key | Folder (bound at `./data/<folder>`) / URL | File the prompt names | Questions |
|---|---|---|---|
| `sox2_holdout` | `imgfmt_sox2/` | `""` (a folder: `sox2_coor.tsv`, `sox2_ann.tsv`, `holdout.tsv`) | T6 |
| `seqfish_det` | `imgfmt_seqfish/` | `""` (a folder: `detections.tsv`, `labels.tsv`) | T7 (dropped; folder kept) |
| `takei_fov4` | `imgfmt_takei2021_fov4/` | `takei2021_fov4.chromdata.zarr` | F2, F3 |
| `atlas_gageseq` | URL `https://uchrom-atlas-r2.u-science.org/gageseq_mouse_cortex.chromdata.zarr` | (served over HTTP) | F1 |
| `atlas_hires_embryo` | URL `https://uchrom-atlas-r2.u-science.org/hires_embryo.chromdata.zarr` | (served over HTTP) | F4 |

## The folders as built

The build records (`build_info.json`) were moved out of the folders the agents see, to
`data_v3/_build_info/<folder>.json` (the F2 / F3 one names the stored results' function and parameters, which F3 asks
the agent to read from the store).

- `imgfmt_sox2/` (3.4 MB; `imgfmt_build.py sox2` in job 47737576 of `imgfmt_setup.sbatch`, which built both T6 / T7
  folders and was then cancelled in the IF download, moved to `imgfmt_outlines.sbatch`): source `huang2021_sox2` (SnapFISH-IMPUTE @7f7f1a7,
  md5-checked), the 129 allele: 1,416 haploids x 41 loci = 58,056 rows, 40,258 detected; `holdout.tsv`: 4,025 of
  them (10 %, `numpy.random.default_rng(0).choice` over the detected rows in file order); every haploid keeps at least
  12 detected loci; locus 26 is never detected.  `sox2_coor.tsv` keeps the true coordinates of the held-out loci (the
  prompt tells the agent to hide them).
- `imgfmt_seqfish/` (24 MB): 316,995 raw detections of 201 cells, 20 chromosomes x 60 loci, with the authors' homolog
  labels (T7, dropped).
- `imgfmt_takei2021_fov4/takei2021_fov4.chromdata.zarr` (1.7 MB; job 47738721, `imgfmt_outlines.sbatch`): sources
  `takei` (4DNFIHF3JCBY), `takei_tables` (4DNFIFINA2U9) and the Zenodo 3735329 IF zip of field of view 4 (2.47 GB,
  md5-checked, read in the job): 31 cells, 1,272 traces, 56,374 spots; `cell_shapes['nucleus']` of all 31 cells (the
  x-y footprint of the authors' segmented nucleus voxels, 103-nm pixels; area vs the authors' `nucleus_area_um2`:
  Pearson 0.99998, median ratio 1.017); results `separation` and `separation_mean` of
  `uchrom.fea.separation_profile` with their provenance.

