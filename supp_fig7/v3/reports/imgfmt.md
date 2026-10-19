# Version 3, family `imgfmt`: report

One Part B question (T6), four Part A questions on reading the open store (F1-F4) and four anchors (A1,
A6, H1, S3; references in `../answers_anchors.json`).  Questions in `../questions/imgfmt.yaml`, references in
`../results/imgfmt.json` (`truth/imgfmt.py` over the runs of `sherlock/imgfmt_truth.sbatch`, job 47751628, 16 CPUs;
data in `../data/imgfmt.md`).

## Summary

| ID | quantity | U-Chrom / chromdata | independent | agreement | grade | decision |
| --- | --- | --- | --- | --- | --- | --- |
| T6 | SnapFISH-IMPUTE on 4,025 held-out Sox2 loci: median / RMS error (nm) | 301.317 / 383.148 | authors' functions (sfimpute 0.0.2): 301.317 / 383.147; authors' `run_impute.py` under `mpiexec`: 302.940 / 378.852 | 2e-6; 0.5 % / 1.1 % | ± 2 % | keep, parity |
| F1 | GAGE-seq atlas store over HTTP, embedded per-cell map of mBC-0716-1_C10,E8, chr7: total, pixel (30, 31), contacts of bin 30 | 8,463 / 5 / 51 | raw Zarr reader written from `spec.md` (no chromdata): same | exact | exact | keep |
| F2 | nuclear outlines: count, area of cell 4_23, largest cell and its area (µm²) | 31 / 119.871 / 4_4, 302.845 | raw: WKB parsed by hand + shoelace: same | exact | ± 1 %, name | keep |
| F3 | provenance of `separation` (function, windows, statistic), stored d_short of 4_23, recomputation on chr7 | `uchrom.fea.separation_profile`; 0.24622; 90 pairs, 0.24069 / 189 pairs, 0.71066 | raw Parquet + numpy from the stated definition: same | exact | ± 1 %, name | keep |
| F4 | HiRES embryo atlas store over HTTP, chr19 Rg: traces, cells, median all / early neurons / blood | 14,935 / 7,469 / 10.784 / 11.191 / 9.511 | raw: one partition + two category tables: same | exact | ± 1 % | keep |

## T6: SnapFISH-IMPUTE

- **Data**: the authors' example tables (Huang et al. 2021 Sox2, allele 129: 1,416 haploids x 41 loci, 40,258
  detected); 10 % of the detected loci held out with a fixed seed.
- **Three implementations** (`truth/imgfmt_t6.py`): U-Chrom's `uchrom.im.impute.impute_snapfish`; the authors'
  functions in the order of their single-process path; the authors' MPI script `run_impute.py` as their README runs it
  (it down-scales the distance matrices before the similarity step, hence the 1 % difference).
- **The authors' environment**: `pip install sfimpute==0.0.2` resolves to pandas 3.0 / numpy 2.5, where the package
  fails (`boxcox_by1d` writes into the read-only array `to_darr()` returns: "assignment destination is read-only").
  With the releases of the package's time (2023-11: numpy 1.26.4, pandas 2.1.4, scipy 1.11.4, scikit-image 0.22.0,
  scikit-learn 1.3.2; `venv-sfimpute2023`) it runs.  An agent without U-Chrom must find this out; the prompt does not
  say it.
- **Runtimes** (16 processes): U-Chrom 162 s, the authors' functions 181 s, `run_impute.py` 131 s.
- **What the method gives here**: plain linear interpolation along the genome is better on this hold-out (median
  206.6 nm, RMS 363.1 nm) than SnapFISH-IMPUTE (301.3 / 383.1).  The question grades a faithful run of the named
  method, not the best imputation.  Linear interpolation lies outside the tolerance: 31 % (median) and 5.2 % (RMS)
  away.
- **Decision**: keep, parity, ± 2 % on both numbers.  The tolerance covers the authors' MPI path.

## F1-F4: the store read without the library

Each value is computed twice, by chromdata / U-Chrom in the agents' environment and by a reader written from
`packages/chromdata/spec.md` with zarr / pyarrow / fsspec / numpy only, in the agents' plain environment
(`venv-plain`, the condition without U-Chrom) — `truth/imgfmt_fmt.py`.  Every graded value is identical.

Times (seconds) and peak memory (MB), chromdata vs raw:

| ID | chromdata | raw |
| --- | --- | --- |
| F1 | 1.5 s, 110 MB | 1.8 s, 107 MB |
| F2 | 0.8 s, 123 MB | 0.1 s, 98 MB |
| F3 | 0.9 s, 183 MB | 0.7 s, 108 MB |
| F4 | 6.7 s, 1,792 MB | 4.1 s, 910 MB |

- **F1**: the embedded map's cell list is not in the order of the cells table (the cell is row 432 of `cells`, index
  323 of the embedded cells), so a reader that assumes the same order reads another cell.
- **F2**: the outlines are the x-y footprint of the authors' segmented nucleus voxels.  Their areas follow the
  authors' `nucleus_area_um2` (Pearson 0.99998, median ratio 1.017).
- **F3**: the prompt states the definition, so the recomputation does not need U-Chrom.  The provenance (function,
  windows, statistic) is in the store, and the build record that also names it was moved out of the agents' folder.
- **F4**: the store is several GB, and a whole-store read does not fit a session's memory.  Both readers read only the
  chr19 partition: 0.9-1.8 GB peak memory.

## Answers stated in the documentation

A scan of the snapshot's docs and skills for every graded reference value (all number formats) found: A1's counts
(383,674 spots, 1,500 traces) and H1's 399 cells in `reference/data_sources.md`, and F4's 7,469 cells in the same page
and a tutorial.  A1 and H1 are anchors and stay as they are (Part A is reported apart from the method contrasts); F4 no longer grades its
cell count.  In other families: BH2's recall / precision and the Mustache numbers BH1's count follows from (a
validation table of `bulk/structures.md`) and R4's purity (the release notes) — removed from the snapshot by
`build_snapshot.py` (`REDACTIONS.json`).

## Dropped: T7 (homolog assignment of raw seqFISH+ detections)

Jie (the authors' aligner, @358f365) and U-Chrom's assignment were compared on Takei 2021 detections with the
authors' homolog labels, chr3, chr7 and chr11 (201 cells; median ARI per cell and chromosome):

| method | median ARI | spots assigned |
| --- | --- | --- |
| plain k-means, 2 clusters per cell and chromosome | 1.0 | 100 % |
| DBSCAN, eps 1 µm | 1.0 | 100 % |
| U-Chrom, the tutorial's settings | 0.725 | 84 % |
| U-Chrom, defaults | 0.508 | 59 % |
| jie, defaults / the authors' notebook settings | 0.0 | 0.5 % / 0.07 % |

Plain clustering separates the two homologs of these 25-kb loci, so the question would not need the method.  The
authors' jie assigns almost no spots on these data.  T7 was dropped.  The exploration runs stay in
`imgfmt_work/explore/` on Sherlock.  U-Chrom's aligner also scores below plain clustering here; that is noted for the
library, not for the benchmark.
