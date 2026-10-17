# U-Chrom reproducibility

The benchmarks, results and figure scripts behind the U-Chrom paper: everything that turns the U-Chrom code
([u-science-lab/U-Chrom](https://github.com/u-science-lab/U-Chrom)) and public data into the paper's numbers,
figures and tables.  One folder per figure or table; the benchmarks that make their numbers, by topic.

```
figures/                    one folder per figure: its script(s), inputs, the rendered figure and a README
  fig1/  fig2/
  extended_data_fig1/  extended_data_fig2/  extended_data_fig3/
  supp_fig1/ ... supp_fig7/
tables/
  supp_table2/              the atlas (generated)
  supp_table3/              the toolkit against reference implementations (sources of each row)
benchmarks/                 the measurements and analyses, each writing results/*.json with provenance
  store/                    ChromData store: sizes, access, streaming, capability audit, linked maps
  browser/                  the web browser's response times
  agents/                   the agent benchmark (version 1) and the handoff session
  atlas/                    the atlas datasets' checks, the hosted browser over HTTP
  cross_assay/              the mouse cortex across single-cell Hi-C, spatial Hi-C and imaging
  spatial_variable/         spatially variable 3D genome features (melanoma, mouse brain)
  toolkit_parity/           the toolkit against the authors' code: hic_structures, sc_structures,
                            arcfish_parity, schicluster_parity, nucdyn (the reconstruction engine)
  _paper.py, _record.py     where results go; results with provenance
shared/                     drawing code and illustrations several figures use
paper_paths.py              where things are: this repository, U-Chrom, the data, the paper checkout
sync_to_paper.py            copies the rendered figures and tables into the paper checkout
env/                        the pinned U-Chrom version
```

## What makes what

| Item | Benchmark → results | Figure / table |
|---|---|---|
| Fig. 1 | (schematic; atlas counts from `benchmarks/atlas`) | `figures/fig1` |
| Fig. 2 | `benchmarks/cross_assay` | `figures/fig2` |
| Extended Data Fig. 1 | `benchmarks/store` | `figures/extended_data_fig1` |
| Extended Data Figs 2, 3 | `benchmarks/spatial_variable` | `figures/extended_data_fig2`, `figures/extended_data_fig3` |
| Supplementary Figs 1, 2 | `benchmarks/store` | `figures/supp_fig1`, `figures/supp_fig2` |
| Supplementary Figs 3-5 | (browser captures of the public atlas) | `figures/supp_fig3` ... `supp_fig5` |
| Supplementary Fig. 6 | `benchmarks/browser`, `benchmarks/atlas` | `figures/supp_fig6` |
| Supplementary Fig. 7 | `benchmarks/agents` | `figures/supp_fig7` |
| Supplementary Table 2 | `benchmarks/atlas` | `tables/supp_table2` |
| Supplementary Table 3 | `benchmarks/toolkit_parity`, `benchmarks/spatial_variable` | `tables/supp_table3` |
| Supplementary Notes | the same results (`benchmarks/store` for the format and streaming numbers) | |

## Running

1. Install U-Chrom: `pip install -r env/uchrom.txt` (the pinned commit), or work against a checkout next to this
   repository (`../U-Chrom`, or `$UCHROM_REPO`; its packages come first when it is there).
2. Data: `uchrom.datasets` fetches and builds every dataset into `$UCHROM_DATA` (default `~/.cache/uchrom`);
   U-Chrom's `apps/atlas/recipes/README.md` documents each one (source, size, licence).
3. A benchmark writes `benchmarks/<topic>/results/*.json` with the revisions of U-Chrom and of this repository,
   the host, the package versions, the data and the parameters.  Heavy runs have Sherlock job scripts
   (`benchmarks/*/sherlock/`: `REPO` is a copy of U-Chrom, `REPRO` a copy of this repository); large outputs are
   archived on Oak (`$OAK/uchrom/results/<name>/`).
4. A figure script reads only those results and writes its figure next to itself, e.g.
   `python figures/supp_fig1/make_supp_access.py`; `make tables` regenerates the generated tables.
5. `make sync` (`python sync_to_paper.py`) copies the figures and tables into the paper checkout (`../u-chrom-paper`,
   or `$UCHROM_PAPER`); build the manuscript there.

Tests of the benchmark helpers: run each folder on its own (`pytest benchmarks/tests`, `pytest benchmarks/store/tests`,
`pytest benchmarks/spatial_variable/tests`; the folders share the package name `tests`).  The committed figures were
rendered on macOS with matplotlib's default backend; `MPLBACKEND=Agg` changes the PNGs' pixel size by one or two
pixels, and scatter panels that subsample points differ in which points are drawn.
