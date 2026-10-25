# U-Chrom reproducibility

The code that turns the U-Chrom software ([u-science-lab/U-Chrom](https://github.com/u-science-lab/U-Chrom)) and
public data into the paper's numbers: the benchmarks and analyses, and the results they wrote.  One folder per figure
or table whose numbers come from data; a benchmark that several items use sits in the folder of the first one.  The
figures themselves are drawn from these results in the paper's repository.

```
fig2/                       the mouse cortex across single-cell Hi-C, spatial Hi-C and imaging (cross-assay analysis)
extended_data_fig1/         the ChromData store: sizes, access, streaming, capability audit, linked maps
                            (also Supplementary Figs 1, 2)
extended_data_fig2/         spatially variable 3D genome features, melanoma and mouse brain (also Extended Data Fig. 3)
supp_fig6/                  browser/: the web browser's response times; atlas/: the atlas datasets' checks and the
                            hosted browser over HTTP (also Supplementary Table 2, Fig. 1)
supp_fig7/                  the agent benchmark and the handoff session
supp_table3/                the toolkit against the authors' code: hic_structures, sc_structures, arcfish_parity,
                            schicluster_parity, nucdyn (the reconstruction engine)
shared/                     paper_paths.py (where things are), _record.py (results with provenance), tests
env/                        the pinned U-Chrom version
```

Each folder has a README (what each script measures, what else uses its results) and a `results/` folder (in
`supp_fig6/` and `supp_table3/`, one per benchmark) with a `MANIFEST.md`.  Fig. 1 and Supplementary Figs 3-5 (a
schematic and browser views of the public atlas) need no code here.

## Running

1. Install U-Chrom: `pip install -r env/uchrom.txt` (the pinned commit), or work against a checkout next to this
   repository (`../U-Chrom`, or `$UCHROM_REPO`; its packages come first when it is there).
2. Data: `uchrom.datasets` fetches and builds every dataset into `$UCHROM_DATA` (default `~/.cache/uchrom`);
   U-Chrom's `apps/atlas/recipes/README.md` documents each one (source, size, licence).
3. Run a benchmark from the repository root, e.g. `python extended_data_fig1/c_access.py`.  It writes
   `results/*.json` with the revisions of U-Chrom and of this repository, the host, the package versions, the data
   and the parameters (`shared/_record.py`; `python shared/_record.py OLD NEW --rtol x` compares two runs).  Heavy runs
   have cluster job scripts (`*/slurm/`: `REPO` is a copy of U-Chrom, `REPRO` a copy of this repository); large
   outputs are archived in long-term storage (`$UCHROM_ARCHIVE/results/<name>/`).

Tests of the helpers: run each folder on its own (`pytest shared/tests`, `pytest extended_data_fig1/tests`,
`pytest extended_data_fig2/tests`; the folders share the package name `tests`).
