# U-Chrom reproducibility

The benchmarks, results and figure scripts behind the U-Chrom papers: everything that turns the U-Chrom code
([u-science-lab/U-Chrom](https://github.com/u-science-lab/U-Chrom)) and public data into the numbers and figures
of the manuscripts ([u-science-lab/u-chrom-paper](https://github.com/u-science-lab/u-chrom-paper), which keeps the
text and the rendered figures only).  Split out of both on 2026-10-17 (U-Chrom `main` 28409ca, u-chrom-paper
`main` ce0ddba); the history before that is in those repositories.

```
benchmarks/               the analyses and measurements that make the papers' numbers
  fig2/                   ChromData store: sizes, access, streaming, capability audit, linked maps
  fig6/                   web browser response times; agent benchmark (agent_bench/, version 1)
  cross_assay/            the mouse cortex across single-cell Hi-C, spatial Hi-C and imaging (nmeth-bc Fig. 2)
  spatial_variable/       spatially variable 3D genome features (nmeth-bc Extended Data Figs 2, 3)
  spatial_composition/    spot compartments against cell-type composition
  _paper.py, _record.py   where results go (this repository); results with provenance
manuscripts/<ms>/figures/ figure scripts, their results (results/*.json) and inputs, as laid out in the paper
manuscripts/nmeth-bc/make_supp_tables.py, manuscripts/full/make_tableS1.py   generated tables
paper_paths.py            where this repository finds U-Chrom, the data and the paper checkout
sync_to_paper.py          copies the rendered figures and generated tables into the paper checkout
env/                      the pinned U-Chrom version
```

## Running

1. Install U-Chrom: `pip install -r env/uchrom.txt` (the pinned commit), or work against a checkout next to this
   repository (`../U-Chrom`, or `$UCHROM_REPO`; its packages come first when it is there).
2. Data: `uchrom.datasets` fetches and builds every dataset into `$UCHROM_DATA` (default `~/.cache/uchrom`);
   U-Chrom's `apps/atlas/recipes/README.md` documents each one (source, size, licence).
3. A benchmark writes `manuscripts/<ms>/figures/<fig>/results/*.json` here, with the revisions of U-Chrom and of this
   repository, the host, the package versions, the data and the parameters.  Heavy runs have Sherlock job scripts
   (`benchmarks/*/sherlock/`: `REPO` is a copy of U-Chrom, `REPRO` a copy of this repository); large outputs are
   archived on Oak (`$OAK/uchrom/results/<name>/`).
4. A figure script reads only those results and writes its figure next to itself:
   `python manuscripts/nmeth-bc/figures/make_supp_access.py`.
5. `make full` (the full manuscript's figure files) and `make tables` (generated tables) where needed, then
   `make sync` (`python sync_to_paper.py`) copies the figures and tables the paper tracks into the paper checkout
   (`../u-chrom-paper`, or `$UCHROM_PAPER`); build the manuscript there.

Tests of the benchmark helpers: run each folder on its own (`pytest benchmarks/tests`, `pytest benchmarks/fig2/tests`,
`pytest benchmarks/spatial_variable/tests`; the folders share the package name `tests`).  The committed figures were
rendered on macOS with matplotlib's default backend; forcing `MPLBACKEND=Agg` changes the PNGs' pixel size by one or
two pixels, and scatter panels that subsample points differ in which points are drawn.

## What makes what (Brief Communication, `manuscripts/nmeth-bc`)

| Item | Benchmark → results | Figure script |
|---|---|---|
| Fig. 1 | schematic; browser captures | `figures/make_fig1.py` |
| Fig. 2 | `benchmarks/cross_assay` → `nmeth-bc/figures/results/` | `make_fig2.py`, `capture_fig2a.py`, `capture_fig2_structures.py` |
| Extended Data Fig. 1 | `benchmarks/fig2` → `full/figures/fig2/results/` | `make_ed_store.py` (`full/figures/fig2/make_fig2.py`) |
| Extended Data Figs 2, 3 | `benchmarks/spatial_variable` → `nmeth-bc/figures/results/` | `make_ed_melanoma.py`, `make_ed_brain.py`, `capture_ed_melanoma.py` |
| Supplementary Fig. 1 | `benchmarks/fig2` (`c_access.py`, `e_contacts.py`, `e_contacts_query.py`) | `make_supp_access.py` |
| Supplementary Fig. 2 | `benchmarks/fig2` (`f_crossmodal.py`, `e_contacts_query.py`) | `make_supp_linked.py` |
| Supplementary Figs 3-5 | whole-page browser captures | `capture_supp_browser.py` |
| Supplementary Fig. 6 | `benchmarks/fig6` (response times); the atlas (U-Chrom `apps/atlas`) | `make_supp_browser_atlas.py` |
| Supplementary Fig. 7 | `benchmarks/fig6/agent_bench` → `full/figures/fig6/results/` | `make_supp_agents.py` |
| Supplementary Table 2 | atlas checks (U-Chrom `apps/atlas/benchmarks`: `validate.py`, `hic_embeddings.py`) | `make_supp_tables.py` |

## Kept in U-Chrom

Benchmarks that validate the library itself stay with the code, run again when it changes, and keep their
results there (`benchmarks/<name>/results/`): the NucDynamics engine against the original program
(`benchmarks/fig3`, which also makes the full manuscript's Fig. 3), EMber's pre-registered test (`screcon`), parity
with the authors' code of the structure callers (`arcfish_parity`, `schicluster_parity`, `hic_structures`,
`sc_structures`), format conversions and streaming (`cd2*_*.py`), and the atlas checks (`apps/atlas/benchmarks`).
The papers cite their results from there (Supplementary Table 3, Supplementary Notes 2 and 3).  The ones that write
a figure's results (`fig3`, `cd21_streaming_validation.py`, `cd22_convert_validate.py`) write them into a checkout
of this repository.
