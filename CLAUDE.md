# uchrom-reproducibility -- notes for agents

* This repository holds only the code that makes the paper's numbers from data, and its results: one folder per
  figure or table whose numbers come from data (`fig2/`, `extended_data_fig1/`, ...; `shared/paper_paths.py` maps
  each benchmark's topic to its folder), each writing `results/*.json`.  No rendered figures, no figure or table
  layout code, no illustrations: those live in the paper repository (`../u-chrom-paper`), whose figure scripts read
  the results here (`$UCHROM_REPRO`).  An item without data code has no folder.  Material outside the paper is kept
  off the repositories.  U-Chrom keeps the library code.
* Commits and pull requests here carry no Claude attribution (no Co-Authored-By trailer, no generated-with footer).
* Numbers come only from the results files.  When a number changes, re-run the benchmark, then redraw the figure and
  update the text in the paper repository.
* A change the library needs is made in U-Chrom first (its PR); this repository then runs against it
  (`$UCHROM_REPO` checkout while developing, `env/uchrom.txt` pin for a submitted version).
* U-Chrom's validation discipline applies: real data, the original papers' datasets and metrics, disappointing
  numbers reported as they are; every dataset registered in U-Chrom (`uchrom.datasets`, recipes README).
* Results record their provenance (revisions of U-Chrom and of this repository, host, versions, data, parameters);
  committed results files are never edited by hand.  Heavy runs go to the cluster (`*/slurm/`, `REPO` = U-Chrom copy,
  `REPRO` = copy of this repository); the laptop only for the timings the paper reports from it, when the user asks.
  Large outputs are archived in long-term storage, never only on `$SCRATCH`.
* Talk to the user in Chinese; write code, commits and documents in English.
