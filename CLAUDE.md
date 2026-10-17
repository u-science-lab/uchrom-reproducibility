# uchrom-reproducibility -- notes for agents

* This repository makes the paper's numbers, figures and tables: `benchmarks/<topic>/` -> `benchmarks/<topic>/results/
  *.json` -> `figures/<item>/` and `tables/<item>/` -> `python sync_to_paper.py` into the paper checkout
  (`../u-chrom-paper`, `$UCHROM_PAPER`).  One folder per figure or table of the current paper; nothing else lives
  here (material outside the paper is kept off the repositories).  The paper repository keeps the text and the
  rendered figures only; U-Chrom keeps the code.
* Commits and pull requests here carry no Claude attribution (no Co-Authored-By trailer, no generated-with footer).
* Numbers come only from the results files.  Never type a number into a figure script, a caption or a manuscript;
  when a number changes, re-run the benchmark, rebuild the figure, sync, and update the text.
* A change the library needs is made in U-Chrom first (its PR); this repository then runs against it
  (`$UCHROM_REPO` checkout while developing, `env/uchrom.txt` pin for a submitted version).
* U-Chrom's validation discipline applies: real data, the original papers' datasets and metrics, disappointing
  numbers reported as they are; every dataset registered in U-Chrom (`uchrom.datasets`, recipes README).
* Results record their provenance (revisions of U-Chrom and of this repository, host, versions, data, parameters).
  Heavy runs go to Sherlock (`benchmarks/*/sherlock/`, `REPO` = U-Chrom copy, `REPRO` = copy of this repository);
  the laptop only for the timings the paper reports from it, when the user asks.  Large outputs are archived on
  Oak, never only on `$SCRATCH`.
* Talk to the user in Chinese; write code, commits and documents in English.
