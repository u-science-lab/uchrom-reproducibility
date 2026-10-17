# Environment

The benchmarks and figure scripts run with the U-Chrom packages (`chromdata`, `uchrom-browser`, `u-chrom`;
the native `uchrom-recon` / `uchrom-maps` where a benchmark says so).

* **While developing**: a U-Chrom checkout next to this repository (`../U-Chrom`, or `$UCHROM_REPO`) and its
  `.venv` (`uv pip install -e "packages/uchrom[web]"` there).  The benchmarks put that checkout's packages ahead of
  installed copies when it is there (`paper_paths.UCHROM_REPO`).
* **For a submitted version**: `pip install -r env/uchrom.txt`, which pins the U-Chrom commit the numbers were
  made with.  Move the pin forward when the results are re-made with a newer version, and record it in the
  manuscript's Code availability section.  Every results file records the revision it was made with
  (`machine.versions.uchrom_git`, and `repro_git` for this repository).

`uchrom.txt` holds U-Chrom `main` at the split of this repository from the code and paper repositories
(2026-10-17); the earlier pin (aac3961, 2026-10-05) named the repository root, which holds no package.
