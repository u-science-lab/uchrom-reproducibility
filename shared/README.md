# Shared helpers

* `paper_paths.py`: where things are -- each benchmark's folder and results by topic (`folder()`, `results()`), the
  U-Chrom checkout (`UCHROM_REPO`), the data directory (`DATA`).
* `_record.py`: results with their provenance (`save`, `load`, `compare`; `python shared/_record.py OLD NEW --rtol x`
  compares a new run with a committed one).
* `tests/`: `pytest shared/tests`.
