# Supplementary Fig. 6: the browser's response times and the atlas

Two benchmarks, each writing its own `results/`:

* `browser/`: the web browser's response times on Takei 2025 FOV 0 replicated 1-30x (`make_scaled.py` builds the
  inputs into `browser/_data/`, not tracked; `bench_server.py`, `bench_frontend.py`; `run_all.sh` runs them on an
  idle machine) -> `browser/results/server.json`, `frontend.json`.  Panel a.
* `atlas/`: the atlas datasets checked against the authors' published numbers (`validate.py` -> `atlas.json`;
  `hic_embeddings.py`, `spatial_hicrna_validation.py`, `spatial_hic_chen2026_reprocessed.py`) and the hosted
  browser reading the atlas stores over HTTP (`hosted.py` -> `hosted.json`).  Panels a-c.

`atlas/results` is also used by Supplementary Table 2 (the atlas table), Fig. 1 (the atlas counts) and
Supplementary Note 3.
