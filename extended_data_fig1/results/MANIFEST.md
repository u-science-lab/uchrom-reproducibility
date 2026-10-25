# extended_data_fig1/results -- results committed here

Most store results are the laptop series (this folder,
the laptop series: Apple M5 Pro).  Runs on another machine are kept here, with `shared/_record.py` provenance,
so that they never merge into the laptop series.

| File | What | Script | Produced |
| --- | --- | --- | --- |
| `c_access_sherlock.json` | single-cell / trace / chromosome reads at 1e6 and 1e7 spots (Takei 2025 cerebellum rep 1 whole-cell subsets, files on node-local disk): chromdata.zarr backed, AnnData read whole (`read_h5ad`), AnnData read lazily (`anndata.experimental.read_lazy` on AnnData Zarr; each query split into the mask and `to_memory`), SpatialData (lazy points + filter); 20 queries x 3 fresh processes | `c_access.py --out-json` via `sherlock/access.sbatch` | 2026-10-13, Sherlock job 47636698 (node sh03-01n43, AMD EPYC 7502, 13 CPUs, 96 GB), revision ember-deep-frozen-216-g5b276e5d114a; library versions in the file (numpy 2.2.6 / pandas 2.3.2 / pyarrow 20.0.0: the nodes' glibc 2.17 has no wheels of the laptop's) |
| `c_access_laptop_zarr23.json` | the chromdata.zarr backed series of Extended Data Fig. 1d measured again on the laptop of the paper's series (Apple M5 Pro; Python 3.11.15, numpy 2.4.3, pandas 2.3.3, pyarrow 23.0.1, zarr 3.1.6, as the first run): 1e4-1e8 spots, ids drawn among the cells and traces with spots; 1e8 = the 1e7 store replicated (`replicate_store.py`, fixed for cell codes without spots).  Merged into the paper's `c_access.json` as series `chromdata.zarr 2.3 (backed, from disk)`, which replaces the 2.2 series in the figure | `c_access.py --backends zarr-backed --out-json` | 2026-10-13, laptop, revision e7f7ed8 + the `replicate_store.py` fix of this commit |
| `anndata_zarr_chunks_sherlock.json` | chunk shapes of the AnnData Zarr stores read lazily (what a one-cell read touches) | `sherlock/access.sbatch` | the same job |
| `clusters_sprite.json`, `clusters_sprite_2M.json` | multi-way contacts in the store (`ChromData.clusters`, format 2.4, U-Chrom PR #153) on Quinodoz 2018 SPRITE mouse ES cells (GSE114242, `mouse.mapq-ge10.clusters.gz`: 188.7 M clusters, 406.3 M reads; `_2M`: its first 2 M clusters): round trip (checksum over clusters: equal), store size by part, parse / build / write / open times, region queries of 40 kb / 1 Mb / 10 Mb windows in memory, backed and over HTTP (local server: requests, bytes) | `g_clusters_sprite.py` via `sherlock/clusters.sbatch` | 2026-10-24, Sherlock job 49446392 (8 CPUs, 400 GB requested; peak RSS 100.6 GB), U-Chrom revision df21e7440ab5; outputs archived on Oak `uchrom/results/clusters_sprite/2026-10-24_df21e7440ab5/` |

Superseded, not committed: job 47631173 (the same run before the chromdata.zarr backend drew its ids only from cells
with spots, and with Linux peak RSS from `ru_maxrss`, which reported the staging parent's peak in every child).
