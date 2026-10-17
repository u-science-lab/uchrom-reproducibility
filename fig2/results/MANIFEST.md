# fig2/results

| File | What |
|---|---|
| `cross_assay.json`, `.md` | Steps 1-2 (`cross_assay.py`): per dataset the class medians, ExN - Oligo and the class ranking with every per-cell value from `uchrom.fea.separation_profile`; Kendall's W |
| `nucdyn/` | The NucDynamics reconstructions (`recon.py`, `recon_summary.py`, `transfer.py`, `biology.py`, `compartments.py`); below |
| `cross_assay_ed_*`, `cross_assay_cells.parquet`, `cross_assay_profiles.parquet`, `cross_assay_figure.json` | Fig. 2's inputs (`figdata.py`, `structures.py`) |

## Run 2026-10-11, Steps 1-2 (`cross_assay.py`)

- **Jobs**: 47401395 (`sherlock/liu2025.sbatch`, 2 min, `sh04-04n05`): `fetch liu2025_mop` (44 s) and the atlas recipe
  `build_liu2025_mop.py` -> `$SCRATCH/uchrom-nd/data/liu2025_mop/` (9,582,069 spots, 188,891 traces, 11,045 cells;
  md5 of both files equal to the registry's, printed by job 47403340).  47403340 (`sherlock/cross_assay.sbatch`,
  8 min 14 s, `sh04-04n05`, 8 CPUs, 10.4 GB peak, `venv-consensus` as above): `seq` 108 s, `tan` 134 s (pairs ->
  `.scool`), `imaging` 117 s (Liu 51 s), `check` 38 s, `analyze` 31 s.
- **Code**: the worktree at ad8aa0f with uncommitted `packages/uchrom/uchrom/fea/separation.py` (md5 `9782acd1`),
  `cross_assay.py` (`10b25b4e`) and `pilot.py` (`bd31a485`), synced through the DTN to `$SCRATCH/uchrom-nd/repo-xa`;
  `provenance.revision` records the same.
- **Check against the pilot** (`stage_logs.check`): the sequencing counts per cell (`n_cis`, `n_ge1`, `n_short`,
  `n_long`, `n_mito`) are equal for every cell of the seven atlas datasets; Tan 2021 differs by at most 3 contacts a
  cell (`pairs_to_scool` bins the 1-based positions as `(pos - 1) // 1 Mb`, the pilot `pos // 1 Mb`), its ExN - Oligo
  is unchanged (+1.449); the imaging pair counts are equal and the medians agree to 1.3e-7 (float32 distances in the
  pilot) with the pilot's closed ranges.
- **Archive**: `$OAK/uchrom/results/cross_assay/2026-10-11_full_47403340/cross_assay_work_results.tar.gz` (97 MB:
  per-cell tables, profiles at 1-Mb offsets, the Tan `.scool`).
- **Figure inputs** (stage `figure`, jobs 47405650 and 47406658, 1 min each, `ARCHIVE=0`): the Extended Data
  figure's `cross_assay_{figure.json,cells.parquet,profiles.parquet}`, written to `$SCRATCH/uchrom-nd/runs/xassay/full/figdata`
  and copied through the DTN to `results/`.  Job 47406658 (`cross_assay.py`
  md5 `b399354d`) added `liu_lineage_pooled` (Liu's OPCs pooled with its oligodendrocytes: 13.8 % OPC; Astro against
  the lineage, Cliff's delta 0.16 [0.12, 0.20]); the per-dataset rows are identical to job 47405650's.
- **Committed code**: `cross_assay.py` as committed is md5 `b399354d` (job 47406658); job 47403340 ran `10b25b4e`, which
  differs only by the later `figure` stage (the functions `_cliff`, `_cliff_ci`, `_per_cell`, `_figdata`,
  `stage_figure`, the option `--figdata`); `fea/separation.py` (`9782acd1`) and `pilot.py` (`bd31a485`) are committed
  as run.

## Figure inputs

- Figure inputs: `figdata.py` (local) -> `results/cross_assay_ed_*`, read by the Fig. 2 script of the paper repository.
- The locus panel of the local runs is `$SCRATCH/uchrom-nd/runs/xassay/liu2025_panel.parquet` (the Liu 2025 store's
  bins); the per-cell scores `$SCRATCH/uchrom-nd/runs/xassay/transfer/scores_*.parquet`.
- Example structures (Fig. 2): `structures.py`, job 47525959 (the first try, 47525588, required 80 %
  of chromosome 2's loci in an imaged copy and found none: detection is about half per copy) ->
  `$SCRATCH/uchrom-nd/runs/xassay/structures/cross_assay_ed_structures.parquet`, copied through the DTN to `results/`.  Cells: imaging ExN 103 / Oligo 99 of 148 loci; Dip-C cortex-p309-bc_153 /
  cortex-p056-cb_054; scHiCAR CCCTACATCGATGGAAGA / CCCTACTTCTCCAACAGA; each at its class median score.
- The locus panel of the figure (`results/cross_assay_ed_panel.parquet`: the Liu 2025 store's bins, chrom /
  start / end) is a copy of `$SCRATCH/uchrom-nd/runs/xassay/liu2025_panel.parquet`.

## Run 2026-10-12/13, NucDynamics reconstructions (`results/nucdyn/`; the paper's Fig. 2)

The reconstructions of the NucDynamics 2017 protocol of the native engine (phased contacts as one chain per
chromosome copy, contacts with an unphased end as ambiguous restraints), analysed as in the README ("Distance
space"); Fig. 2 uses this run (`figdata.py`).

| File | What | Run |
|---|---|---|
| `nucdyn/recon/recon_<dataset>.json`, `loci_<dataset>.json`, `loci_all_<dataset>.json`, `contacts_<dataset>.json` | Stages analyze, loci, loci_all, contacts as above | reconstructions: GPU arrays 47565733 (tan2021), 47565734 (hires_brain), 47565747 (tan2021_maps), 47565752 (gageseq), 47565754 (dschic), 47565755 (droplet), 47565756 (schicar), 40 tasks each (`recon.sbatch`, `METHOD=nucdyn`); collect + stages: jobs 47634791, 47634797, 47634801, 47634803, 47634804, 47634806, 47634808 (`recon_collect.sbatch`) |
| `nucdyn/recon_summary.{json,md}` | `recon_summary.py` over `nucdyn/recon/` | local |
| `nucdyn/transfer.json` | `transfer.py` with the NucDynamics signatures | job 47634836 |
| `nucdyn/compartments.json` | `compartments.py` against the NucDynamics consensus | job 47634837 |
| `nucdyn/biology.json` | `biology.py` on the NucDynamics locus tables | local |

- Code: `recon.py` etc. at 78e7994 with the `--method` and contacts changes uncommitted at run time (provenance
  `78e7994+nucdyn-contacts-uncommitted`), committed here unchanged: md5 of the run's repo copy (`repo-xa`) and of
  the commit agree (`recon.py` `94d2b9de`, `transfer.py` `a1a0dde0`, `compartments.py` `18408f74`, `biology.py`
  `25e43716`).  The GPU arrays ran on H100 nodes (one GPU, 4 CPUs per task).
- Archives: `$OAK/uchrom/results/cross_assay/2026-10-13_recon_nucdyn_<dataset>_<job>/recon_<dataset>.tar.gz`
  (76-222 MB each: models, logs, results, cell and locus tables).
- An earlier CPU attempt for Droplet Hi-C (array 47558438, 14 tasks) is superseded by the GPU arrays.
