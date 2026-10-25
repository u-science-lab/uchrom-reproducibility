# Version 3, family `bulk`: report

Seven Part B questions (`../questions/bulk.yaml`), reference values in `../results/bulk.json` (`truth/bulk.py`; data
in `../data/bulk.md`).  Every question was computed with U-Chrom (the agents' environment `venv-uchrom`, snapshot
fb43ac4) and with an independent implementation that does not import `uchrom` (checked at the end of each run):
the authors' / Open2C's code in `venv-hicref` (`$HOME/uchrom/agentbench/venv-hicref`, built inside
job 47730778 in 281 s: mustache-hic 1.3.3, hic-straw 1.3.1, cooltools 0.7.1, cooler 0.10.4, scikit-learn 1.9.1,
numpy 1.26.4, pandas 2.3.3, scipy 1.14.1 — freeze in `venv-hicref.freeze.txt`), or numpy written from the papers.

cluster job (16 CPUs; logs `$SCRATCH/uchrom/agentbench/logs/`): `ab3-bulk-setup` 47730778
(failed: a cooler balance statistic is an array; venv-hicref built), 47731273 (failed after the 462-MB Zenodo download:
no pyarrow in `venv-verify`; Rao / Su folders built), 47734059 (the RNA table, in `venv-uchrom`, 23 s);
`ab3-bulk-truth` 47734061 (all seven questions, both implementations, 39 min, 51 GB peak), 47740189 (P1 again after
the definition fix below, 17 min; U-Chrom's P1 numbers reproduced to the last digit).

Runtimes below are wall-clock seconds on one cluster node (16 cores, `OMP_NUM_THREADS=16`); "independent" is what
an agent following the prompt with the named public code would spend computing (the environment build, ~5 min, comes
on top for H1-H3).

## Summary

| ID | quantity | U-Chrom | independent | agreement | grade | decision |
| --- | --- | --- | --- | --- | --- | --- |
| H1 | Mustache loops, IMR90 chr21 10 kb | 157 | 165 (mustache-hic 1.3.3 CLI) | 4.9 % apart; 91.5 % of the reference's loops within 2 bins of ours | 161 ± 4 % (154.6-167.4) | keep with tolerance |
| H2 | HiCCUPS recall / precision vs Rao 2014 loops | 0.761 / 0.745 (94 loops) | 0.750 / 0.697 (cooltools dots, 99 loops) | both above the floors | floors 0.70 / 0.64 | keep (floor) |
| H3 | APA of the 52 published loops (O/E, mean) | 4.563482789397527 | 4.563482789397528 (cooltools pileup) | identical (1e-15) | ± 2 % | keep |
| H4 | MDS of Su 2020 Hi-C vs imaging medians (Spearman) | 0.8741 | 0.8732 (scikit-learn MDS, 3 seeds ± 0.0003) | 0.001 apart | floor 0.82 | keep (floor), with a caveat |
| XM1 | imaging contact frequency (< 500 nm) vs Hi-C (Spearman) | 0.85711 / 0.93429 | 0.85711 / 0.93429 (numpy) | 3e-8 | ± 2 % | keep |
| P1 | R² of short/long ratio on 5 copy-number PCs, MSR null (bands of 20) | R² 0.343917; null median 0.0337, q95 0.0857, P 1/1001 (10 runs) | R² 0.343917; same null within Monte Carlo (numpy, 10 runs) | R² 1e-16; null runs interleave | ± 10 % (Monte Carlo spread ≤ 6.5 %) | keep with tolerance |
| P2 | partial Spearman total RNA vs median trace Rg | 0.68734 (Spearman 0.82665) | 0.68734 / 0.82665 (numpy) | 4e-16 | ± 2 % | keep |

## H1 — Mustache loops (Roayaei Ardakany et al. 2020) on Rao 2014 IMR90 chr21, 10 kb

*Prompt*: the balanced 10-kb `.cool`; Mustache with sigma0 1.6, 2 octaves, 10 scales, FDR (`-pt`) 0.1, sparsity
(`-st`) 0.88, anchors at most 2 Mb apart; the number of loops.

- U-Chrom (`uchrom.strc.loop.call_loops_mustache`, `MustacheParams(fdr=0.1, sparsity=0.88, max_dist=2 Mb)`): **157**
  loops, 2.9 s.
- mustache-hic 1.3.3, its command line (`mustache -f IMR90_chr21_10kb_hg19.cool -r 10000 -ch chr21 -pt 0.1 -st 0.88
  -sz 1.6 -oc 2 -i 10 -d 2000000 -p 8`): **165** loops, 11.4 s.  Same numbers as `hic_structures/mustache_validate.json`
  (157 / 165) — the map is built by the same recipe.
- Agreement: 86 % of the reference's loops have a U-Chrom loop at the same summit pixel, 91.5 % within 2 bins; 96 % of
  U-Chrom's loops are within 2 bins of a reference loop; recall of the published HiCCUPS loops (anchors within 25 kb)
  0.717 vs 0.728.  The difference is the documented one: U-Chrom fits the exponential null and the BH correction once
  per chromosome, the reference per 2,000-bin block (3 blocks on chr21), which moves borderline candidates across the
  FDR threshold.
- Grade: the reference value is the midpoint, 161, with a relative tolerance of 0.04 (154.6-167.4): it covers both
  implementations plus ~1 % of slack, and stays far from what the command line's own default FDR gives (`-pt 0.2`:
  268 loops with mustache-hic, 264 with U-Chrom) — an agent that drops the stated FDR fails.
- Feasibility: mustache-hic 1.3.3 installs on Python 3.12 only with numpy < 2 and scipy < 1.15 (it imports
  `scipy.ndimage.filters`) and builds hic-straw from source (works in the python:3.12-bookworm image); `pip install
  mustache-hic` alone pulls numpy 2 and fails at import, so the agent has to pin.  Runtime 11 s.
- **Decision: keep with tolerance** (band covering both).

## H2 — HiCCUPS (Rao et al. 2014) recall of the published chr21 loops

*Prompt*: HiCCUPS on the 10-kb map with p = 2, w = 5, per-diagonal expected (two diagonals ignored), λ-chunked BH FDR
0.1 per kernel, all four kernels, ≤ 10 Mb, 20-kb clustering, the Rao 2014 filter; a published loop (92 on chr21) is
recovered when a call has both anchor midpoints within 25 kb; report the number of loops, recall and precision.

| implementation | loops | recall | precision | time |
| --- | --- | --- | --- | --- |
| U-Chrom `call_loops_hiccups` (HiCCUPS p-value P(X ≥ obs), greedy HiCCUPS clustering) | 94 | 0.761 (70/92) | 0.745 | 1.7 s |
| cooltools 0.7.1 `dots` (defaults at 10 kb: w 5, p 2; Birch clustering 20 kb) | 99 | 0.750 (69/92) | 0.697 | 5.2 s |
| U-Chrom, loops ≤ 2 Mb | 92 | 0.750 | 0.750 | 0.6 s |
| cooltools, ≤ 2 Mb | 96 | 0.739 | 0.708 | 2.0 s |

- Floors (numbers_min): recall ≥ **0.70**, precision ≥ **0.64** = the smaller of the two implementations minus 0.05
  (about five loops of 92), rounded down: the two formulations differ in the p-value convention and the clustering
  (U-Chrom's module docstring), and an agent may clip at 2 Mb; both implementations pass in both settings.  Precision
  is graded too so that calling everything does not pass.  The published list was made from merged 5/10/25-kb calls,
  so a single 10-kb map cannot reach recall 1.
- Feasibility: cooltools 0.7.1 (the current release) `dots` fails on pandas 3 (`benchmarks/hic_structures/README.md`);
  an agent must pin pandas < 3.  Runtime 5 s.
- **Decision: keep (floor).**

## H3 — APA (Rao et al. 2014) of the published chr21 loops

*Prompt*: loops at the 10-kb bins of their anchor midpoints, j - i ≥ 22; 21 x 21 windows of balanced O/E (expected =
per-diagonal mean over valid pixels); windows averaged pixel by pixel (NaN ignored); APA = centre / mean of the
lower-left 6 x 6 corner (rows i+5..i+10, columns j-10..j-5).

- U-Chrom `uchrom.fea.pileup(..., flank=100 kb, aggregate="mean", corner=6, min_separation=22)`: 52 loops,
  **APA 4.563482789397527**, 0.13 s.  (The pooled-ratio aggregate, U-Chrom's default, gives 4.215: the prompt names the
  mean of the O/E windows.)
- cooltools 0.7.1 `pileup` with `expected_cis(ignore_diags=2)` 'balanced.avg' and `np.nanmean` over the stack: 52
  loops, **4.563482789397528**, 1.0 s.  Identical to 1e-15 (same as `hic_structures/results/pileup.json`, depth 100 %).
- Grade: numbers [n_loops, apa], tol 0.02.  The definition is fully stated; the likely agent errors (pooled sum
  instead of the mean of windows: 4.21; observed instead of O/E; another corner) fall outside ± 2 %.
- **Decision: keep.**

## H4 — bulk 3-D reconstruction (metric MDS) of Su 2020 chr21 from Hi-C, against imaging

*Prompt*: Hi-C_contacts_chromosome21.tsv (651 loci, raw counts); diagonal ignored; distance-decay prior c' = 0.95 c +
0.05 E(|i-j|) (E = mean of the non-zero counts at each index separation; 9.5 % of the pairs are 0, every separation
has non-zero pairs, so every c' > 0); d = c'^(-1/4); metric MDS (raw stress, SMACOF) to convergence; Spearman of the
structure's distances with the imaging median distances (copies with both loci detected) over all 211,575 pairs.

- U-Chrom (`uchrom.io.read_su2020_hic` → `write_locus_cool` → `uchrom.recon.bulk.mds.reconstruct_mds(balance=False,
  alpha=4, weight=0.05, method="smacof")`, imaging medians `uchrom.fea.distance_map`): **ρ = 0.8741**, 1.6 s for the
  reconstruction (176 s with reading the 260-MB table and the medians).
- Independent: the prior and distances in numpy (equal to U-Chrom's `apply_distance_decay_prior` + `contact_to_distance`
  to 2e-16, checked locally), scikit-learn 1.9.1 `MDS(metric="precomputed", init="random", n_init=4, max_iter=3000,
  eps=1e-9)` with seeds 0/1/2: **ρ = 0.8731, 0.8731, 0.8734** (19 s each); imaging medians in numpy (≈ 2 min).
- Floor: **0.82** = min(0.8741, 0.8732) - 0.05.  Sensitivity (U-Chrom): α = 3 gives 0.896; the ICE-balanced map 0.870;
  zero pairs left out of the stress instead of the prior 0.599 (with U-Chrom's diagonal SMACOF update) — so the floor
  separates the stated handling of zeros from leaving them out.
- Caveat: the target distances themselves, without any embedding, already correlate 0.862 with the imaging medians;
  the floor therefore checks the input preparation, not that a structure was built.  A floor above 0.862 would sit
  within 0.01 of both implementations, too tight for other solvers; the prompt asks for the structure's distances, and
  the judge can check that a structure was reconstructed.
- **Decision: keep (floor 0.82)**, caveat recorded.

## XM1 — imaging contact frequency vs Hi-C on Su 2020 chr21

*Prompt*: per pair i < j, the fraction of copies (with both loci detected) closer than 500 nm; Spearman with the Hi-C
counts over all 211,575 pairs and over the 24,815 pairs whose midpoints are ≤ 2 Mb apart.

- U-Chrom (`uchrom.fea.distance.contact_frequency(cd, 0.5)` on `read_su2020_tracing`, µm): **0.857108 / 0.934286**,
  96 s.
- numpy (nm, threshold 500): **0.857108 / 0.934286**, 121 s; differences 3e-8 (distances within float rounding of
  500 nm in nm vs µm).
- Grade: numbers [rho_all, rho_le_2mb], tol 0.02.  Likely agent errors that the tolerance catches: Pearson instead of
  Spearman, ≤ instead of < (negligible), counting undetected loci in the denominator, dropping zero Hi-C pairs.
- **Decision: keep.**

## P1 — Moran spectral randomisation (Wagner & Dray 2015) for the melanoma copy-number share

*Prompt*: Guo 2026 melanoma, section human_melanoma_r1 (6,268 spots); y = log2 short / long contacts (bin offset 1 vs
offsets 19-53, + 0.5, autosomes); copy-number profile = log2 of each chromosome's share of contact legs (chr1-22, X)
over its median share; five PCs of the z-scored profiles; observed R² of y on the PCs; null = MSR of the PCs on the
8-neighbour spot grid (doubly centred binary weights, eigenvectors orthogonal to 1, increasing eigenvalue, bands of 20,
random Gaussian directions of the same length per band and component), 1,000 draws; report the observed R², the null
median and 95th percentile, P.

- U-Chrom path: y = `uchrom.fea.separation_profile` on the embedded per-spot map; legs from `chromdata.embedded.export_scool`;
  spot graph `uchrom.spatial.spatial_neighbors(method="radius", radius=1.5)` on (array_col,
  array_row) (48,724 links); the null with the benchmark's own MSR code (`_msr_basis`, `msr` of
  `benchmarks/spatial_variable/spatial_variable.py`, PR #127 / 52c3f77, copied into `truth/bulk.py` as
  `pr127_msr_basis` / `pr127_msr`; the basis by its shift of the constant).  **R² = 0.343917**; 10 runs of 1,000 draws
  (seeds 0-9): null median 0.0315-0.0357, q95 0.0824-0.0899, P = 1/1001 in every run.  820 s (654 s reading the maps
  through `export_scool`, 11 s basis, 15 s per 1,000 draws).
- Independent numpy: y and legs in one pass over the embedded map's partitions (chromdata's reader, no coolers), PCs
  by SVD, the basis from an orthonormal complement of the constant (QR of [1 | I], eigenvectors of Q'WQ), bands, OLS.
  **R² = 0.343917** (equal to 1e-16), 10 runs (seeds 10,000-10,009): null median 0.0319-0.0358, q95 0.0827-0.0891,
  P = 1/1001.  153 s in all (30 s reading, 44 s basis, 8 s per 1,000 draws): an agent needs ~2 min of compute.
- Graded values (`results/bulk.json`): R² 0.343917; null median **0.03367** and q95 **0.08569** = the means of the 20
  runs; P = 0.000999.  Monte Carlo spread over the 20 runs: median s.d. 3.8 % (largest deviation 6.5 %), q95 s.d.
  2.5 % (5.0 %); the seed is the agent's, so the tolerance is **0.10** relative (2.6 s.d. of the median, 4 of the q95;
  per-number tolerances r2 0.02 / median 0.10 / q95 0.08 / P 0.10 in `tol_by_number` for a grader that takes them).
- A definition fix during the build: the first run computed y in the independent path from bin midpoints (the shorter
  last bin of a chromosome moves a few pairs out of the ranges), U-Chrom from bin offsets x 1 Mb; R² 0.34355 vs
  0.34392, median y 0.2158 vs 0.2217.  The prompt now states bin offsets, and both agree.
- Finding: copy number explains 34 % of the spot-to-spot variance of the short / long ratio; random smooth covariates
  with the same spatial spectrum explain 3.4 % (95th percentile 8.6 %), P = 0.001 — the copy-number share is not a
  generic effect of spatial smoothness (no RNA or depth adjustment here, unlike stage cnv_null of PR #127).
- **Decision: keep with tolerance.**

## P2 — partial Spearman, total RNA vs median trace Rg, Takei 2025 FOV 0 (cell as unit)

*Prompt*: traces with ≥ 50 spots (1,393 of 1,500), Rg = RMS distance of all spots from their centroid, cell value =
median over its traces; RNA = total transcripts of the 60-gene table; Spearman and partial Spearman given cell type
(ranks minus the type mean of the ranks, Pearson), 47 cells, 6 types.

- U-Chrom (`uchrom.fea.distance.radius_of_gyration`; RNA totals from `uchrom.io.read_seqfish_transcriptome` on the
  Zenodo zip, identical to the question's table; scipy ranks + least squares on the type indicators): **partial
  0.687337**, **Spearman 0.826645**, 2.6 s.
- numpy / pandas (Rg from the coordinates, the table, ranks minus group means): **0.687337 / 0.826645**, 0.4 s;
  difference 4e-16.
- Grade: numbers [spearman, partial_spearman], tol 0.02 (± 0.017 / 0.014).  The answer moves beyond the tolerance
  with every choice the prompt fixes (measured on the same data, Spearman / partial): traces with ≥ 2 spots 0.852 /
  0.717, ≥ 100 spots 0.761 / 0.631, mean instead of median over traces 0.856 / 0.748, autosomes only 0.774 / 0.605,
  Pearson on raw values 0.597 / 0.408 (partial on group-centred values), no cell-type control 0.827 vs 0.687.
- Note: the store itself holds no RNA counts; the table is built from the registered `takei2025_rna` archive (its only
  cerebellum mRNA table is the 60-gene non-barcoded panel, mostly cell-type markers — hence the strong raw
  correlation, mostly carried by cell type).
- **Decision: keep.**

## Open issues

- IDs H1 / H2 collide with Part A's H1 / H2 (HiRES atlas): rename at the merge.
- P1's grade uses one relative tolerance for all four numbers (the grader takes a scalar `tol`); the deterministic
  observed R² would deserve 0.02 (see `tol_by_number` in the YAML, for a grader that supports it).
