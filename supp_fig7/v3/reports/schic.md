# Family schic (version 3, Part B): reports

## Summary

| ID | Method (paper) | Graded quantity | U-Chrom | Independent | Agreement | Authors' runtime (8 CPUs) | Decision |
|---|---|---|---|---|---|---|---|
| R1 | NucDynamics (Stevens 2017), Cell 1 | Spearman vs published model 1, chr1, 1 Mb | 0.9890 | nuc_dynamics @ 00eaf0b: 0.9868 (default anneal, 2 seeds), 0.9855 / 0.9858 (5x / 20x shorter) | all above the floor 0.939 | 1,887-1,967 s default; 387 / 134 s shorter | keep, floor (prompt: anneal length free) |
| R2 | NucDynamics, Cell 2 | same | 0.9656 | nuc_dynamics: 0.9617 (default anneal), 0.9633 / 0.9620 (5x / 20x shorter) | all above the floor 0.915 | 1,496 s default; 305 / 92 s shorter | keep, floor |
| R3 | RWR imputation (scHiCluster, Zhou 2019) | Pearson of two cells' imputed chr1 maps | 0.694665 | impute_cpu (schicluster 1.3.2): 0.694666 | 8e-7 | 0.04 s | keep, parity, tol 0.005 |
| R4 | scHiCluster embedding | kNN purity of cell types (k = 30) | 0.2433 (fb43ac4); 0.240191 after the fix | impute_cpu + hicluster_cpu: 0.240191 (variants 0.236-0.247; no imputation 0.213) | fixed U-Chrom identical; both above the floor 0.230 (b6d8873) | 289 s | keep, floor (k changed 10 -> 30) |
| R5 | scHiCluster compartment score (scA/B) | partial r with the pseudo-bulk E1 given CpG | 0.26951 | single_chrom_compartment v1.3.4 + cooltools: 0.26951 | 1e-10 | 76 s | keep (reworded), parity, tol 0.05 |
| R6 | cell-cycle phasing (Nagano 2017) | median near / mitotic fraction, 5 phase counts | 0.789556 / 0.093392; 1/38/125/214/17 | numpy and the authors' published values: identical | exact (per cell 2e-5; published 5e-10) | 63 s (numpy) | keep, parity, tol 0.02 |
| R7 | SnapHiC (Yu 2021), 50 cells, chr19 | number of loops | 17 | SnapHiC @ 782dc2b: fails on Python 3.12 (pandas 3: interaction step; pandas 2.3: post-processing) | none | > 833 s (failed) | drop |

## How the references were computed

- **U-Chrom** (`truth/schic.py`, the agents' environment `venv-uchrom`, commit fb43ac4, image py312full.sif,
  `sherlock/schic_truth.sbatch`, 16 CPUs).
- **Independent implementations**, each in its own environment without U-Chrom, built inside a job the way an agent
  would build it on Python 3.12 / Debian bookworm (`sherlock/schic_setup.sbatch`), and run on 8 CPUs (an agent's
  session):
  - `venv-schicluster`: `pip install schicluster==1.3.2` (PyPI's newest) **fails** (its dependency rpy2 needs R:
    "rpy2 in API mode cannot be built without R"); `pip install --no-deps schicluster==1.3.2` + numpy, scipy,
    scikit-learn, pandas, cooler, h5py, xarray, torch (CPU) works (389 s, torch is 200 MB; `schicluster.schicluster`
    imports torch at module level).  The GitHub tag v1.3.4 (a248677) is checked out for the compartment score (PyPI's
    1.3.2 lacks its row normalisation).  `truth/schic_ref_schicluster.py`.
  - `venv-nucdyn`: `tjs23/nuc_dynamics` @ 00eaf0b; its Cython module does **not** compile with the newest numpy 2.5 /
    Cython 3.3 (`dyn_util.pyx`, CompileError); with `numpy<2` and `cython<3` (numpy 1.26.4, Cython 0.29.37, wheels for
    3.12) it builds in 48 s.  Its `open_file` opens plain text with mode `'rU'`, which Python 3.11+ rejects
    (`ValueError: invalid mode: 'rU'`): the NCC input has to be gzipped (gzip files are opened with `'rt'`).
    `truth/schic_ref_nucdyn.py`.
  - `venv-verify` (cooltools 0.7.1): the E1 of R5 and the numpy implementation of R6
    (`truth/schic_ref_numpy.py`).
  - `venv-snaphic` (R7): `HuMingLab/SnapHiC` @ 782dc2b (its requirements.txt pins numpy 1.19 / pandas 0.24: no
    wheels for 3.12).  With the newest packages (pandas 3.0) it fails in the interaction step
    (`interaction_caller.py:253`, `TypeError: Invalid value ... for dtype 'int32'`) after 237 s; with `pandas<3`
    (2.3.3) in the post-processing (`'eta_cluster' is both an index level and a column label`) after 833 s.
- `truth/schic_collect.py` merges everything into `results/schic.json` (format of `answers_v3.json`).

## R1 / R2: NucDynamics reconstruction of Stevens et al. 2017 Cells 1 and 2 vs the published structures

**Prompt.**  R1 on Cell 1; R2 the same on Cell 2 (GSM2219498, 79,569 contacts).  Stated precisely:
particle sizes 8, 4, 2, 1 Mb as an example, the bin rule of the published particles ([k, k+1) Mb by the PDB's last
column) and the pair set (each unordered pair of distinct chr1 bins present in both, once).  Spearman over the 18,528
chr1 bin pairs (193 bins).  Floor: min(rho of U-Chrom's runs with seeds 0 and 1, rho of the published model 1
vs model 2) - 0.05.

| | Cell 1 (R1) | Cell 2 (R2) |
|---|---|---|
| U-Chrom `reconstruct_nucdyn` (protocol nuc_dynamics_2017, CPU, 8/4/2/1 Mb), seeds 0 / 1 | 0.98908 / 0.98895 | 0.96536 / 0.96590 |
| published model 1 vs model 2 | 0.99943 | 0.99776 |
| **floor** | **0.939** | **0.915** |
| authors' nuc_dynamics @ 00eaf0b, `-m 1 -s 8 4 2 1 -ran <seed>` (default anneal: 500 temperatures x 100 steps per stage), seeds 0 / 1 | 0.98679 / 0.98684 | 0.96170 (seed 0) |
| authors' nuc_dynamics, shorter anneal `-temps 100` (seed 0; 5x fewer steps) | 0.98554 | 0.96334 |
| authors' nuc_dynamics, shorter anneal `-temps 50 -dyns 50` (seed 0; 20x fewer steps) | 0.98583 | 0.96200 |
| U-Chrom runtime per run (16 CPUs) | 150 s | 122 s |
| authors' runtime per run (8 CPUs; the code uses one core per model) | default 1,887 / 1,967 s; `-temps 100` 387 s; `-temps 50 -dyns 50` 134 s | default 1,496 s; `-temps 100` 305 s; `-temps 50 -dyns 50` 92 s |

**Decision: keep both, floor (numbers_min on rho; R1 0.939, R2 0.915).**  The authors' code reaches both floors, even
with an anneal 20 times shorter than its default (0.986 and 0.962 in 134 s and 92 s), so the floors reject failed or
wrong reconstructions, not other correct ones.  The default anneal of the authors' code is far slower (see the table)
than an agent's 1,800-s session allows for a first attempt plus a fix, so the prompts now say that the length of the
annealing schedule is up to the agent.  Getting the authors' code to run at all costs an agent two fixes on Python 3.12
(`numpy<2 cython<3` to compile `dyn_util.pyx`; a gzipped NCC file because of the `'rU'` file mode) plus writing the NCC
file from the GEO table.  Cell 2 (79,569 contacts) is reconstructed less faithfully than Cell 1 by both implementations
(0.966 / 0.962 vs 0.989 / 0.987), which is why its floor is lower.

## R3: random walk with restart imputation of scHiCluster, two scHiCAR cells, chr1 at 1 Mb

**Prompt.**  Cells GTGAAATCTACGAGATAC (Oligo) and TAGACCGGTATGTACCCG (L6CT): the cells with the median chr1 contacts of
their types (1,013 and 1,234 contacts; `truth/schic_explore.py`).  The imputation exactly as the authors' code
(`schicluster.schicluster.impute_cpu`, PyPI 1.3.2): A = log2(U + U^T + 1) of the stored upper triangle (the diagonal
doubled, as the code does it), 3 x 3 mean filter with zeros outside, diagonal set to 0 then 1 for empty columns,
column-normalised P, Q <- 0.5 I + 0.5 Q P from I until the Frobenius change < 1e-6 (<= 30 steps), symmetrised; Pearson
r over the 19,110 pixels i < j.

| | r |
|---|---|
| U-Chrom (`uchrom.fea.impute_rwr.convolve_maps` + `rwr_dense(self_loops=False)`, the code's diagonal rule applied in between) | **0.694665** |
| authors' `impute_cpu` (schicluster 1.3.2) | **0.694666** (relative difference 8e-7) |
| U-Chrom's `uchrom.emb.hic.schicluster_impute` as is, on the counts | 0.6874 (-1.05 %) |
| U-Chrom's `schicluster_impute` as is, on log2(U + U^T + 1) | 0.7025 (+1.13 %) |
| authors' Q not symmetrised (upper triangle of Q) | 0.6702 (-3.5 %) |
| diagonal included | 0.9948 |
| raw counts, no imputation | 0.086 |

**Runtime** of the authors' code: 0.04 s for the two cells (plus the install: `pip install schicluster` fails on rpy2;
`--no-deps` + numpy/scipy/scikit-learn/pandas/cooler/torch, 6.5 min with torch's download).

**Decision: keep, parity, tol 0.005 (relative).**  The two implementations agree to 1e-6; every deviation from the
stated recipe that an implementation might make (self-loops instead of the code's diagonal rule, no log, no
symmetrisation, diagonal pixels) moves r by >= 1 %, so a 0.5 % tolerance separates the recipe from its variants.
Note for the paper: U-Chrom's `schicluster_impute` deviates from the authors' code in the diagonal handling (identity
added to every bin instead of zeroing the diagonal and adding 1 only to empty bins) and leaves the log transform to the
caller; an agent that calls it unchanged misses this question (r 0.687 / 0.702).  U-Chrom's lower-level kernels
reproduce the authors' value exactly.

## R4: scHiCluster embedding of the 5,313 scHiCAR cells, k-nearest-neighbour purity of the cell types

**Prompt.**  1-Mb maps, chr1-19; pad 1, RWR restart 0.5, top 20 % binarisation per cell and chromosome, 20 PCs per
chromosome, concatenation (380), PCA to 20; purity of `cell_type` (22 labels, `data_v3/schicar_cell_types`, identical
to the store's cell table) among each cell's **30** nearest neighbours.  The authors' `hicluster_cpu` sets the PCs per
chromosome itself (int(0.2 x min(cells, n^2)) - 1 = 1,061 here) and runs 200 k-means restarts; the question states 20
PCs (the `ndim=20` of the function's signature, which its loop overrides), and the reference runs the authors'
`impute_cpu` body and `hicluster_cpu`'s binarisation with 20 + 20 PCs.  The details the paper leaves open (log2 of the
counts, the diagonal, binarising the n x n matrix or its upper triangle) are not stated: the floor covers them.

| | k = 30 (graded) | k = 5 | k = 10 | k = 15 | L2/3 IT pooled (k = 10) |
|---|---|---|---|---|---|
| U-Chrom `uchrom.emb.hic.embed_contacts` + `uchrom.emb.metrics.knn_purity` at fb43ac4 (no log, self-loops, upper triangle) | 0.2433 | 0.2591 | 0.2545 | 0.2522 | 0.2868 |
| U-Chrom after the scHiCluster fix (branch `Nanguage/schicluster-parity`: the authors' computation), `n_per_chrom=20` | 0.240191 | 0.2533 | 0.2506 | 0.2472 | |
| authors' `impute_cpu` (schicluster 1.3.2) + `hicluster_cpu` binarisation, PCA 20 + 20 | 0.2402 | 0.2533 | 0.2506 | 0.2472 | 0.2832 |
| authors' code without the log2 (the paper's description) | 0.2468 | 0.2627 | 0.2594 | 0.2554 | 0.2910 |
| authors' code, binarised over the upper triangle (i < j) | 0.2363 | 0.2515 | 0.2467 | 0.2447 | 0.2795 |
| baseline: no imputation (pad 0, no random walk; log2 counts binarised) | **0.2134** | 0.2440 | 0.2347 | 0.2273 | 0.2652 |
| chance (label frequencies) | 0.0718 | | | | |

(The first two rows were rerun with k = 30 as the graded value, jobs 47738050 / 47738057; deterministic: the same
embeddings.  A brute-force numpy kNN on the stored embeddings reproduces every purity above to all digits.)

**Runtime**: the authors' path 289 / 304 s on 8 processes (5,313 cells x 19 chromosomes; two runs); U-Chrom's `embed_contacts` 827 / 864 s
on 16 threads at fb43ac4 (one process reading the cells, dense batched random walks).  After the fix (the authors'
computation, one map per thread with single-threaded BLAS) 325 s on 16 threads, purity identical to the authors' to
all digits (job 47754119); the agents' snapshot will be rebuilt with it.

**Decision: keep, floor at k = 30 = U-Chrom's value minus 0.01 (0.233 at fb43ac4; 0.230 with the fixed code, b6d8873).**  Hi-C at 1 Mb separates these 22 types
weakly (3.4x chance), and the imputation adds little over the binarised raw maps, so the floor has to sit between the
correct implementations and the no-imputation baseline: at k = 30 the four correct variants span 0.236-0.247 and the
baseline is 0.213, at k = 10 the gap (0.247 vs 0.235) was within the variants' spread — hence k = 30 (changed from the
draft's 10; the gap keeps growing with k — U-Chrom / authors / baseline at k = 50: 0.236 / 0.232 / 0.199, at k = 100:
0.222 / 0.221 / 0.175 — but 30 neighbours stay within the smaller types: CLA has 42 cells).  The margin 0.01 keeps
every correct variant above the floor (lowest 0.2363) and the baseline 0.02 below it.  The 0.05 margin of R1's floor rule would have accepted no imputation at all.  With U-Chrom fixed (0.2402, the
authors' value; snapshot b6d8873, job 47799338) the same rule (U-Chrom's value minus 0.01, `schic_collect.py`) gives
the floor **0.230**: every correct variant stays above it (lowest 0.2363), the no-imputation baseline below (0.2134).

## R5: single-cell compartment scores (scHiCluster's scA/B), all scHiCAR cells, vs the pseudo-bulk E1

**Choice of definition and data.**  U-Chrom's `single_cell_compartments` implements two published forms: Dip-C's
(Tan et al. 2021: all contacts, rank-normalised per cell) and scHiCluster's compartment score
(`ScABParams(contacts="cis", normalize="none")`).  Dip-C's code is Python 2.7 (`benchmarks/sc_structures` validated it
with a Python 2.7 environment), which an agent cannot be expected to run on the benchmark's image; scHiCluster's is
Python 3 (`single_chrom_compartment`, GitHub v1.3.4; PyPI's 1.3.2 still multiplies the unnormalised map, i.e. a
contact-weighted sum).  The scHiCAR cells (5,313, 1 Mb, shared data) need no new download; the CpG ratio of the bins is
given (`data_v3/mm10_cpg`, scHiCluster's definition) so that no genome FASTA is needed in a session.

**Why the partial correlation.**  The first version of the question asked for the Pearson r of the mean score with E1:
0.7144 (authors' scHiCluster + cooltools), but the CpG ratio alone correlates 0.7038 with E1 — the score is a smoothed
CpG ratio, and a shortcut that skipped the cells entirely would have passed any tolerance or floor.  The question
therefore asks for the partial correlation given the CpG ratio (the correlation of the residuals of the mean score and
of E1 after a linear fit on the CpG ratio), which is what the contacts add.

| | partial r given CpG | plain r | r(CpG, E1) |
|---|---|---|---|
| U-Chrom (`single_cell_compartments`, cis / none, NaN -> own CpG; `call_compartments_eig` 1 Mb) | **0.2695** | 0.7144 | 0.7038 |
| authors' `single_chrom_compartment` (v1.3.4) + cooltools 0.7.1 `eigs_cis` | **0.2695** | 0.7144 | 0.7038 |

The two agree to 1e-10: the mean scores differ by at most 6.7e-11, the E1s by 1.3e-14.
Sensitivity of the partial correlation (U-Chrom, same E1): bins without a contact in a cell left out of the mean
instead of taking their own CpG ratio 0.2391 (-11 %); median over the cells instead of the mean 0.3142 (+17 %); each
cell's scores rank-normalised first (Dip-C's form) 0.2630 (-2.4 %); the CpG ratio alone 0 (by construction).  The
score itself correlates 0.892 with the CpG ratio.

**Runtime** of the authors' path: 71 s for the 5,313 cells' scores (8 processes) + 5 s for the E1; U-Chrom 10.5 s +
2.8 s.

**Decision: keep (reworded), parity, tol 0.05 (relative) on the partial correlation.**  Parity rather than the floor
the plan suggested: a floor on a correlation with E1 is passed by a shortcut that correlates the CpG ratio itself
(0.704 vs 0.714), and the partial correlation is a deterministic quantity of a fully stated recipe on which the two
implementations agree to 1e-10.  The 5 % tolerance (+-0.013) admits per-cell rank normalisation (-2.4 %) but not the
recipe's other choices (empty bins left out -11 %, median +17 %).

## R6: cell-cycle phasing of Nagano et al. 2017, 395 haploid 2i cells

**Prompt.**  The 395 cells of condition `2i_all` with `passed_qc == 1` (`data_v3/nagano2017_hap_2i`: the authors' adj
tables, `GATC.fends` on mm9, no published features).  Each fend pair once, cis on chr1-19 + X, d = |coord1 - coord2|,
among d > 2^14.5: near = (2^14.5, 2^21], mitotic band = (2^21, 2^23.5]; the schic2 "newer" gates (post-M, pre-M, G1,
early-S, mid-S/G2) stated with their thresholds; report the two medians and the five phase counts.  The definitions are
those of the authors' code (`tanaylab/schic2` @ ab593d2: `analyzeScHiC.r` `sch_cis_breaks` = 2^(10 ... 27.75 by 1/8),
`sch_remove_near_dist` = 2^14.5, `.sch_gen_mat_decay_res` with R's right-closed `cut`; `hap_serum_2i_params.r`
`sch_near_dists` = 2^(14.5, 21), `sch_mitotic_dists` = 2^(21, 23.5) and the gate thresholds; `ordering.r`
`group_cells_by_metrics`).

| | median near | median mitotic | post-M | G1 | early-S | mid-S/G2 | pre-M |
|---|---|---|---|---|---|---|---|
| U-Chrom (`uchrom.io.ContactStats`, 1/8-octave bins from 2^10, + `uchrom.fea.contact_features`, bands by bin) | 0.789556 | 0.093392 | 1 | 38 | 125 | 214 | 17 |
| numpy from the definition (`truth/schic_ref_numpy.py`) | 0.789556 | 0.093392 | 1 | 38 | 125 | 214 | 17 |
| the authors' published per-cell values (GEO GSE94489 feature table, `f_near_band`, `f_mitotic_band`) through the same gates | 0.789556 | 0.093392 | 1 | 38 | 125 | 214 | 17 |

The numpy implementation reproduces the authors' published `f_near_band` and `f_mitotic_band` of every one of the 395
cells to within 5e-10 (the rounding of the published table): the stated definition is exactly the authors'.  (The
published `group` column differs from the gates — late-S/G2 224, early-S 121, G1 38, pre-M 11, post-M 1 — because the
authors then reassign cells by k-means clusters of their decay profiles; the question asks for the gates only.)
Sensitivity: 13 cells lie within 0.005 of the 0.77 gate and 15 within 0.005 of the pre-M line, so implementations that
deviate from the definition (counting the count column, both orientations of a pair, other chromosomes) move several
cells across gates.  Exactly one contact of the 395 cells sits at d = 2^21 (the only edge convention that could matter).

U-Chrom and numpy agree cell by cell to 2.2e-5 (one cell: its single contact at exactly d = 2^21 falls into the
mitotic band under U-Chrom's half-open bins [lo, hi), into the near band under the paper's right-closed ones); the
medians agree to all digits and every count is equal.  **Runtime**: numpy 63 s (reading the
395 adj tables and the 253-MB fends table), U-Chrom 37 s; no authors' code is needed (the published per-cell values are
the authors' output, and the definition is three lines).

**Decision: keep, parity, tol 0.02 (relative).**  Three independent sources agree exactly (U-Chrom, numpy, the authors'
published values).  The tolerance leaves the medians +-2 % and the large phase counts a few cells (mid-S/G2 +-4,
early-S +-2); the small counts (G1 38, pre-M 17, post-M 1) must be exact, which a correct implementation achieves.

## R7 (conditional): SnapHiC loops, 50 Nagano 2017 diploid cells, chr19 at 10 kb

**Prompt.**  50 of SnapHiC's 742 mES cells (seed-0 draw), their chr19 10-kb maps as SnapHiC input (one line per
non-empty pixel at the bin centre), SnapHiC's defaults (whole-chromosome RWR "inverse", rp 0.05, 2-Mb band, candidates
100 kb - 2 Mb, FDR 0.1, t 3, outlier 1.96, the mm10 filter regions), number of loop summits.

| | loops | runtime |
|---|---|---|
| U-Chrom `call_loops_snaphic` (defaults, `exclude` = filter regions; numpy path, uchrom-maps not installed) | 17 | 175 s (8 threads) |
| authors' SnapHiC @ 782dc2b, newest packages (pandas 3.0) | fails in the interaction step | 237 s until the failure (RWR of the 50 cells: 175 s on 8 processes) |
| authors' SnapHiC @ 782dc2b, `pandas<3` (2.3.3) | fails in the post-processing (`ValueError: 'eta_cluster' is both an index level and a column label, which is ambiguous`) | 833 s until the failure: RWR 3 min, local background 10.5 min (8 processes) |

(A second attempt with `pandas<3` first resumed from the previous run's output folder on the node and failed in the
combine step; SnapHiC reuses an existing output folder.  The run reported above started from an empty folder.)

**Decision: drop.**  The authors' SnapHiC does not run to the end on Python 3.12 without patching its code (pandas 3:
dtype error in the interaction step; pandas 2.3: groupby ambiguity in the post-processing; its pinned pandas 0.24 /
numpy 1.19 have no wheels for 3.12), and even its steps before the failure take 14 min on 8 processes for 50 cells
and one small chromosome — beyond the ~10 min an agent's session can give it.  No loop count of the authors' code
exists to set or check a tolerance against U-Chrom's 17 loops.  (U-Chrom's SnapHiC was validated against the authors'
code in `benchmarks/sc_structures` in the authors' own Python 3.9 environment: identical loops with
`sliding_window_holes=True`.)

## What the references showed about U-Chrom (for the paper; no library change made here)

- `uchrom.emb.hic.schicluster_impute` is not the authors' imputation: it adds the identity to every bin (self-loops)
  where the authors' code zeroes the smoothed diagonal and puts 1 only on empty bins, and it leaves the log2 transform
  of the code to the caller.  On R3's two cells this moves the answer by -1.05 % (counts) / +1.13 % (log2 input); the
  lower-level kernels (`convolve_maps`, `rwr_dense(self_loops=False)`) reproduce the authors' value to 1e-6.
- `uchrom.emb.hic.embed_contacts` took 827 s on 16 threads for the 5,313 scHiCAR cells, the authors' code 289 s on 8
  processes (R4).
- `single_cell_compartments(..., ScABParams(contacts="cis", normalize="none"))` gives NaN where a cell has no cis
  contact at a bin; scHiCluster gives the bin's own CpG ratio.  With the NaN filled that way the mean scores over the
  cells agree with the authors' to 7e-11; left as NaN they do not (R5: partial r 0.239 instead of 0.270).
- `ContactStats` bins are half-open [lo, hi), the Nagano et al. bins right-closed (lo, hi]: one contact of the 395 cells
  sits on an edge (R6; no effect on the answer).

## Open issues

- R1 / R2: the authors' default anneal takes 25-33 min per run (one core per model); the prompts leave the schedule length to the
  agent, and 5x / 20x shorter anneals still clear the floors.  Running the authors' code at all needs two fixes on
  Python 3.12 (pinned `numpy<2 cython<3`; gzipped NCC input).
- R3: at fb43ac4, U-Chrom's `schicluster_impute` did not follow the authors' code (r 0.687 / 0.702 against 0.6947).
  Fixed in U-Chrom (PR #130, in the snapshot b6d8873): it returns the authors' `Q` (not symmetrised; r 0.6702 over i < j),
  and symmetrising it as the prompt says gives the reference 0.694665.
- R4: the cell types are weakly separable by 1-Mb Hi-C (purity 0.24 at k = 30, chance 0.07), so
  the floor's window above the no-imputation baseline is narrow (0.230 vs 0.213).  U-Chrom's `embed_contacts` took ~14
  min at fb43ac4; with the fix 325-443 s on 16 threads (jobs 47754119, 47799338), the authors' code 289-304 s on 8.
- R5: the question was reworded from the plan's plain correlation with E1 to the partial correlation given CpG (see
  above), and graded by parity rather than a floor.
- R7: dropped (no run of the authors' SnapHiC on Python 3.12 without patching its pandas code; > 10 min anyway).
- `schicluster==1.3.2` cannot be pip-installed with its dependencies without R (rpy2); agents need `--no-deps` (or
  `RPY2_CFFI_MODE=ABI`, not tried).

## Sherlock jobs (partition xiaojie,normal; logs in `$SCRATCH/uchrom-nd/agentbench/logs/`)

| Job | Script | What |
|---|---|---|
| 47730908 | `schic_setup.sbatch` | venv-schicluster (full install fails on rpy2; --no-deps works), first venv-nucdyn (fails to compile), data_v3 folders, inspection |
| 47731948 | `schic_setup.sbatch` (`ENV_ONLY=1`) | venv-nucdyn with `numpy<2 cython<3` |
| 47731459, 47732465 | `schic_truth.sbatch` | U-Chrom R1; R2-R6 |
| 47733961 | `schic_truth.sbatch` (`ONLY=R5 BUILD=types`) | U-Chrom R5 with the partial correlation; `schicar_cell_types` |
| 47738050 | `schic_truth.sbatch` (`ONLY=R4`) | U-Chrom R4 with k = 30 |
| 47732468, 47733962 | `schic_ref.sbatch` | authors' scHiCluster R3-R5; cooltools / numpy R5, R6 (rerun with the partial correlation) |
| 47734274, 47735651, 47738057 | `schic_ref.sbatch` | R4: no-imputation baseline; the two variants; the authors' path with k = 30 |
| 47731949 | `schic_nucdyn.sbatch` | failed: the `'rU'` file mode (fixed by gzipping the NCC file) |
| 47732961, 47736279 | `schic_nucdyn.sbatch` | authors' NucDynamics, default anneal: Cell 1 seeds 0, 1 (the job cancelled once Cell 1 was written); Cell 2 seed 0 (`CELLSEL=R2 TAG=_def`) |
| 47735245, 47735247 | `schic_nucdyn.sbatch` | authors' NucDynamics, `-temps 100` and `-temps 50 -dyns 50`, both cells, seed 0 |
| 47732313, 47734204, 47734957 | `schic_snaphic.sbatch` | R7 input + U-Chrom R7; authors' SnapHiC with the newest packages (fails), with `pandas<3` (stale output folder; fails in post-processing) |
| 47732306, 47732308 | | cancelled before starting (sbatch `--export` split the comma lists; resubmitted with ":") |

