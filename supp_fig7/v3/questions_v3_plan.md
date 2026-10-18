<!-- Recovered 2026-10-18 from the session record (the working copy was lost): the plan as approved on 2026-10-13,
     with its one later edit.  The questions actually run are those of questions_v3.yaml (Part A extended by F1-F4). -->

# Agent benchmark, version 3: question plan (draft, 2026-10-13)

Lessons of version 2 (`results/agent_bench_v2.json`): 18 of 43 questions were passed by every model in every
condition; most other data questions discriminated through choices their wording left open (mean or median, the
binning, which radial measure, how to fit), not through U-Chrom; B7 and B8 named U-Chrom's own callers and could not
be answered without it; M1 and M2 were solved everywhere with cooltools.  Version 3 therefore has two parts:

- **Part A, reading the open store** (8 questions kept from version 2, fully specified): expected near the ceiling in
  every condition; reported on their own, as evidence that the store is readable without U-Chrom.
- **Part B, published methods** (9 + 2 conditional): each names the method, its paper, its parameters and the input
  preparation, never U-Chrom.  The reference is computed with U-Chrom **and** with the authors' implementation; a
  question enters only if the two agree within its tolerance (otherwise the tolerance covers both, or it is dropped).

Agent: Codex CLI 0.156 (`module load codex`), models gpt-6-astra, gpt-6-sol, gpt-6-luna; conditions none / bare /
docs / skills as in version 2; 3 repeats; 1,800-s limit (method questions run reconstructions).  Judge: a model of
another family (Claude Opus 5.5), numbers tested in code.

## Part A: reading the open store (kept)

| ID | Data | Question (short) | Grade | v2 (min over 12 cells) |
|---|---|---|---|---|
| A1 | Takei 2025 FOV 0 (seqFISH+) | cells, traces, spots | exact | 1.00 |
| A6 | scHiCAR | largest cell type; its median contacts per cell | name; ±1 % | 1.00 |
| B6 | scHiCAR (linked pseudo-bulk map) | fraction of chr1 contacts > 10 Mb apart | ±2 % | 1.00 |
| B4 | NucDynamics cell (store + linked .mcool) | Spearman of contacts vs 3-D distance, chr1, 1 Mb, contacting pairs | ±5 %; sign | 1.00 |
| H1 | HiRES atlas store over HTTP | cells, traces, median spots per cell | exact | 1.00 |
| H2 | HiRES atlas store over HTTP | Rg of each chr7 copy of cell ValaB8w3043 | ±2 % | 1.00 |
| S3 | spatial Hi-C-RNA, melanoma | copy number from contact legs: most variable chromosome, its s.d. | name; ±5 % | 0.83 |
| X1 | scHiCAR + Liu 2025 DNA-MERFISH | short/long proximity, ExN vs Oligo in both; do they agree | ±5 %; yes/no | 1.00 |

## Part B: published methods (new or reworded)

| ID | Method (paper) | Data (registered) | Question (short) | Reference / independent check | Grade |
|---|---|---|---|---|---|
| T1 | ArcFISH TAD caller, arc p-value (Yu et al. 2025) | Takei 2025 FOV 0, chr15 | average each trace into 1-Mb bins; window 4 Mb, FDR 0.1, other defaults; domains at the finest level and their boundaries (replaces B7) | U-Chrom (after the fix) = ArcFISH GitHub a21d889 | exact count; boundaries ±1 bin |
| T2 | ArcFISH TAD caller | Bintu 2018 IMR90 chr21:18.6-20.6 Mb, 30 kb (`bintu_imr90`) | defaults (window 100 kb, FDR 0.1); boundary positions | same | boundaries ±1 bin |
| T3 | ArcFISH loop caller, axis-wise F test | Su 2020 IMR90 chr21, 50 kb (`su2020_chr21`) | defaults (FDR 0.1, 100 kb-1 Mb, gap 50 kb); number of loops, strongest loop's anchors | same | count exact (or ±5 %); anchors ±1 bin |
| T4 | ArcFISH A/B caller, axes PC | Su 2020 IMR90 chr21 | min compartment 1 Mb; A = cluster with the smaller mean pairwise distance (ArcFISH's rule without a reference genome); fraction of loci in A | same | ±2 % |
| T5 | normalised insulation score on median distances (Su et al. 2020) | Bintu 2018 IMR90 | window 150 kb, peaks with prominence 0.1; boundary positions | ArcFISH `method="insulation"` and U-Chrom | ±1 bin |
| R1 | NucDynamics, hierarchical annealing to 1-Mb particles (Stevens et al. 2017) | Stevens 2017 Cell 1 (contacts + published PDB) | reconstruct; Spearman with the published model on chr1 at 1 Mb (= M3) | U-Chrom; floor from published model 1 vs 2 | at least the floor |
| R2 | NucDynamics | Stevens 2017 Cell 2 | same protocol, another cell | same | at least the floor |
| R3 | random walk with restart imputation (scHiCluster, Zhou et al. 2019) | scHiCAR, two named cells, chr1 1 Mb | pad 1, restart 0.5, tolerance 1e-6; Pearson between the two imputed maps (upper triangle) | U-Chrom = `schicluster` (PyPI) | ±2 % |
| H3 | Mustache loop caller (Roayaei Ardakany et al. 2020) | Rao 2014 IMR90 chr21, 10 kb (`rao2014_imr90_chr21`) | defaults (pt 0.1, st 0.88, sz 1.6, oc 2, 2 Mb); number of loops | U-Chrom 157, mustache-hic 1.3.3 165 (`hic_structures/mustache_validate.json`) | range covering both ±5 % |
| R4 (cond.) | scA/B, CpG density of 3-D neighbours (Tan et al. 2018, Dip-C) | HiRES cell ValaB8w3043 (HTTP) | per-particle CpG frequency averaged within a stated radius; a summary of chr7 | U-Chrom vs Dip-C `color` | enters only if both agree |
| S1 (cond.) | SnapHiC (Yu et al. 2021) | Nagano 2017, ~30 cells, chr19 | defaults; number of loops | U-Chrom vs authors' SnapHiC (`sc_structures/snaphic_reference.json`) | enters only if the authors' code finishes within the limit and both agree |

Dropped from version 2 (33 of 43; A1, A6, B6, B4, H1, H2, S3, X1 kept, M3 = R1, B7 reworded as T1): B8 (U-Chrom's
own caller; T4 replaces it on other data); M1, M2 (cooltools, solved everywhere); 10 more passed everywhere (B1, B5,
C2, C3, C5, C10, D4, D5, D8, E8); 20 that discriminated through what their wording left open (B2, B3, C1, C4, C6-C9,
D1-D3, D6, D7, E1-E5, S1, S2).
