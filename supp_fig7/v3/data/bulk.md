# Version 3, family `bulk`: data

Questions H1-H4, XM1, P1, P2 (`../questions/bulk.yaml`).  Three new folders under
`$SCRATCH/uchrom-nd/agentbench/data_v3/`, built by `bulk_build.py` (this folder; cooler / pandas / pyarrow / numpy only,
never imports `uchrom`; run in `venv-uchrom` for its pyarrow) in the job `../sherlock/bulk_setup.sbatch` (Sherlock jobs
47731273 + 47734059, 2026-10-13) from registered sources; two shared stores reused.  A scratch folder
`data_v3/bulk_work/` (call tables of H1) and `data_v3/bulk_build_report.json` (counts, md5s, balancing) are written
next to them; neither is given to the agents.

| data key | folder / file the prompt names | questions | built from (registered source) |
| --- | --- | --- | --- |
| `rao_chr21_10kb` | `data_v3/rao2014_imr90_chr21_10kb/IMR90_chr21_10kb_hg19.cool` and `GSE63525_IMR90_HiCCUPS_looplist.txt.gz` | H1, H2, H3 | `rao2014_imr90_chr21` (the 5-kb chr21 slice of `GSE63525_IMR90_combined.hic`, Rao et al. 2014), `rao2014_annotations` |
| `su2020_tsv` | `data_v3/su2020_chr21_tsv/chromosome21.tsv`, `Hi-C_contacts_chromosome21.tsv`, `README_August_2020.txt` | H4, XM1 | `su2020_chr21` (Zenodo 3928890, Su et al. 2020) |
| `takei_rna` | `data_v3/takei2025_fov0_rna/takei2025_fov0_rna_counts.tsv` | P2 (with `takei`) | `takei2025_rna` (Zenodo 7693825 `transcriptomic-data.zip`, Takei et al. 2025) |
| `takei` (shared) | `data/takei2025_fov0/takei2025_fov0.chromdata.zarr` | P2 | `takei2025_fov0` |
| `melanoma` (shared) | `data/hicrna_human_melanoma/hicrna_human_melanoma.chromdata.zarr` | P1 | atlas store `hicrna_human_melanoma` (Guo et al. 2026, GEO GSE311199) |

## `rao2014_imr90_chr21_10kb/`

- `IMR90_chr21_10kb_hg19.cool`: the registered 5-kb slice (`$SCRATCH/uchrom-nd/data/rao2014_imr90/IMR90_chr21_5kb_hg19.cool`,
  raw counts read from the GEO `.hic` over HTTP by `uchrom.datasets._builders.slice_hic`) summed into 10-kb bins with
  `cooler.coarsen_cooler(..., factor=2)` and ICE-balanced with `cooler.balance_cooler(clr, cis_only=True, store=True)`
  (= `cooler balance --cis-only`, its defaults: two diagonals ignored, MAD filter 5, at least 10 non-zero pixels per
  bin); weights in the bin column `weight`.  The same recipe as `benchmarks/hic_structures/mustache_validate.py`
  (`balanced()`), so H1's U-Chrom / mustache-hic numbers can be compared with `mustache_validate.json`.
- `GSE63525_IMR90_HiCCUPS_looplist.txt.gz`: copied, md5-checked (`326c6dbafcb860ad1e726c544f71a5a6`); chr21 has 92
  loops (column `chr1` = `21`, hg19).
- Size: 2.8 MB (4,813 bins, 1,487 masked by the balancing, 1.67 M pixels, 10.27 M contacts; balancing converged).  Build report (counts, masked bins, balancing convergence): `data_v3/bulk_build_report.json`.

## `su2020_chr21_tsv/`

The three files of the registered source `su2020_chr21` as published (md5-checked copies of
`$SCRATCH/uchrom-nd/data/su2020_imr90/`): `chromosome21.tsv` (Z, X, Y in nm of 651 loci of 50 kb for each of 7,591
imaged chromosome copies; undetected loci empty), `Hi-C_contacts_chromosome21.tsv` (651 x 651 Rao 2014 IMR90 read
counts summed into the same loci, symmetric, hg38; 9.5 % of the off-diagonal pairs are 0, every index separation has
a non-zero pair), `README_August_2020.txt`.  ~262 MB.  (Other families may use the `su2020_chr21` ChromData store;
this folder gives the raw tables so that H4 / XM1 name no store.)

## `takei2025_fov0_rna/`

`takei2025_fov0_rna_counts.tsv`: RNA seqFISH transcripts per cell (rows, `cell_id` = `rep_fov_cellID` as in the
`takei` store, all 47 cells) and gene (60 columns), counted from the mRNA spot table of cerebellum replicate 1 in the
Zenodo `transcriptomic-data.zip`: the members whose name contains `cerebellum`, `replicate1` and `mrna` — for the
cerebellum that is one table, `cerebellum-replicate1-mRNA-seqFISH-non-barcoded-ch1.csv` (the 60-gene sequential
smFISH panel; the archive's other cerebellum table is intron seqFISH); a spot of fov `f` and `cellID` `c` belongs to
cell `1_<f>_<c>`.  Every one of the 47 cells has transcripts (total per cell: median 161, range 42-4,649); md5 of the
table `eda7365618087d0fce337a2a1ee95146`.  The zip (462 MB, md5 `c29594d9ed0f4dfea67c0140c56050a1`) is downloaded
inside the job into `$SCRATCH/uchrom-nd/agentbench/dl/takei2025_rna/` (not part of the agents' data).
`uchrom.io.read_seqfish_transcriptome(zip, dataset="cerebellum", assay="mrna", replicate=1)` gives the same totals
(checked in `truth/bulk.py`, P2).

## Registration to add

Nothing new: every input comes from a registered source (`rao2014_imr90_chr21`, `rao2014_annotations`,
`su2020_chr21`, `takei2025_rna`, `takei2025_fov0`, the atlas store `hicrna_human_melanoma`).  For
`apps/atlas/recipes/README.md`, add to the entries of these sources a line "users: agent benchmark v3 (family bulk:
`supp_fig7/v3/data/bulk_build.py`)"; `takei2025_rna` gains its first user outside
`uchrom.io.load_takei2025_cerebellum`.
