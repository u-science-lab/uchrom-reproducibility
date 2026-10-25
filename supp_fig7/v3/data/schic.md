# Data of the schic family (version 3)

All inputs are copies or derivations of registered sources; nothing new is downloaded.  New folders sit under
`$SCRATCH/uchrom/agentbench/data_v3/` and are built by `v3/truth/schic_build_data.py`
(`v3/slurm/schic_setup.sbatch` job 47730908: stevens2017_cell2, mm10_cpg, nagano2017_hap_2i;
`schic_snaphic.sbatch` job 47732313: nagano2017_dip_chr19; `schic_truth.sbatch` with `BUILD=types` job 47733961:
schicar_cell_types) from the shared data directory (`$SCRATCH/uchrom/data`, the files
of `python -m uchrom.datasets fetch ...`) and the shared data root (`$SCRATCH/uchrom/agentbench/data`).

| Data key | Folder (+ file the prompt names) | Content | From |
|---|---|---|---|
| `stevens` (shared) | `data/stevens2017_cell1/` | `GSM2219497_Cell_1_contact_pairs.txt.gz`, `GSM2219497_Cell_1_genome_structure_model.pdb.gz` | fetch `stevens2017` (unchanged) |
| `stevens_cell2` (new) | `data_v3/stevens2017_cell2/` | `GSM2219498_Cell_2_contact_pairs.txt.gz` (79,569 contacts), `GSM2219498_Cell_2_genome_structure_model.pdb.gz` (10 models, 100-kb particles) | fetch `stevens2017` (`stevens2017_mesc/`), Cell 2 copied unchanged |
| `schicar` (shared) | `data/schicar_mop/` (`contacts_cells_1Mb.scool`, `contacts_all.mcool`, `schicar_mop.chromdata.zarr`) | 5,313 scHiCAR cells, per-cell 1-Mb maps, pseudo-bulk, cell table (`cell_type`) | `build_schicar_mop.py` |
| `mm10_cpg_1mb` (new) | `data_v3/mm10_cpg/mm10_cpg_1Mb.tsv` | `chrom, start, end, cpg_ratio` for every 1-Mb bin of the scHiCAR maps (the `.scool`'s bins; chromosome ends shorter): scHiCluster's CpG ratio (`bedtools nuc -pattern CG -C`): CG occurrences (any case) inside the bin over the bin's length minus its N bases, 0 where the bin has no base | **derived** from fetch `mm10_genome` (UCSC `mm10.fa.gz`) |
| `schicar_types` (new) | `data_v3/schicar_cell_types/schicar_cell_types.tsv` | `cell, cell_type` of the 5,313 scHiCAR cells (cell = the cell name in `contacts_cells_1Mb.scool`), the store's cell table as a plain table | the shared store `schicar_mop.chromdata.zarr` |
| `nagano_dip_chr19` (new, R7, dropped) | `data_v3/nagano2017_dip_chr19/` (`contacts/<cell>.contacts.txt`, `cells.txt`) | 50 of the 742 diploid serum mES cells SnapHiC analysed (Nagano et al. 2017), drawn with `numpy.random.default_rng(0).permutation` of the sorted ids with `cells.snaphic_top742`: per cell SnapHiC's input format (`cell chr19 pos chr19 pos`), one line per non-empty 10-kb pixel of chr19 of the loader's per-cell map (fend pairs lifted mm9 -> mm10, binned at 10 kb), positions at the bin centres (2,077-12,338 pixels per cell, median 5,881) | loader `nagano2017_dip_serum` (fetch `nagano2017_dip_serum`: `nagano2017_dip/`), its `nagano2017_dip_serum_10kb.scool` |
| `nagano_hap_2i` (new) | `data_v3/nagano2017_hap_2i/` (`adj/<cell>.adj`, `GATC.fends`, `cells.tsv`) | the 395 haploid mES cells of Nagano et al. 2017 in condition `2i_all` with `passed_qc == 1` (GEO feature table): the authors' fend-pair tables unchanged (`fend1 fend2 count`), the mm9 GATC fragment ends (`fend chr coord`), `cells.tsv` (cell, condition, rows of the adj table); the published cell-cycle features are **not** included | fetch `nagano2017_hap` (`nagano2017_hap/adj_2i/<cell>/adj`, `GATC.fends`, `GSE94489_haploids_features_table.txt.gz` for the selection) |

The build log with md5 sums: `data_v3/schic_build_log.json` (`mm10_cpg_1Mb.tsv` 83826705a1a05d1b4ae6720b01bbcb32,
2,737 bins of which 65 without sequence; `schicar_cell_types.tsv` 0c298a2045e6dbbcb536bde144df810d;
`nagano2017_hap_2i/cells.tsv` cd314c7bcc0b5c48610270498efb265c, `GATC.fends` 604e877b4c87c9f00ffc15fbff9aca46; the adj
row counts match GEO's `total_contacts` to within 2; Stevens Cell 2 files keep their registered md5s).

## Registration to add

No new external source: every folder above derives from entries already in
`packages/uchrom/uchrom/datasets/_sources.py` (`stevens2017`, `mm10_genome`, `nagano2017_hap`; scHiCAR via the shared data root).  For
`apps/atlas/recipes/README.md`, section `agentbench/data/`, add a version-3 table:

```
### `agentbench/data_v3/` (schic family): inputs of the agent benchmark, version 3

| Folder | Content | From |
|---|---|---|
| `stevens2017_cell2/` | `GSM2219498_Cell_2_contact_pairs.txt.gz`, `GSM2219498_Cell_2_genome_structure_model.pdb.gz` (5.5 MB) | fetch `stevens2017` (`stevens2017_mesc/` above), Cell 2 only |
| `mm10_cpg/mm10_cpg_1Mb.tsv` | CpG ratio (scHiCluster: CG count / non-N bases, 0 without bases) of the 1-Mb bins of `schicar_mop/contacts_cells_1Mb.scool` | **derived** from fetch `mm10_genome` by `supp_fig7/v3/truth/schic_build_data.py` |
| `schicar_cell_types/schicar_cell_types.tsv` | cell id and `cell_type` of the 5,313 scHiCAR cells (the shared store's cell table) | `schicar_mop/` (shared) by the same script `--only types` |
| `nagano2017_hap_2i/` | 395 Nagano 2017 haploid cells (2i_all, passed_qc): `adj/<cell>.adj`, `GATC.fends` (mm9), `cells.tsv` (668 MB) | fetch `nagano2017_hap` (`nagano2017_hap/` above), subset by the same script |
| `nagano2017_dip_chr19/` | (question R7, dropped: only if R7 is revived) 50 Nagano 2017 diploid serum cells (seed-0 draw from SnapHiC's 742): chr19 10-kb pixels in SnapHiC's input format | loader `nagano2017_dip_serum` (`nagano2017_dip/` above), the same script `--only snaphic` |
```

Code the questions point to (not data; pinned in the questions and the report): scHiCluster (PyPI `schicluster==1.3.2`;
GitHub `zhoujt1994/scHiCluster` tag v1.3.4 = a248677867fa4b250bc5b2c5f858871f001e3bb1, MIT), NucDynamics
(`tjs23/nuc_dynamics` @ 00eaf0b57385365c645eee0fda9dd5a2055550a0, LGPL-3.0), schic2 (`tanaylab/schic2` @ ab593d2, the
definitions of R6 only).
