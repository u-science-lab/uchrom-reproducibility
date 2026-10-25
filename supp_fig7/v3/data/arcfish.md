# Version 3 data of T1-T5 (ArcFISH methods)

Folders under `$SCRATCH/uchrom/agentbench/data_v3/` (built by `v3/slurm/arcfish_data.sbatch`) and the data key
of each question (the runner binds the folder read-only at `./data/<folder>`).

| Key | Folder / file the prompt names | Content | From |
|---|---|---|---|
| `takei` | `takei2025_fov0/takei2025_fov0.chromdata.zarr` (the shared data root) | Takei 2025 cerebellum FOV 0 store | `build_takei2025_fov0.py` |
| `takei2021_fofct` | `takei2021_fofct/` (a folder): `4DNFIHF3JCBY_replicate1.csv`, `4DNFIQXONUUH_replicate2.csv` | the two 4DN FOF-CT core tables of Takei et al. 2021 mESC (25 kb, 1,200 loci), unchanged | fetch `takei`, `takei_rep2` (registered in the U-Chrom branch Nanguage/arcfish-parity) |
| `su2020_chr2p` | `su2020_chr2p/su2020_chr2p_first100.csv` | FOF-CT core table (nm) of the first 100 traces (file order) of Su et al. 2020's chromosome 2 p-arm replicate | fetch `su2020_chr2_parm` → `uchrom.io.read_su2020_tracing(unit="nm", extras=False)` → `to_fofct` (`benchmarks/arcfish_parity/inputs.py`) → first 100 traces |
| `bintu2018` | `bintu2018/bintu2018_IMR90_chr21.csv` | FOF-CT core table (nm) of Bintu et al. 2018 IMR90 chr21:18.6-20.6 Mb, 30 kb, 1,277 traces | fetch `bintu_imr90` → `read_bintu_tracing(chrom="chr21", start_bp=18_627_714, resolution_bp=30_000)` → `to_fofct` |

References: `v3/results/arcfish.json` (`v3/truth/arcfish.py` from the parity runs 47732259 and 47738660 of
`benchmarks/arcfish_parity`, both implementations; T5 at prominence 0.01 from the identical per-locus insulation
scores of both, to be recomputed with the final snapshot).
