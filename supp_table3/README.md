# Supplementary Table 3: the analysis toolkit against reference implementations

Each row of the table (written in the manuscript) comes from one of these checks, each writing its own `results/`:

* `hic_structures/`: compartments, insulation, HiCCUPS, Mustache on Hi-C maps against cooltools and the authors' code;
* `sc_structures/`: scA/B, imputed-map boundaries, SnapHiC against the authors' code;
* `arcfish_parity/`: ArcFISH loops and domains against the authors' package;
* `schicluster_parity/`: the scHiCluster embedding against the authors' package;
* `nucdyn/`: the NucDynamics engine against the original program.

The spatial statistics row comes from `extended_data_fig2/results/validate_stats.json`.  Heavy runs: each check's
`sherlock/`.
