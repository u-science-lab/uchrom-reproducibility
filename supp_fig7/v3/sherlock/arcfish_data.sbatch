#!/bin/bash
# Version 3 data folders of T2-T5 (v3/data/arcfish.md), from the registered sources (Takei 2021 4DN FOF-CT replicates
# 1 and 2) and the FOF-CT tables the ArcFISH parity check wrote (benchmarks/arcfish_parity/inputs.py of U-Chrom).
#SBATCH --job-name=ab3-arcfish-data
#SBATCH --partition=xiaojie,normal
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
set -euo pipefail
D=$SCRATCH/uchrom-nd/data
IN=$SCRATCH/uchrom-nd/arcfish_parity/inputs
OUT=$SCRATCH/uchrom-nd/agentbench/data_v3
mkdir -p "$OUT/takei2021_fofct" "$OUT/su2020_chr2p" "$OUT/bintu2018"
cp "$D/4DNFIHF3JCBY.csv" "$OUT/takei2021_fofct/4DNFIHF3JCBY_replicate1.csv"
cp "$D/4DNFIQXONUUH.csv" "$OUT/takei2021_fofct/4DNFIQXONUUH_replicate2.csv"
cp "$IN/bintu2018.csv" "$OUT/bintu2018/bintu2018_IMR90_chr21.csv"
# the first 100 traces (in file order) of the chr2 p-arm table, header kept
python3 - "$IN/su2020_chr2p.csv" "$OUT/su2020_chr2p/su2020_chr2p_first100.csv" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
keep, n = [], 0
with open(src) as fh, open(dst, "w") as out:
    for line in fh:
        if line.startswith("#"):
            out.write(line); continue
        tid = line.split(",")[1]
        if tid not in keep:
            if len(keep) == 100:
                continue
            keep.append(tid)
        out.write(line); n += 1
print(len(keep), "traces", n, "spots")
PY
md5sum "$OUT"/takei2021_fofct/* "$OUT"/su2020_chr2p/* "$OUT"/bintu2018/*
