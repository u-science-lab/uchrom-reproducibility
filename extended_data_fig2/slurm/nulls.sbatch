#!/bin/bash
# The null stages of spatial_variable.py on the analysis cache ($OUT/cache, read only): cnv_null (both tissues)
# and patch_null (melanoma).  Results go to $W/runs/sv-nulls/results (copied to extended_data_fig2/results
# from the laptop).    REPO=$SCRATCH/uchrom/repo-<name> sbatch --export=ALL nulls.sbatch [n_random]
#SBATCH --job-name=sv-nulls
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=%x-%j.out
source "${REPRO:-$SCRATCH/uchrom/repro}/extended_data_fig2/slurm/_env.sh"
cd "$REPRO"
RES=$W/runs/sv-nulls/results
mkdir -p "$RES"
for tissue in melanoma brain; do
    $PY -u extended_data_fig2/spatial_variable.py --tissue $tissue --out "$OUT/cache" --results "$RES" \
        --stages cnv_null --random "${1:-1000}"
done
$PY -u extended_data_fig2/spatial_variable.py --tissue melanoma --stats "$STATS" --out "$OUT/cache" \
    --results "$RES" --stages patch_null --random "${1:-1000}"
