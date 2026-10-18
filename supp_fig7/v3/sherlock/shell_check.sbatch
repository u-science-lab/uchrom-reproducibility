#!/bin/bash
# The agents shell inside the container (login shell, as Codex runs commands): which python, its packages, the job.
#SBATCH --job-name=ab3-shell-check
#SBATCH --partition=xiaojie,normal,owners
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/wzxu/uchrom-nd/agentbench/logs/%x-%j.out
GH=$GROUP_HOME/$USER/uchrom-nd/agentbench
T=$L_SCRATCH/sc; mkdir -p $T/home $T/work $T/tmp
printf "export PATH=/opt/venv/bin:/opt/codex/bin:\$PATH\nexport VIRTUAL_ENV=/opt/venv\n" > $T/home/.bash_profile
for v in venv-plain venv-uchrom-v3; do
  apptainer exec --containall --cleanenv --home $T/home:/home/agent --pwd /work --bind $T/work:/work --bind $T/tmp:/tmp \
    --bind $GH/$v:/opt/venv:ro --env HOME=/home/agent --env SLURM_JOB_ID=$SLURM_JOB_ID \
    --env PATH=/opt/venv/bin:/usr/local/bin:/usr/bin:/bin $GH/py312full.sif bash -lc \
    "echo $v; which python; echo job=\$SLURM_JOB_ID; python -c \"import importlib.util as u; print({m: bool(u.find_spec(m)) for m in (\x27uchrom\x27, \x27zarr\x27, \x27cooler\x27)})\""
done
