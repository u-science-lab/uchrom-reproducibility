#!/bin/bash
# One-off: the references the spatial statistics are checked against —
#   * SPARK (xzhoulab/SPARK @a8b4bf2, GPL-3; SPARK-X) and its CRAN dependencies into $BASE/ref/Rlib (R/4.4.2);
#   * esda + libpysal (PySAL; Moran's I, Geary's C) into venv-sc, binary wheels only, numpy / pandas / scipy pinned
#     at the versions already there (pip installs and compiling belong in a job, not on a login node).
#   sbatch extended_data_fig2/slurm/setup.sbatch
#SBATCH --job-name=sv-setup
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=%x-%j.out
source "${REPRO:-$SCRATCH/uchrom/repro}/extended_data_fig2/slurm/_env.sh"
mkdir -p "$R_LIBS_USER"
Rscript -e 'want <- c("Rcpp", "RcppArmadillo", "foreach", "doParallel", "CompQuadForm", "matlab", "pracma");
             need <- setdiff(want, rownames(installed.packages(lib.loc = Sys.getenv("R_LIBS_USER"))));
             if (length(need)) install.packages(need, repos = "https://cloud.r-project.org",
                                                lib = Sys.getenv("R_LIBS_USER"), Ncpus = 4)'
SHA=a8b4bf27b804604dfda53da42992f100b8e4e727
curl -sSL -o "${TMPDIR:-/tmp}/SPARK.tar.gz" "https://github.com/xzhoulab/SPARK/archive/$SHA.tar.gz"
tar -xzf "${TMPDIR:-/tmp}/SPARK.tar.gz" -C "${TMPDIR:-/tmp}"
# the repository ships object files built elsewhere (SPARK.so linked to a system liblapack): build from source
rm -f "${TMPDIR:-/tmp}/SPARK-$SHA"/src/*.o "${TMPDIR:-/tmp}/SPARK-$SHA"/src/*.so
# its Makevars asks for C++11; RcppArmadillo 15 needs C++14 or later: use R's default standard
sed -i '/^CXX_STD/d' "${TMPDIR:-/tmp}/SPARK-$SHA"/src/Makevars
R CMD INSTALL -l "$R_LIBS_USER" "${TMPDIR:-/tmp}/SPARK-$SHA"
Rscript -e 'library(SPARK); cat("SPARK", as.character(packageVersion("SPARK")), "\n")'
$PY -m pip freeze | grep -iE "^(numpy|pandas|scipy|scikit-learn)==" > "${TMPDIR:-/tmp}/pins.txt"
cat "${TMPDIR:-/tmp}/pins.txt"
$PY -m pip install --only-binary=:all: -c "${TMPDIR:-/tmp}/pins.txt" "esda==2.10.0" "libpysal==4.15.0"
$PY -c "import esda, libpysal; print('esda', esda.__version__, 'libpysal', libpysal.__version__)"
du -sh "$BASE/ref"
