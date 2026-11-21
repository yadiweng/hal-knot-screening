#!/bin/bash
#SBATCH --job-name=hal_gskew
#SBATCH --partition=general
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --array=51-170%4
#SBATCH --output=logs/gskew_%A_%a.out
#SBATCH --error=logs/gskew_%A_%a.err
set -eo pipefail   # 'set -u' only after module load (Lmod scripts may use unset variables)
ROOT=/work/users/y/a/yadiw/hal_analysis
cd "$ROOT"
module purge
module load "${R_MODULE:-r}"
set -u
export R_LIBS_USER="$ROOT/.R/library"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# 17 groups x 10 independent Monte Carlo runs. Array ID: 1..170.
# Keep the order identical to R/simulation_config.R (checked at start).
# Only submit tasks 51..170 for the NEW gskew groups. IDs 1..50 are running separately under audited-v2.
groups=(regular corr_0.1 corr_0.5 corr_0.9 skewed \
        gskew_0.5_corr_0.1 gskew_0.5_corr_0.5 gskew_0.5_corr_0.9 \
        gskew_1_corr_0.1 gskew_1_corr_0.5 gskew_1_corr_0.9 \
        gskew_2_corr_0.1 gskew_2_corr_0.5 gskew_2_corr_0.9 \
        gskew_3_corr_0.1 gskew_3_corr_0.5 gskew_3_corr_0.9)
# NREPS repeats per group; to go beyond 10 see README ("More repeats").
NREPS=${NREPS:-10}
if ! [[ "$NREPS" =~ ^[1-9][0-9]*$ ]]; then
  echo "Invalid NREPS=$NREPS; expected positive integer" >&2
  exit 2
fi
# This concurrent add-on keeps the existing v2 array IDs 1..50 unchanged.
# Restrict to the already-running 10-repeat layout to prevent mis-mapped IDs.
if (( NREPS != 10 )); then
  echo "Concurrent add-on expects NREPS=10 (v2 tasks 1..50 already running); got $NREPS" >&2
  exit 2
fi
n_tasks=$(( ${#groups[@]} * NREPS ))
if (( SLURM_ARRAY_TASK_ID < (5 * NREPS + 1) || SLURM_ARRAY_TASK_ID > n_tasks )); then
  echo "gskew task $SLURM_ARRAY_TASK_ID outside valid range $((5*NREPS+1))..$n_tasks" >&2
  exit 2
fi
i=$(( (SLURM_ARRAY_TASK_ID - 1) / NREPS ))
repeat=$(( (SLURM_ARRAY_TASK_ID - 1) % NREPS + 1 ))
Rscript -e 'source("R/simulation_config.R"); g <- commandArgs(TRUE); if (!identical(g, simulation_groups$group)) stop("slurm group list != R/simulation_config.R")' "${groups[@]}"
Rscript scripts/run_simulation.R "${groups[$i]}" "$repeat"
