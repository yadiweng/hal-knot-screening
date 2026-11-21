#!/bin/bash
#SBATCH --job-name=hal_sim
#SBATCH --partition=general
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --array=1-130%10
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
set -eo pipefail   # 'set -u' only after module load (Lmod scripts may use unset variables)
ROOT=/work/users/y/a/yadiw/hal_analysis
cd "$ROOT"
module purge
module load "${R_MODULE:-r}"
set -u
export R_LIBS_USER="$ROOT/.R/library"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# 13 groups x 10 independent Monte Carlo runs. Array ID: 1..130.
groups=(regular corr_0.1 corr_0.5 corr_0.9 skewed skew_0.5_corr_0.1 \
        skew_0.5_corr_0.5 skew_0.5_corr_0.9 skew_2_corr_0.1 \
        skew_2_corr_0.5 skew_2_corr_0.9 skew_3_corr_0.1 skew_5_corr_0.1)
# NREPS repeats per group; to go beyond 10 see README ("More repeats").
NREPS=${NREPS:-10}
if ! [[ "$NREPS" =~ ^[1-9][0-9]*$ ]]; then
  echo "Invalid NREPS=$NREPS; expected positive integer" >&2
  exit 2
fi
n_tasks=$(( ${#groups[@]} * NREPS ))
if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > n_tasks )); then
  echo "Array task $SLURM_ARRAY_TASK_ID outside valid range 1..$n_tasks" >&2
  exit 2
fi
i=$(( (SLURM_ARRAY_TASK_ID - 1) / NREPS ))
repeat=$(( (SLURM_ARRAY_TASK_ID - 1) % NREPS + 1 ))
Rscript scripts/run_simulation.R "${groups[$i]}" "$repeat"
