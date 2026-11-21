#!/bin/bash
#SBATCH --job-name=hal_real
#SBATCH --partition=general
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=1-90%6
#SBATCH --output=logs/real_%A_%a.out
#SBATCH --error=logs/real_%A_%a.err
set -eo pipefail   # 'set -u' only after module load (Lmod scripts may use unset variables)
# Repository root: the directory you submit from (or set HAL_ROOT).
ROOT="${HAL_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$ROOT"
if [ ! -f R/common.R ]; then echo "Submit from the repository root or set HAL_ROOT" >&2; exit 2; fi
module purge
module load "${R_MODULE:-r}"
set -u
export R_LIBS_USER="$ROOT/.R/library"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# 9 datasets x 10 independently generated CV partitions.
datasets=(cpu automobile communities_crime pima oecdpanel bike_sharing \
          liver_disorders laheart real_estate)
# NREPS repeats per group (10 in the paper); see README ("More repeats").
NREPS=${NREPS:-10}
if ! [[ "$NREPS" =~ ^[1-9][0-9]*$ ]]; then
  echo "Invalid NREPS=$NREPS; expected positive integer" >&2
  exit 2
fi
n_tasks=$(( ${#datasets[@]} * NREPS ))
if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > n_tasks )); then
  echo "Array task $SLURM_ARRAY_TASK_ID outside valid range 1..$n_tasks" >&2
  exit 2
fi
i=$(( (SLURM_ARRAY_TASK_ID - 1) / NREPS ))
repeat=$(( (SLURM_ARRAY_TASK_ID - 1) % NREPS + 1 ))
Rscript scripts/run_real_data.R "${datasets[$i]}" "$repeat"
