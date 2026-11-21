# Execution history

This folder records how the results in the paper were actually produced.

`EXECUTED_SOURCE.sha256` lists SHA-256 checksums of the R, script and Slurm
files as they were on Longleaf when the paper's results were computed. Run
`sha256sum -c history/EXECUTED_SOURCE.sha256` from the repository root to see
which files differ from the executed ones:

- `scripts/make_overleaf_tables.R`: only the table captions were shortened to
  the form used in the paper, and one label changed: "≈ 2" became "2" for the
  skewness of $X_1$ in the independent skewed design. Table contents are
  unchanged.
- `scripts/install_packages.R`: the default fork reference is the commit used for the paper (`a81af69b28639221b4146ba8cfdcfcc83bb68aa2`, as reported by `packageDescription("hal9001")$RemoteSha` on Longleaf) instead of the branch name.
- `slurm/*.sbatch`: the project root is taken from the submission directory
  (or `HAL_ROOT`) instead of a hard-coded Longleaf path.
- `R/simulation_config.R.before_gskew`, `slurm/gskew_only.sbatch`: kept here as
  `R/simulation_config_13groups.R` and `slurm/gskew_only_executed.sbatch`.

## Simulation arrays

1. **Array tasks 1–50 (designs `regular`, `corr_0.1`, `corr_0.5`, `corr_0.9`,
   `skewed`).** These were submitted with
   `slurm/simulations_13groups_executed.sbatch` while
   `R/simulation_config_13groups.R` was the active configuration.
   - That configuration also listed eight legacy normal-mixture groups
     (`skew_*`). These were replaced before they were used in the paper.
   - Their output, if any, is not part of `results/`.
2. **Array tasks 51–170 (the 12 `gskew_*` designs).** These were submitted
   with `slurm/gskew_only_executed.sbatch` and the current
   `R/simulation_config.R`.

The first five designs occupy positions 1–5 in both configurations, with the
same generator and test size. Because seeds are computed from the position of
a design (`seed_for(group_id, repeat, setting)`), their results are identical
to those of one 170-task run with the current `slurm/simulations.sbatch`.

All output rows carry `analysis_version = "2026-11-16-audited-v2"`.

## Real data

The 90 tasks (9 data sets × 10 repeats) were run with
`slurm/real_data_executed.sbatch`. This is identical to `slurm/real_data.sbatch`
apart from the project root.
