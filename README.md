# Knot point screening for the highly adaptive lasso

Code and results for the paper

> Y. Weng and S. Qiu. *A comparison of knot point screening algorithms for highly adaptive lasso.* (arXiv link to be added)

The paper compares three outcome-blind knot screening rules for the highly
adaptive lasso (HAL):

- quantile-based discretization (the default in `hal9001`)
- $k$-means clustering
- partitioning around medoids (PAM)

The comparison covers 17 simulation designs (10 independent replicates each)
and 9 real data sets (10 repeats of 10-fold cross-validation).

This repository contains:
- the exact analysis code that produced the paper's results;
- the input data;
- the summary results;
- the scripts that turn those results into every table and quoted number in the paper.

---

## Repository layout

```
R/
  common.R               shared helpers: seeds, HAL fit wrapper, clustering knot budget
  data_prep.R            loading and preprocessing of the 9 real data sets
  simulation_config.R    the 17 simulation designs (order defines the seeds)
  generators/            data-generating mechanisms
    sim_data.R             regular design (Benkeser & van der Laan, 2016)
    correlation.R          correlated design (Gaussian copula, uniform margins)
    skewness_sim_nofactor.R  skewed design with independent variables
    skew_copula.R          skewed and correlated designs (Gaussian copula, standardized Gamma margins)
    skewness.R             legacy normal-mixture design; not used for the paper,
                           only by checks in scripts/preflight_audit.R
scripts/
  install_packages.R     installs the HAL fork and dependencies into ./.R/library
  check_environment.R    verifies that the fork exposes k-means/PAM screening
  preflight_audit.R      fast data/DGP checks; writes results/tables/dgp_diagnostics.csv
  dgp_distribution_check.R  realized skewness/correlation of each design
  run_simulation.R       one simulation replicate of one design (9 settings x 3 methods)
  run_real_data.R        one repeat of 10-fold CV for one data set (9 knot budgets x 3 methods)
  summarize_results.R    combines per-replicate output into results/tables/*.csv
  check_precision.R      optional: were 10 replicates enough?
  make_overleaf_tables.R all result tables of the paper (LaTeX) + optional figures
  paper_numbers.R        every number quoted in the text, and the Discussion table
  optimal_r2.R           optimal ("oracle") test R^2 of every simulation setting
slurm/                   Slurm array scripts (UNC Longleaf)
data/raw/                input data (see data/README.md)
results/
  simulations/           per-replicate simulation output (one CSV per design x replicate)
  real_data/             fold-level real-data output (one CSV per data set x repeat)
  tables/                summary tables used in the paper
history/                 scripts exactly as executed for the paper, and their checksums
docs/MANUSCRIPT_MAP.md   which script and file produce each table/figure of the paper
docs/CHANGES_FROM_LEGACY_CODE.md  corrections relative to the 2025 code
```

## Software

The results were produced with R 4.5.0 on UNC Longleaf (RHEL 9.7, OpenBLAS).
Package versions are listed in `results/tables/sessionInfo.txt`; the main ones are:

| Package | Version | Note |
|---|---|---|
| hal9001 | 0.4.5 | **fork** `tq21/hal9001`, branch `kmeans-knot-screen`, commit `a81af69` (`a81af69b28639221b4146ba8cfdcfcc83bb68aa2`); adds `screen_knots` and `algorithm = "kmeans"/"pam"` to `enumerate_basis` |
| glmnet | 4.1-8 | |
| cluster | 2.1.8.1 | |
| MASS | 7.3-65 | |
| rlang | 1.1.6 | |
| ggplot2 | any recent | needed by `make_overleaf_tables.R` (it also draws two summary figures); install with `install.packages("ggplot2")` |

The CRAN version of `hal9001` does **not** contain the k-means/PAM screening.

Install into a project library:

```bash
export R_LIBS_USER="$PWD/.R/library"
Rscript scripts/install_packages.R      # installs the fork at the commit used for the paper
Rscript scripts/check_environment.R
```

## Quick start: reproduce the tables from the included results (no HAL needed)

```bash
Rscript scripts/make_overleaf_tables.R   # -> results/overleaf/tables/*.tex (Tables 1-4, 6-19, 21-30)
Rscript scripts/paper_numbers.R          # -> results/tables/paper_numbers.txt, discussion_table.csv (Table 31)
Rscript scripts/optimal_r2.R             # -> results/tables/optimal_r2.csv (about 10 min)
Rscript scripts/dgp_distribution_check.R # -> results/tables/dgp_distribution_check.csv
```

The corrected repeated-CV tests (the ⁺/⁻ markers in the real-data tables, and
the cpu statements in Section 6.10) need the fold-level files in
`results/real_data/`. The unstable-setting statement in Section 5.4 needs the
per-replicate files in `results/simulations/`.

`docs/MANUSCRIPT_MAP.md` lists which script and which file produce each table.

## Full reproduction (Slurm)

Run all commands from the repository root. The scripts use the submission
directory as the project root; alternatively, set `HAL_ROOT`.

```bash
module load r                                   # R 4.5.0 on Longleaf
export R_LIBS_USER="$PWD/.R/library"
mkdir -p logs
Rscript scripts/preflight_audit.R               # seconds; data and generator checks
Rscript scripts/run_simulation.R regular 1      # optional smoke test (one replicate)

sbatch slurm/simulations.sbatch   # 17 designs x 10 replicates = 170 array tasks
sbatch slurm/real_data.sbatch     # 9 data sets x 10 CV repeats = 90 array tasks

# after all tasks have finished
Rscript scripts/summarize_results.R 10
Rscript scripts/make_overleaf_tables.R
Rscript scripts/paper_numbers.R
```

- **Output:** each array task writes its own CSV and checkpoints after every
  setting. Resubmitting a task skips completed settings.
- **Versioning:** every output row carries `analysis_version`
  (`R/common.R`). The scripts refuse to mix output from different versions.
- **Resources:** the `#SBATCH` lines (partition `general`, 12–16 GB, 1–2 days)
  are what was used on Longleaf; adjust them to your cluster.
- **Different R module:** `sbatch --export=ALL,R_MODULE=r/4.5.0 slurm/...`.

### More replicates

Seeds depend only on the design (or data set), the replicate number and the
setting. Replicates 11–20 can therefore be added without changing 1–10, e.g.:

```bash
sbatch --array=1-340%10 --export=ALL,NREPS=20 slurm/simulations.sbatch
Rscript scripts/summarize_results.R 20
```

## Design summary

**Simulations**
- `run_simulation.R` uses $n = 1000$ training observations per replicate.
- The test set has 10,000 observations; for the skewed design with independent
  variables it has 1,000.
- Settings: $p \in \{1,3,5\}$ crossed with the smooth, jump and sinusoidal
  regression functions of Benkeser & van der Laan (2016) (listed in Appendix A
  of the paper), with $Y = f_0(X) + \varepsilon$, $\varepsilon \sim N(0,1)$.
- HAL settings: zero-order basis, maximum interaction degree 3.
- Quantile knot budget: (50, 15, 5). The k-means and PAM budgets are derived
  from the quantile basis (`clustering_budget()` in `R/common.R`).
- Seed of each fit: `seed_for(group, replicate, setting)`.

| Groups | Generator | Designs |
|---|---|---|
| `regular` | `sim_data.R` | independent covariates |
| `corr_{0.1,0.5,0.9}` | `correlation.R` | latent correlation 0.1/0.5/0.9, uniform $X_1, X_2$ |
| `skewed` | `skewness_sim_nofactor.R` | independent shifted-Gamma covariates |
| `gskew_{0.5,1,2,3}_corr_{0.1,0.5,0.9}` | `skew_copula.R` | target skewness × latent correlation |

**Real data** (`run_real_data.R`, `R/data_prep.R`)
- 10 repeats of 10-fold CV. Within a repeat, all methods and knot budgets use
  the same folds.
- First-degree knot budget from 20 to 100 in steps of 10, with $(B_2, B_3) = (10, 5)$.
- For automobile and communities and crime, predictors are screened by
  `cv.glmnet` inside each training fold.
- Preprocessing of each data set is described in `data/README.md` and in
  Section 6 of the paper.

## How the paper's results were produced

The paper's results were produced with the code in this repository. Checksums
of the executed files are in `history/EXECUTED_SOURCE.sha256`.

- **Identical to the executed versions:** all files in `R/`, and all scripts
  in `scripts/` except `make_overleaf_tables.R`.
- **Changed for this repository:**
  - `scripts/make_overleaf_tables.R`: table captions only. The table contents
    are unchanged.
  - `slurm/*.sbatch`: the project root is no longer hard-coded.
- **Added:** `scripts/paper_numbers.R` and `scripts/optimal_r2.R`.

The simulations were run in two arrays:

1. Designs 1–5 (`regular`, `corr_*`, `skewed`; array tasks 1–50) were run with
   an earlier 13-group configuration. These are
   `history/slurm/simulations_13groups_executed.sbatch` and
   `history/R/simulation_config_13groups.R`.
2. The 12 `gskew_*` designs (tasks 51–170) were run with
   `history/slurm/gskew_only_executed.sbatch` and the current
   `R/simulation_config.R`.

Designs 1–5 have the same position, generator and test size in both
configurations, so their seeds and results are identical to those of a single
170-task run with `slurm/simulations.sbatch`. The legacy normal-mixture groups
of the 13-group configuration are not part of the paper.

## Results files

| File | Content |
|---|---|
| `results/tables/simulation_summary.csv` | per design × setting × method: mean, SD and 95% t-interval of test $R^2$; total and active basis counts |
| `results/tables/simulation_paired_differences.csv` | paired differences (k-means/PAM − quantile) in test $R^2$ with 95% t-intervals |
| `results/tables/real_data_repeat_means.csv` | per data set × CV repeat × knot budget × method: CV MSE (mean of the 10 fold MSEs) and basis counts |
| `results/tables/real_data_summary.csv` | mean and SD over the 10 CV repeats |
| `results/tables/dgp_diagnostics.csv`, `dgp_distribution_check.csv` | realized properties of the simulated covariates |
| `results/tables/optimal_r2.csv` | optimal test $R^2$ of every simulation setting |
| `results/tables/paper_numbers.txt`, `discussion_table.csv` | numbers quoted in the text; Table 31 |

The SD over CV repeats describes sensitivity to the fold assignment. It is not
a standard error of the population MSE. For comparisons, the paper uses the
corrected repeated-CV t-test (Nadeau & Bengio, 2003; Bouckaert & Frank, 2004),
implemented in `make_overleaf_tables.R` and `paper_numbers.R`.

## Data

See `data/README.md` for sources, the outcome modeled for each data set,
preprocessing, and redistribution terms.

## License

Code: MIT (see `LICENSE`). The data in `data/raw/` are third-party data and are
**not** covered by this license; see `data/README.md`.

## Citation

If you use this code, please cite the paper above, Benkeser & van der Laan
(2016) for HAL, and Hejazi, Coyle & van der Laan (2020) for `hal9001`.
