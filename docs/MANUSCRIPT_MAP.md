# Where each table, figure and number of the paper comes from

Table numbers refer to the arXiv version of the paper. All paths are relative
to the repository root.

| Paper | Content | Source data | Produced by | Output |
|---|---|---|---|---|
| Table 1 | regular design | `results/tables/simulation_summary.csv`, `simulation_paired_differences.csv` (group `regular`) | `scripts/make_overleaf_tables.R` | `results/overleaf/tables/sim_regular.tex` |
| Tables 2–4 | correlated designs | same files, groups `corr_0.1`, `corr_0.5`, `corr_0.9` | same | `sim_corr_*.tex` |
| Table 5 | sample correlation matrices of one simulated data set per `corr` | taken over from the earlier manuscript; not regenerated | — | — |
| Table 6 | skewed design, independent variables | group `skewed` | `make_overleaf_tables.R` | `sim_skewed.tex` |
| Tables 7–18 | skewed and correlated designs | groups `gskew_*_corr_*` | `make_overleaf_tables.R` | `sim_gskew_*.tex` |
| Table 19 | summary by design | `simulation_paired_differences.csv` | `make_overleaf_tables.R` | `sim_overview.tex` |
| Table 20 | data sets ($n$, $d$, outcome) | `data/raw/`, `R/data_prep.R` | $n$ and $d$: `load_dataset()` (also listed in `real_overview.tex`) | entered in the manuscript |
| Tables 21–29 | real-data results per data set | `results/tables/real_data_summary.csv`; ⁺/⁻ markers need `results/real_data/*.csv` | `make_overleaf_tables.R` | `real_<dataset>.tex` |
| Table 30 | comparison across data sets (Wilcoxon, sign test, Spearman) | `real_data_summary.csv`, `data/raw/` | `make_overleaf_tables.R` | `real_across.tex` |
| Table 31 | Discussion table: share of knot budgets with lower MSE; skewness and correlation of the continuous predictors | `real_data_summary.csv`, `data/raw/` | `scripts/paper_numbers.R` | `results/tables/discussion_table.csv` |
| Figure 1 | memory of the HAL design matrix | — | plot code not part of this repository | `mem_plot.pdf` in the Overleaf project |
| Appendix A | the nine regression functions | `R/generators/*.R` (identical in all four generators) | — | — |
| Section 4.4 | ratio of basis counts (clustering vs quantile) | summary tables | `paper_numbers.R` | `results/tables/paper_numbers.txt` |
| Section 5 | optimal test $R^2$ (0.80 under the regular design; ranges for other designs) | generators | `scripts/optimal_r2.R` | `results/tables/optimal_r2.csv` |
| Section 5.3 | Spearman/Pearson correlation and skewness of the copula design | generators | `scripts/dgp_distribution_check.R`, `preflight_audit.R` | `dgp_distribution_check.csv`, `dgp_diagnostics.csv` |
| Section 5.4 | counts of significant differences; trends by regression function and correlation; active bases; unstable setting | summary tables; unstable setting needs `results/simulations/gskew_3_corr_0.1_repeat_*.csv` | `paper_numbers.R` | `paper_numbers.txt` |
| Section 6.10 | average MSE differences; cpu details; number of significant corrected tests (0 of 162) | summary tables; cpu fold shares and corrected intervals need `results/real_data/*.csv` | `paper_numbers.R` | `paper_numbers.txt` |
| Section 7 | Spearman correlation of skewness and MSE difference | as Table 30 | `make_overleaf_tables.R` | `real_across.tex` |

`make_overleaf_tables.R` also writes `real_overview.tex` and two figures
(`sim_delta_r2.pdf`, `real_rel_mse.pdf`). These are not used in the current
version of the paper.
