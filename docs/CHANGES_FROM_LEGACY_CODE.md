# Changes relative to the 2025 analysis code

The analysis was re-implemented from the authors' 2025 code (`k-means-hal-sim`).
The changes below affect results; each is also marked `FIXED (2026-11)` in the
source where applicable. Results in earlier drafts of the manuscript are not
comparable with the current results.

## Design of the study

| Area | 2025 code | This repository |
|---|---|---|
| Simulation replicates | 50 runs for some designs, single runs for others | 10 independent replicates per design, seeded by `seed_for()` |
| Real-data evaluation | one 10-fold CV | 10 repeats of 10-fold CV with different fold assignments; same folds for all methods and budgets within a repeat |
| Output | summary objects only | one CSV per replicate/repeat, checkpointed, tagged with `analysis_version` |
| Skewed and correlated designs | normal-mixture noise added to correlated covariates (`skewness.R`). This produced \|skewness\| ≤ 0.33, inflated the spread of X, made X3 continuous and removed most of the correlation. | Gaussian copula with standardized Gamma margins (`skew_copula.R`): target skewness 0.5/1/2/3 × latent correlation 0.1/0.5/0.9. X3 stays binary; X1, X2 and X5 keep the mean and variance of the correlated design. |

## Corrections

1. **Univariate jump function** (`correlation.R`, `skewness.R`). A line break
   ended the expression early, so Y had only three of the five jumps and no
   noise. This inflated test $R^2$ to about 0.97 for $p = 1$.
2. **PAM on the real data.** The real-data scripts called k-means twice. They
   now call `algorithm = "pam"`.
3. **Automobile and communities and crime.** Predictor screening by
   `cv.glmnet` used all observations before cross-validation (outcome
   leakage). It is now done inside each training fold. If no predictor is
   selected, the first predictor(s) to enter the lasso path are kept.
4. **Communities and crime.**
   - The five non-predictive columns (state, county, community, community name,
     fold) are now removed *before* removing incomplete observations
     (n = 319 instead of 123).
   - The community-name column had been one-hot encoded into ID dummies; with
     the columns removed, this no longer happens.
5. **Bike sharing.** `casual` and `registered`, whose sum is the outcome `cnt`,
   are excluded.
6. **Liver disorders.** An undefined object `X` in the 2025 script is replaced
   by the intended predictors (columns 1–5). The outcome is `drinks`
   (column 6), and the `selector` column is dropped.
7. **Knot budget.** If a training fold contains only binary predictors,
   `clustering_budget()` uses $d_c = d$ and warns, instead of stopping.

## Not changed (described in the paper)

- **Skewed design with independent variables.** Training and test data come
  from the same distribution. The earlier manuscript text described this
  design as a distribution shift with uniform test data, which the 2025 code
  never implemented. The paper now describes the design as implemented.
- **Clustering knot budgets.** These are derived from the quantile basis
  counts, as in the 2025 code. In the simulations, $d_c = d$.
- **Encoding of categorical predictors.** `model.matrix` is applied to the
  complete data set (unsupervised).
- **Outcome definitions.** The outcomes are those of the 2025 code:
  - automobile: `symboling`
  - oecdpanel: `humancap`
  - other data sets: see `data/README.md`
