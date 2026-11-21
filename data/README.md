# Input data

The files in `raw/` are the inputs used for the paper. They are byte-identical
to the files in the authors' 2025 code archive. `checksums.sha256` lists their
SHA-256 checksums; run `sha256sum -c data/checksums.sha256` from the repository
root to verify them.

All preprocessing is in `R/data_prep.R`:

- observations with any missing value are removed;
- constant columns are dropped;
- predictors stored as text become indicator variables (`model.matrix`);
- predictors stored as integer codes are used as numeric.

For automobile and communities and crime, predictors are additionally screened
by a cross-validated lasso inside each training fold
(`select_predictors()`). $n$ and $d$ below are after preprocessing; the number
of continuous predictors (more than two distinct values) is in parentheses.

| File | Data set | Source | $n$ | $d$ | Outcome modeled | Data-set specific steps |
|---|---|---|---|---|---|---|
| `cpu.csv` | Computer Hardware | UCI #29, doi:10.24432/C5830D | 209 | 6 (6) | `prp` (published relative performance), first column | — |
| `automobile/imports-85.data` | Automobile | UCI #10, doi:10.24432/C5B01C | 159 | 55 (15) | `symboling` (insurance risk rating, −3 to +3), first column; treated as continuous | constant `engine-location` dropped; 9 categorical attributes → 40 indicators |
| `communities_crime/communities.data` | Communities and Crime | UCI #183, doi:10.24432/C53W3X | 319 | 122 (122) | `ViolentCrimesPerPop` (normalized to 0–1 in the source), last column | 5 non-predictive columns removed before complete-case filtering; 22 LEMAS police variables are missing for 1,675 of 1,994 communities |
| `liver_disorders/bupa.data` | Liver Disorders | UCI #60, doi:10.24432/C54G67 | 345 | 5 (5) | `drinks` (half-pint equivalents per day), column 6 | `selector` (column 7, a train/test indicator) dropped |
| `bike_sharing/day.csv` | Bike Sharing (daily) | UCI #275, doi:10.24432/C5W894 | 731 | 11 (8) | `cnt` (daily rentals) | `instant`, `dteday`, `casual`, `registered` removed (`cnt = casual + registered`) |
| `real_estate.csv` | Real Estate Valuation | UCI #477, doi:10.24432/C5J30W | 414 | 6 (6) | house price per unit area, first column | — |
| `laheart.csv` | Los Angeles heart data | not UCI; data set `laheart` of Benkeser & van der Laan (2016, Table 1), who cite Afifi & Azen (1979) | 200 | 11 (11) | `death` (0/1), first column | the file ends with one empty line, which is dropped as incomplete (Benkeser & van der Laan list n = 201) |
| `oecdpanel.csv` | OECD cross-country growth panel | not UCI; data set `oecdpanel` of the R package `np` (Hayfield & Racine, 2008), from Liu & Stengos (1999); also used by Benkeser & van der Laan (2016) | 616 | 6 (5) | `humancap` (average secondary school enrolment rate), first column | — |
| `pima.csv` | Pima Indians diabetes (complete cases) | not UCI; data set `pima` of Benkeser & van der Laan (2016, Table 1), who cite Smith et al. (1988) | 392 | 7 (7) | `diabetes` (0/1), first column | this version has no diabetes pedigree function |

The binary outcomes (`laheart`, `pima`) are fitted with the squared-error loss.

`laheart`, `oecdpanel` and `pima` are not part of the UCI Machine Learning
Repository. Together with `cpu`, they are four of the five real data sets
analysed in the original HAL paper (Benkeser & van der Laan, 2016, Table 1).

References:
- Afifi, A. A. and Azen, S. P. (1979). *Statistical Analysis: A Computer Oriented Approach.*
- Benkeser, D. and van der Laan, M. J. (2016). The highly adaptive lasso estimator.
  *2016 IEEE International Conference on Data Science and Advanced Analytics*, 689–696.
  doi:10.1109/DSAA.2016.93
- Hayfield, T. and Racine, J. S. (2008). Nonparametric econometrics: the np package.
  *Journal of Statistical Software* 27(5), 1–32. doi:10.18637/jss.v027.i05
- Liu, Z. and Stengos, T. (1999). Non-linearities in cross-country growth regressions:
  a semiparametric approach. *Journal of Applied Econometrics* 14, 527–538.
- Smith, J. W., Everhart, J. E., Dickson, W. C., Knowler, W. C. and Johannes, R. S. (1988).
  Using the ADAP learning algorithm to forecast the onset of diabetes mellitus.
  *Proceedings of the Annual Symposium on Computer Application in Medical Care*, 261–265.

## Redistribution

- The UCI data sets are distributed under CC BY 4.0, as stated on the UCI pages
  checked for Computer Hardware and Real Estate Valuation. Please confirm the
  license on the other UCI pages before publishing.
- `oecdpanel.csv` corresponds to `oecdpanel` in the R package `np`: the same
  616 rows and 7 variables, with `humancap` moved to the first column. The
  first rows were checked to be identical. An equivalent file can be created
  with:

  ```r
  data("oecdpanel", package = "np")
  oecdpanel$oecd <- as.integer(as.character(oecdpanel$oecd))
  oecdpanel <- oecdpanel[, c("humancap", "growth", "oecd", "year", "initgdp", "popgro", "inv")]
  write.csv(oecdpanel, "data/raw/oecdpanel.csv", row.names = FALSE)
  ```

  The bytes of such a file, and therefore its checksum, differ from the
  bundled file, which quotes the `oecd` column.
- `laheart.csv` and `pima.csv` are copies from the authors' earlier code. They
  are the versions analysed by Benkeser & van der Laan (2016). Their
  redistribution terms have not been verified.
- If redistribution is not permitted, remove the file and give download
  instructions instead; the checksums allow users to verify a downloaded copy.
- The data are not covered by the code license of this repository.
