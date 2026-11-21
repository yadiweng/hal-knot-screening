# Data preparation reproduces original input conventions unless explicitly noted.
# All supervised screening is done INSIDE each outer CV training fold.
load_dataset <- function(name) {
  root <- file.path(project_root(), "data", "raw")
  simple <- c(cpu = "cpu.csv", laheart = "laheart.csv",
              pima = "pima.csv", oecdpanel = "oecdpanel.csv",
              real_estate = "real_estate.csv")
  if (name %in% names(simple)) {
    d <- read.csv(file.path(root, simple[[name]]), na.strings = "?")
    d <- na.omit(d)
    y <- as.numeric(d[[1L]])
    x <- as.matrix(d[, -1L, drop = FALSE])
    storage.mode(x) <- "double"
    return(list(X = x, Y = y, selection = FALSE,
                notes = "Response in first column, as in the legacy analysis"))
  }

  if (name %in% c("automobile", "communities_crime")) {
    filepath <- if (name == "automobile")
      file.path(root, "automobile", "imports-85.data") else
      file.path(root, "communities_crime", "communities.data")
    d <- read.csv(filepath, header = FALSE, na.strings = "?")
    # Communities: remove NONPREDICTIVE identifiers BEFORE complete-case filtering.
    # Otherwise county/community NA values incorrectly exclude ~196 usable rows.
    # With the supplied data, this changes complete cases from 123 to 319.
    if (name == "communities_crime") d <- d[, -(1:5), drop = FALSE]
    d <- na.omit(d)
    keep <- vapply(d, function(z) length(unique(z)) > 1L, logical(1))
    d <- d[, keep, drop = FALSE]
    encoded <- model.matrix(~ . - 1, data = d)
    target_index <- if (name == "automobile") 1L else ncol(encoded)
    y <- as.numeric(encoded[, target_index])
    x <- encoded[, -target_index, drop = FALSE]
    # Legacy: glmnet variable selection before CV; corrected to within each fold.
    return(list(X = x, Y = y, selection = TRUE,
                notes = "Legacy target and encoding; glmnet selection moved inside CV"))
  }

  if (name == "bike_sharing") {
    d <- read.csv(file.path(root, "bike_sharing", "day.csv"))
    d <- na.omit(d)
    # Original removed instant/dteday, but included casual and registered,
    # whose sum is exactly the response cnt: direct target leakage.
    d <- d[, setdiff(names(d), c("instant", "dteday", "casual", "registered")), drop = FALSE]
    if (!"cnt" %in% names(d)) stop("Bike sharing dataset missing response cnt")
    y <- as.numeric(d$cnt)
    x <- model.matrix(~ . - 1, data = d[, setdiff(names(d), "cnt"), drop = FALSE])
    return(list(X = x, Y = y, selection = FALSE,
                notes = "Removed leakage predictors casual and registered"))
  }

  if (name == "liver_disorders") {
    d <- read.csv(file.path(root, "liver_disorders", "bupa.data"),
                  header = FALSE, na.strings = "?", fileEncoding = "ISO-8859-1")
    d <- na.omit(d)
    d <- d[, -7L, drop = FALSE] # Historical script excludes dataset's selector column.
    encoded <- model.matrix(~ . - 1, data = d)
    y <- as.numeric(encoded[, ncol(encoded)])
    x <- encoded[, -ncol(encoded), drop = FALSE]
    # Historical script used an undefined X; this fixes it to X_selected.
    return(list(X = x, Y = y, selection = FALSE,
                notes = "Fixed legacy undefined X; response is sixth raw column"))
  }
  stop("Unknown dataset: ", name)
}

# The legacy analyses used cv.glmnet for automobile and communities_crime.
# To avoid target leakage, fit it on the outer TRAINING fold only.
select_predictors <- function(X_train, y_train, use_glmnet, seed) {
  if (!use_glmnet) return(seq_len(ncol(X_train)))
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install glmnet")
  set.seed(seed)
  fit <- glmnet::cv.glmnet(X_train, y_train, family = "gaussian",
                          nfolds = min(10L, nrow(X_train)))
  keep <- which(as.numeric(coef(fit, s = "lambda.min"))[-1L] != 0)
  if (length(keep) == 0L) {
    # FIXED (2026-11): falling back to ALL predictors (~120 for communities)
    # makes the 3-way HAL basis explode (memory/time-outs). Instead keep the
    # first variable(s) to enter the lasso path.
    beta <- fit$glmnet.fit$beta
    # glmnet stores decreasing lambda: the first column with a nonzero
    # coefficient corresponds to the earliest entrant(s) on the path.
    entrant_cols <- which(colSums(beta != 0) > 0)
    if (length(entrant_cols)) {
      keep <- which(beta[, entrant_cols[1L]] != 0)
    } else {
      # Rare degenerate path: guard against missing/constant columns and zero
      # covariance. This fallback is deterministic, but must be reported.
      corrs <- suppressWarnings(as.numeric(cor(X_train, y_train)))
      corrs[!is.finite(corrs)] <- 0
      keep <- which.max(abs(corrs))
    }
    if (!length(keep) || anyNA(keep)) {
      stop("glmnet fallback could not identify even one valid predictor")
    }
    warning("glmnet selected no variables at lambda.min; retained ", length(keep),
            " variable(s) using the path-entry fallback")
  }
  keep
}
