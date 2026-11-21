# Shared helpers: no packages required. Run all scripts from repository root.
# Increment whenever methods or data-generation change; old files must not be reused.
analysis_version <- "2026-11-16-audited-v2"

assert_checkpoint_version <- function(rows, path) {
  if (is.null(rows) || nrow(rows) == 0L) return(invisible(TRUE))
  if (!("analysis_version" %in% names(rows)) ||
      anyNA(rows$analysis_version) ||
      any(rows$analysis_version != analysis_version)) {
    stop("Stale checkpoint: ", path, ". Back it up to a different directory; ",
         "do NOT mix old and updated results.")
  }
  invisible(TRUE)
}

project_root <- function() {
  required <- file.path(getwd(), "R", "common.R")
  if (!file.exists(required)) stop("Run from the repository root: hal_analysis/")
  normalizePath(getwd(), mustWork = TRUE)
}

check_hal <- function() {
  if (!requireNamespace("hal9001", quietly = TRUE)) {
    stop("Missing HAL fork. See README.md: remotes::install_github('tq21/hal9001@kmeans-knot-screen')")
  }
  invisible(TRUE)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- tempfile(pattern = "._write_", tmpdir = dirname(path), fileext = ".csv")
  write.csv(x, temp, row.names = FALSE, na = "")
  if (!file.rename(temp, path)) {
    unlink(temp)
    stop("Could not write: ", path)
  }
}

# Matching the original code: use the *quantile* basis count to set the
# clustering knot budget. The denominator treats selected columns as axes.
# Missing interaction orders are NA because HAL ignores degrees > dimension.
clustering_budget <- function(quant_basis, n_continuous, n_columns) {
  # FIXED (2026-11): if fold-wise screening keeps only binary predictors the
  # old code stopped, and the array task could never finish. Fall back to
  # treating all columns as axes and warn.
  if (n_columns < 1L) stop("No selected predictors")
  if (n_continuous == 0) {
    warning("Only binary predictors in this fold; approximating knot budget ",
            "using all selected axes. Validate clustering behavior in this case.")
    n_continuous <- n_columns
  }
  degrees <- vapply(quant_basis, function(b) length(b$cols), integer(1))
  counts <- tabulate(degrees, nbins = 3L)
  ans <- rep(NA_real_, 3L)
  for (k in seq_len(min(3L, n_columns))) {
    denom <- if (k == 1L) n_continuous else choose(n_columns, k)
    if (is.finite(denom) && denom > 0) ans[k] <- counts[k] / denom
  }
  ans
}

make_balanced_folds <- function(n, v = 10L, seed = 1L) {
  if (n < v) stop("Not enough observations for ", v, " folds")
  set.seed(seed)
  sample(rep(seq_len(v), length.out = n))
}

# Do not overwrite completed analysis in a partially rerun array task.
completed_key <- function(df, keys, values, n_methods = 3L) {
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  ok <- rep(TRUE, nrow(df))
  for (j in seq_along(keys)) ok <- ok & df[[keys[j]]] == values[j]
  sum(ok) == n_methods && length(unique(df$method[ok])) == n_methods
}

load_previous <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

append_checkpoint <- function(current, new_rows, path) {
  x <- if (is.null(current)) new_rows else rbind(current, new_rows)
  atomic_write_csv(x, path)
  x
}

seed_for <- function(group_id, repeat_id, setting_id = 1L, extra = 0L) {
  # Stable across independent Slurm array tasks and checkpoint restarts.
  as.integer(202611L + group_id * 1000000L + repeat_id * 10000L + setting_id * 100L + extra)
}

fit_hal_screened <- function(train_x, train_y, test_x, method, max_degree, knot_vector) {
  screen <- method != "quantile"
  args <- list(x = train_x, max_degree = max_degree,
               smoothness_orders = 0L, num_knots = knot_vector,
               screen_knots = screen)
  if (screen) args$algorithm <- method
  bases <- do.call(hal9001::enumerate_basis, args)
  fit <- hal9001::fit_hal(X = train_x, Y = train_y,
                         family = "gaussian", basis_list = bases)
  prediction <- as.numeric(predict(fit, new_data = test_x))
  list(prediction = prediction,
       total_bases = length(fit$basis_list),
       active_bases = sum(fit$coefs[-1] != 0),
       basis_list = bases)
}
