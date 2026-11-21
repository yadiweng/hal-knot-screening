#!/usr/bin/env Rscript
# Recomputes the numbers quoted in the text of the manuscript (Sections 4.4,
# 5.4, 6.10 and 7) and the Discussion table, from the summary tables written by
# scripts/summarize_results.R.
#   Rscript scripts/paper_numbers.R
# Optional inputs, used when present:
#   results/simulations/*.csv  per-replicate results (unstable setting, Sec. 5.4)
#   results/real_data/*.csv    fold-level results (cpu fold-level statements and
#                              the corrected repeated-CV intervals, Sec. 6.10)
# Outputs: results/tables/paper_numbers.txt, results/tables/discussion_table.csv
source("R/common.R")
source("R/data_prep.R")
invisible(project_root())
tab <- "results/tables"
ss <- read.csv(file.path(tab, "simulation_summary.csv"), stringsAsFactors = FALSE)
pd <- read.csv(file.path(tab, "simulation_paired_differences.csv"), stringsAsFactors = FALSE)
rs <- read.csv(file.path(tab, "real_data_summary.csv"), stringsAsFactors = FALSE)
rm_ <- read.csv(file.path(tab, "real_data_repeat_means.csv"), stringsAsFactors = FALSE)

out <- character()
say <- function(...) out <<- c(out, sprintf(...))
f3 <- function(x) formatC(x, format = "f", digits = 3)

# ---------------------------------------------------------------- simulations
pd$method <- sub("_minus_quantile", "", pd$metric)
pd$sig <- ifelse(pd$lower_95 > 0, 1L, ifelse(pd$upper_95 < 0, -1L, 0L))
fam <- function(g) ifelse(g %in% c("regular", "corr_0.1", "corr_0.5", "corr_0.9"), "uniform X1,X2",
                    ifelse(g == "skewed", "skewed independent", "skewed + correlated"))
pd$family <- fam(pd$group)
pd$skew <- suppressWarnings(as.numeric(sub("^gskew_([0-9.]+)_corr_.*$", "\\1", pd$group)))
pd$corr <- suppressWarnings(as.numeric(sub("^.*corr_([0-9.]+)$", "\\1", pd$group)))

say("== Section 4.4: total number of basis functions, k-means / quantile ==")
tb <- reshape(ss[ss$metric == "total_bases", c("group", "dimension", "scenario", "method", "mean")],
              idvar = c("group", "dimension", "scenario"), timevar = "method", direction = "wide")
tb$ratio <- tb$mean.kmeans / tb$mean.quantile
for (d in c(1, 3, 5)) say("simulations p=%d: ratio %s to %s", d, f3(min(tb$ratio[tb$dimension == d])),
                         f3(max(tb$ratio[tb$dimension == d])))
rb <- reshape(rs[rs$metric == "total_bases", c("dataset", "knot_start", "method", "mean")],
              idvar = c("dataset", "knot_start"), timevar = "method", direction = "wide")
rb$ratio <- rb$mean.kmeans / rb$mean.quantile
for (ds in sort(unique(rb$dataset))) say("real %-18s ratio %s to %s", ds,
  f3(min(rb$ratio[rb$dataset == ds])), f3(max(rb$ratio[rb$dataset == ds])))

say("")
say("== Section 5.4: paired differences in test R^2 (method - quantile) ==")
for (fm in unique(pd$family)) for (m in c("kmeans", "pam")) {
  z <- pd[pd$family == fm & pd$method == m, ]
  say("%-20s %-6s better %3d, worse %d, of %d settings; average gain %s", fm, m,
      sum(z$sig == 1), sum(z$sig == -1), nrow(z), f3(mean(z$mean)))
}
u <- pd[pd$family == "uniform X1,X2", ]
avg_u <- aggregate(mean ~ group + method, u, mean)
say("uniform X1,X2: design-average gain between %s and %s; largest single gain %s",
    f3(min(avg_u$mean)), f3(max(avg_u$mean)), f3(max(u$mean)))
g <- u[u$sig == 1, ]
say("uniform X1,X2: %d of %d significant gains are jump/sin with p >= 3", sum(g$dimension >= 3 &
    g$scenario %in% c("jump", "sin")), nrow(g))
l <- u[u$sig == -1, ]
say("uniform X1,X2: %d significant losses, largest %s (p/scenario: %s)", nrow(l), f3(max(-l$mean)),
    paste(unique(paste0(l$dimension, "/", l$scenario)), collapse = ", "))
cop <- pd[pd$family == "skewed + correlated", ]
say("copula: k-means losses at p = %s; non-significant settings all p=1 or jump: %s",
    paste(unique(cop$dimension[cop$sig == -1 & cop$method == "kmeans"]), collapse = ","),
    all(cop$dimension[cop$sig == 0] == 1 | cop$scenario[cop$sig == 0] == "jump"))
sk <- aggregate(mean ~ group + method, pd[pd$family != "uniform X1,X2", ], mean)
w <- reshape(sk, idvar = "group", timevar = "method", direction = "wide")
say("PAM average gain >= k-means in %d of %d skewed designs", sum(w$mean.pam >= w$mean.kmeans), nrow(w))
say("")
say("Gain by regression function (copula; averaged over p and corr):")
for (m in c("pam", "kmeans")) {
  z <- cop[cop$method == m, ]
  sm <- tapply(z$mean[z$scenario == "smooth"], z$skew[z$scenario == "smooth"], mean)
  js <- tapply(z$mean[z$scenario != "smooth"], z$skew[z$scenario != "smooth"], mean)
  say("  %-6s smooth: %s | jump+sin: %s (skewness %s)", m, paste(f3(sm), collapse = ", "),
      paste(f3(js), collapse = ", "), paste(names(sm), collapse = ", "))
}
say("")
say("Correlation effect, p in {3,5} (p = 1 does not depend on corr):")
z <- cop[cop$dimension > 1, ]
for (m in c("kmeans", "pam")) for (s in c(2, 3)) {
  v <- tapply(z$mean[z$method == m & z$skew == s], z$corr[z$method == m & z$skew == s], mean)
  say("  %-6s skewness %g: corr %s -> %s", m, s, paste(names(v), collapse = "/"), paste(f3(v), collapse = " / "))
}
zc <- pd[grepl("^corr_", pd$group) & pd$dimension > 1, ]
for (m in c("kmeans", "pam")) {
  v <- tapply(zc$mean[zc$method == m], zc$corr[zc$method == m], mean)
  say("  %-6s correlated designs: corr %s -> %s", m, paste(names(v), collapse = "/"), paste(f3(v), collapse = " / "))
}
say("")
ab <- aggregate(mean ~ method, ss[ss$metric == "active_bases" & ss$dimension == 5, ], mean)
say("Active basis functions, p=5, averaged over designs: %s",
    paste(sprintf("%s %.1f", ab$method, ab$mean), collapse = "; "))

simfiles <- list.files("results/simulations", "^gskew_3_corr_0\\.1_repeat_[0-9]+\\.csv$", full.names = TRUE)
if (length(simfiles)) {
  r <- do.call(rbind, lapply(simfiles, read.csv, stringsAsFactors = FALSE))
  r <- r[r$dimension == 3 & r$scenario == "smooth", ]
  q <- r$r2[r$method == "quantile"]
  say("Unstable setting (skewness 3, corr 0.1, p=3, smooth): quantile R^2 mean %s, median %s, negative: %s",
      f3(mean(q)), f3(median(q)), paste(f3(sort(q[q < 0])), collapse = ", "))
  say("  medians: k-means %s, PAM %s", f3(median(r$r2[r$method == "kmeans"])), f3(median(r$r2[r$method == "pam"])))
} else say("Unstable setting: per-replicate files not found in results/simulations (skipped)")

# ---------------------------------------------------------------- real data
say("")
say("== Section 6.10: real data ==")
m <- reshape(rs[rs$metric == "mse", c("dataset", "knot_start", "method", "mean")],
             idvar = c("dataset", "knot_start"), timevar = "method", direction = "wide")
m$rel_k <- 100 * (m$mean.kmeans / m$mean.quantile - 1)
m$rel_p <- 100 * (m$mean.pam / m$mean.quantile - 1)
m$either <- m$mean.kmeans < m$mean.quantile | m$mean.pam < m$mean.quantile
agg <- aggregate(cbind(rel_k, rel_p, either) ~ dataset, m, mean)
for (i in seq_len(nrow(agg))) say("%-18s avg MSE diff: k-means %+.1f%%, PAM %+.1f%%; K-means or PAM lower at %3.0f%% of budgets",
                                  agg$dataset[i], agg$rel_k[i], agg$rel_p[i], 100 * agg$either[i])
cpu <- m[m$dataset == "cpu", ]
say("cpu quantile MSE: budget 20 = %.0f; other budgets %.0f to %.0f", cpu$mean.quantile[cpu$knot_start == 20],
    min(cpu$mean.quantile[cpu$knot_start > 20]), max(cpu$mean.quantile[cpu$knot_start > 20]))
rr <- reshape(rm_[rm_$dataset == "cpu", c("repeat_id", "knot_start", "method", "mse")],
              idvar = c("repeat_id", "knot_start"), timevar = "method", direction = "wide")
say("cpu: both clustering methods lower than quantile in %d of %d (repeat, budget) pairs",
    sum(rr$mse.kmeans < rr$mse.quantile & rr$mse.pam < rr$mse.quantile), nrow(rr))

foldfiles <- list.files("results/real_data", "_repeat_[0-9]+\\.csv$", full.names = TRUE)
if (length(foldfiles)) {
  fo <- do.call(rbind, lapply(foldfiles, read.csv, stringsAsFactors = FALSE))
  corrected <- function(ds, kn, meth) {
    z <- fo[fo$dataset == ds & fo$knot_start == kn, ]
    w <- reshape(z[, c("repeat_id", "fold", "method", "mse", "n_train", "n_valid")],
                 idvar = c("repeat_id", "fold"), timevar = "method", direction = "wide")
    x <- w[[paste0("mse.", meth)]] - w$mse.quantile; J <- length(x)
    h <- qt(0.975, J - 1) * sqrt((1 / J + mean(w$n_valid.quantile) / mean(w$n_train.quantile)) * var(x))
    c(share_lower = mean(x < 0), est = 100 * mean(x) / mean(w$mse.quantile),
      lo = 100 * (mean(x) - h) / mean(w$mse.quantile), hi = 100 * (mean(x) + h) / mean(w$mse.quantile))
  }
  for (kn in c(20, 50)) if (any(fo$dataset == "cpu" & fo$knot_start == kn)) {
    v <- corrected("cpu", kn, "kmeans")
    say("cpu budget %d, k-means: lower fold MSE in %.0f%% of folds; corrected CI %.1f%% [%.1f%%, %.1f%%]",
        kn, 100 * v[["share_lower"]], v[["est"]], v[["lo"]], v[["hi"]])
  }
  n_sig <- 0L; n_all <- 0L
  for (ds in unique(fo$dataset)) for (kn in unique(fo$knot_start[fo$dataset == ds])) for (meth in c("kmeans", "pam")) {
    v <- corrected(ds, kn, meth); n_all <- n_all + 1L
    n_sig <- n_sig + (v[["hi"]] < 0 || v[["lo"]] > 0)
  }
  say("corrected repeated-CV test: %d significant of %d comparisons", n_sig, n_all)
} else say("Fold-level files not found in results/real_data (cpu fold statements and corrected tests skipped)")

# ---------------------------------------------------------------- Discussion table
skew3 <- function(x) { n <- length(x); x <- x - mean(x); (mean(x^3) / mean(x^2)^1.5) * ((n - 1) / n)^1.5 }
datasets <- c("cpu", "automobile", "communities_crime", "pima", "oecdpanel",
              "bike_sharing", "liver_disorders", "laheart", "real_estate")
disc <- do.call(rbind, lapply(datasets, function(ds) {
  X <- load_dataset(ds)$X; X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
  cont <- apply(X, 2, function(v) length(unique(v)) > 2)
  sk <- apply(X[, cont, drop = FALSE], 2, skew3)
  cm <- cor(X[, cont, drop = FALSE]); cm <- cm[upper.tri(cm)]
  data.frame(dataset = ds, kmeans_or_pam_better_pct = round(100 * agg$either[agg$dataset == ds]),
             average_skewness = round(mean(sk), 2), pct_abs_skew_gt_0.5 = round(100 * mean(abs(sk) > 0.5)),
             pct_abs_corr_gt_0.5 = round(100 * mean(abs(cm) > 0.5)))
}))
atomic_write_csv(disc, file.path(tab, "discussion_table.csv"))
say("")
say("== Discussion table (also written to results/tables/discussion_table.csv) ==")
op <- options(width = 160)
out <- c(out, capture.output(print(disc, row.names = FALSE)))
options(op)
writeLines(out, file.path(tab, "paper_numbers.txt"))
cat(out, sep = "\n")
