#!/usr/bin/env Rscript
# Overleaf project with all result tables (same columns as the original xtable
# tables) and the two overview figures.
#   Rscript scripts/make_overleaf_tables.R        # simulations: mean (SD)  [default]
#   Rscript scripts/make_overleaf_tables.R ci     # simulations: mean [95% CI]
# Inputs : results/tables/*.csv (from summarize_results.R) and, for the real-data
#          confidence intervals, the fold-level files results/real_data/*_repeat_*.csv
# Output : results/overleaf/  (main.tex, tables/*.tex, figures/*.pdf)
#          -> zip the folder and use Overleaf "New project > Upload project".
#
# Variability shown (no extra columns)
#   Simulations: test R^2 mean (SD over the 10 independent replicates), or the
#     95% t-CI with argument `ci`. k-means/PAM cells carry ^+ / ^- when the PAIRED
#     95% CI of the difference from Quantile (same data) lies above / below 0.
#   Real data: CV-MSE mean (SD over the 10 CV repeats). ^+ / ^- = significantly
#     lower / higher MSE than Quantile by the corrected repeated-CV t-test
#     (Nadeau & Bengio 2003), which accounts for overlapping training sets;
#     needs the fold-level files, otherwise no markers are printed.
source("R/common.R")
source("R/simulation_config.R")
source("R/data_prep.R")
suppressPackageStartupMessages(library(ggplot2))
project_root()
style <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[[1L]] else "sd"
if (!style %in% c("sd", "ci")) stop("argument must be 'sd' or 'ci'")
tab_dir <- "results/tables"; out <- "results/overleaf"
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
methods <- c("quantile", "kmeans", "pam")
mlab <- c(quantile = "Quantile", kmeans = "$k$-means", pam = "PAM")
shapes <- c("smooth", "jump", "sin")

# sign = TRUE: always show the sign of a difference; a value that rounds to 0 keeps
# its sign ($-$0.000 / $+$0.000) so a CI bound just below/above 0 stays readable.
num <- function(x, k, sign = FALSE) {
  s <- formatC(abs(x), format = "f", digits = k)
  tiny <- round(abs(x), k) == 0
  pre <- ifelse(x < 0 & (!tiny | sign), "$-$", ifelse(x > 0 & sign, "$+$", ""))
  ifelse(is.finite(x), paste0(pre, s), "--")
}
bold_if <- function(txt, flag) ifelse(flag, paste0("\\textbf{", txt, "}"), txt)
# significance marker vs Quantile: better / worse / none
marker <- function(better, worse) if (isTRUE(better)) "$^{+}$" else if (isTRUE(worse)) "$^{-}$" else ""
diff_cell <- function(m, lo, hi, k, pct = FALSE) {
  if (!is.finite(lo)) return(paste0(num(m, k, TRUE), " {\\footnotesize[--]}"))
  star <- if (lo > 0 || hi < 0) "$^{*}$" else ""
  paste0(num(m, k, TRUE), star, " {\\footnotesize[", num(lo, k, TRUE), ", ", num(hi, k, TRUE), "]}")
}
write_table <- function(file, caption, label, align, head, rows_blocks) {
  # [H] (float package) keeps each table under its heading
  lines <- c("\\begin{table}[H]", "\\centering", paste0("\\caption{", caption, "}"),
             paste0("\\label{", label, "}"),
             "\\resizebox{\\ifdim\\width>\\textwidth\\textwidth\\else\\width\\fi}{!}{%",  # shrink only if too wide
             paste0("\\begin{tabular}{", align, "}"), "\\toprule", head, "\\midrule")
  for (b in seq_along(rows_blocks)) {
    lines <- c(lines, paste0(vapply(rows_blocks[[b]], paste, character(1), collapse = " & "), " \\\\"))
    if (b < length(rows_blocks)) lines <- c(lines, "\\midrule")
  }
  writeLines(c(lines, "\\bottomrule", "\\end{tabular}}", "\\end{table}"), file.path(out, "tables", file))
  paste0("tables/", sub("\\.tex$", "", file))
}

# ================================ simulations ===================================
ss <- read.csv(file.path(tab_dir, "simulation_summary.csv"), stringsAsFactors = FALSE)
pd <- read.csv(file.path(tab_dir, "simulation_paired_differences.csv"), stringsAsFactors = FALSE)
pd$method <- sub("_minus_quantile", "", pd$metric)
pick <- function(df, ...) { f <- list(...); i <- rep(TRUE, nrow(df)); for (k in names(f)) i <- i & df[[k]] == f[[k]]; df[i, ] }
sim_caption <- function(g) {
  # Short standard caption; column definitions are in the text (Section 5).
  sp <- simulation_groups[simulation_groups$group == g, ]
  desc <- if (g == "regular") "regular design" else if (g == "skewed")
    "skewed design with independent variables" else if (is.na(sp$skewness))
    sprintf("correlated design, $corr=%s$", sp$correlation) else
    sprintf("skewed and correlated design, $skewness=%s$, $corr=%s$", sp$skewness, sp$correlation)
  paste0("Simulation results, ", desc, ". Bold: highest mean; $^{+}$/$^{-}$: significantly higher/lower ",
         "test $R^2$ than Quantile. Columns are described in Section~\\ref{sec:sim}.")
}
sim_input <- character()
for (g in simulation_groups$group) {
  blocks <- lapply(c(1, 3, 5), function(d) lapply(shapes, function(sh) {
    r2 <- sapply(methods, function(m) pick(ss, group = g, dimension = d, scenario = sh, method = m, metric = "r2"), simplify = FALSE)
    means <- sapply(r2, `[[`, "mean")
    cells <- mapply(function(r, best, m) {
      z <- if (m == "quantile") NULL else pick(pd, group = g, dimension = d, scenario = sh, method = m)
      paste0(bold_if(num(r$mean, 3), best), marker(z$lower_95 > 0, z$upper_95 < 0),
             if (style == "ci") paste0(" {\\footnotesize[", num(r$lower_95, 3), ", ", num(r$upper_95, 3), "]}")
             else paste0(" (", num(r$sd, 3), ")"))
    }, r2, means == max(means), methods)
    tb <- sapply(methods, function(m) num(pick(ss, group = g, dimension = d, scenario = sh, method = m,
                                               metric = "total_bases")$mean, 1))
    c(if (sh == "smooth") d else "", sh, cells, tb)
  }))
  sim_input[g] <- write_table(paste0("sim_", g, ".tex"), sim_caption(g), paste0("tab:sim-", g),
    "rlcccccc",
    c(paste0(" & & \\multicolumn{3}{c}{Avg. test $R^2$ ", if (style == "ci") "[95\\% CI]" else "(SD)", "} & ",
             "\\multicolumn{3}{c}{Avg. total number of bases} \\\\"),
      "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
      "$p$ & Scenario & Quantile & $k$-means & PAM & Quantile & $k$-means & PAM \\\\"),
    blocks)
}

# overview over the 9 settings of each design (settings independent; Welch df)
ov <- do.call(rbind, lapply(split(pd, list(pd$group, pd$method)), function(z) {
  v <- z$sd^2 / z$n_repeats; h <- qt(0.975, sum(v)^2 / sum(v^2 / (z$n_repeats - 1))) * sqrt(sum(v)) / nrow(z)
  data.frame(group = z$group[1], method = z$method[1], mean = mean(z$mean), lo = mean(z$mean) - h,
             hi = mean(z$mean) + h, better = sum(z$lower_95 > 0), worse = sum(z$upper_95 < 0))
}))
ov_blocks <- lapply(list(c("regular", "corr_0.1", "corr_0.5", "corr_0.9"), "skewed",
                         grep("^gskew", simulation_groups$group, value = TRUE)), function(gs) lapply(gs, function(g) {
  sp <- simulation_groups[simulation_groups$group == g, ]
  c(if (g == "regular") "Regular" else if (g == "skewed") "Skewed (indep.)" else if (is.na(sp$skewness)) "Correlated" else "Skewed + corr.",
    if (is.na(sp$skewness)) (if (g == "skewed") "2" else "0") else sp$skewness,
    if (is.na(sp$correlation)) "0" else sp$correlation,
    unlist(lapply(c("kmeans", "pam"), function(m) { r <- ov[ov$group == g & ov$method == m, ]
      c(diff_cell(r$mean, r$lo, r$hi, 3), paste0(r$better, " / ", r$worse)) })))
}))
sim_overview <- write_table("sim_overview.tex",
  paste0("Simulation results by design: difference in test $R^2$ from Quantile, averaged over the nine ",
         "settings of each design. $^{*}$: 95\\% CI excludes 0. Columns are described in ",
         "Section~\\ref{sec:sim-summary}."),
  "tab:sim-overview", "lllcccc",
  c(" & & & \\multicolumn{2}{c}{$k$-means $-$ Quantile} & \\multicolumn{2}{c}{PAM $-$ Quantile} \\\\",
    "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
    "Design & Skewness & Corr. & Avg. $\\Delta R^2$ [95\\% CI] & Better / worse & Avg. $\\Delta R^2$ [95\\% CI] & Better / worse \\\\"),
  ov_blocks)
fd <- merge(ov[ov$group != "skewed", ], simulation_groups, by = "group")
fd$skew_x <- ifelse(is.na(fd$skewness), 0, fd$skewness)
fd$corr_lab <- ifelse(is.na(fd$correlation), "independent", paste("latent corr", fd$correlation))
fd$method_lab <- factor(fd$method, c("kmeans", "pam"), c("k-means - Quantile", "PAM - Quantile"))
g1 <- ggplot(fd, aes(skew_x, mean, colour = corr_lab, group = corr_lab)) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(position = position_dodge(width = 0.15)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.25, position = position_dodge(width = 0.15)) +
  facet_wrap(~ method_lab) + labs(x = "Covariate skewness (0 = uniform margins)",
  y = expression("Average paired " * Delta * R^2 * " (95% CI)"), colour = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "top")
ggsave(file.path(out, "figures", "sim_delta_r2.pdf"), g1, width = 8, height = 3.6)

# ================================= real data ====================================
rs <- read.csv(file.path(tab_dir, "real_data_summary.csv"), stringsAsFactors = FALSE)
datasets <- c("cpu", "automobile", "communities_crime", "pima", "oecdpanel",
              "bike_sharing", "liver_disorders", "laheart", "real_estate")
dlab <- c(cpu = "cpu", automobile = "Automobile", communities_crime = "Communities and Crime",
          pima = "Pima", oecdpanel = "oecdpanel", bike_sharing = "Bike Sharing",
          liver_disorders = "Liver Disorders", laheart = "laheart", real_estate = "Real Estate")
fold_files <- list.files("results/real_data", "_repeat_[0-9]+\\.csv$", full.names = TRUE)
folds <- if (length(fold_files)) do.call(rbind, lapply(fold_files, read.csv, stringsAsFactors = FALSE)) else NULL

# corrected repeated-CV t-interval for mean(fold MSE_method - fold MSE_quantile), in % of quantile MSE
corrected_ci <- function(ds, kn, m) {
  if (is.null(folds)) return(c(NA, NA, NA))
  z <- folds[folds$dataset == ds & folds$knot_start == kn, ]
  w <- reshape(z[, c("repeat_id", "fold", "method", "mse", "n_train", "n_valid")],
               idvar = c("repeat_id", "fold"), timevar = "method", direction = "wide")
  if (nrow(w) < 100 || any(!is.finite(w[[paste0("mse.", m)]]))) return(c(NA, NA, NA))
  x <- w[[paste0("mse.", m)]] - w$mse.quantile; J <- length(x)
  ratio <- mean(w$n_valid.quantile) / mean(w$n_train.quantile)
  h <- qt(0.975, J - 1) * sqrt((1 / J + ratio) * var(x)); q <- mean(w$mse.quantile)
  100 * c(mean(x), mean(x) - h, mean(x) + h) / q
}
digits_for <- function(x) max(0, min(4, 3 - floor(log10(median(abs(x))))))
real_input <- character(); rel_rows <- list()
for (ds in datasets) {
  z <- rs[rs$dataset == ds, ]; k <- digits_for(z$mean[z$metric == "mse"])
  rows <- lapply(sort(unique(z$knot_start)), function(kn) {
    mse <- sapply(methods, function(m) pick(z, knot_start = kn, method = m, metric = "mse"), simplify = FALSE)
    means <- sapply(mse, `[[`, "mean")
    cis <- lapply(c(quantile = "quantile", kmeans = "kmeans", pam = "pam"), function(m) {
      if (m == "quantile") return(c(NA, NA, NA))
      ci <- corrected_ci(ds, kn, m)
      if (!is.finite(ci[1])) ci[1] <- 100 * (means[[m]] - means[["quantile"]]) / means[["quantile"]]
      rel_rows[[length(rel_rows) + 1]] <<- data.frame(dataset = ds, knot_start = kn, method = m,
                                                     rel = ci[1], lo = ci[2], hi = ci[3])
      ci })
    cells <- mapply(function(r, best, ci) paste0(bold_if(num(r$mean, k), best),
                                                 marker(ci[3] < 0, ci[2] > 0), " (", num(r$sd, k), ")"),
                    mse, means == min(means), cis)
    tb <- sapply(methods, function(m) num(pick(z, knot_start = kn, method = m, metric = "total_bases")$mean, 1))
    c(paste0(kn, ", 10, 5"), cells, tb)
  })
  real_input[ds] <- write_table(paste0("real_", ds, ".tex"),
    paste0("Real data results, ", tolower(dlab[[ds]]), ". Bold: lowest mean; $^{+}$/$^{-}$: significantly ",
           "lower/higher MSE than Quantile. Columns are described in Section~\\ref{sec:real}."),
    paste0("tab:real-", ds), "lcccccc",
    c(paste0(" & \\multicolumn{3}{c}{Avg. MSE (SD)} & ",
             "\\multicolumn{3}{c}{Avg. total number of bases} \\\\"),
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
      "Number of knots & Quantile & $k$-means & PAM & Quantile & $k$-means & PAM \\\\"),
    list(rows))
}
rel <- do.call(rbind, rel_rows)

skew3 <- function(x) { n <- length(x); x <- x - mean(x); (mean(x^3) / mean(x^2)^1.5) * ((n - 1) / n)^1.5 }
chars <- do.call(rbind, lapply(datasets, function(ds) {
  dat <- load_dataset(ds); X <- dat$X[, apply(dat$X, 2, sd) > 0, drop = FALSE]
  cont <- apply(X, 2, function(v) length(unique(v)) > 2)
  sk <- apply(X[, cont, drop = FALSE], 2, skew3); cm <- cor(X[, cont, drop = FALSE]); cm <- cm[upper.tri(cm)]
  data.frame(dataset = ds, n = nrow(X), d = ncol(X), d_cont = sum(cont), avg_skew = mean(sk),
             pct_skew = 100 * mean(abs(sk) > 0.5), pct_corr = 100 * mean(abs(cm) > 0.5))
}))
chars <- chars[order(-chars$avg_skew), ]
real_overview <- write_table("real_overview.tex",
  paste0("Real data sets: characteristics of the continuous predictors and the CV-MSE difference vs Quantile ",
         "(\\% of the Quantile MSE) averaged over the 9 knot settings. Sig. better / worse: number of knot ",
         "settings whose corrected 95\\% CI lies below / above 0 (-- if fold-level results were not available)."),
  "tab:real-overview", "lrrrrrcccc",
  c(paste0(" & & & \\multicolumn{3}{c}{Continuous predictors} & \\multicolumn{2}{c}{$k$-means vs Quantile} & ",
           "\\multicolumn{2}{c}{PAM vs Quantile} \\\\"),
    "\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}",
    paste0("Data set & $n$ & $d$ (cont.) & Avg. skewness & \\% $|$skew$|>0.5$ & \\% $|$corr$|>0.5$ & ",
           "Avg. diff. & Sig. better / worse & Avg. diff. & Sig. better / worse \\\\")),
  list(lapply(chars$dataset, function(ds) { a <- chars[chars$dataset == ds, ]
    c(dlab[[ds]], a$n, paste0(a$d, " (", a$d_cont, ")"), num(a$avg_skew, 2), num(a$pct_skew, 0), num(a$pct_corr, 0),
      unlist(lapply(c("kmeans", "pam"), function(m) { r <- rel[rel$dataset == ds & rel$method == m, ]
        c(paste0(num(mean(r$rel), 1, TRUE), "\\%"),
          if (all(is.finite(r$lo))) paste0(sum(r$hi < 0), " / ", sum(r$lo > 0)) else "--") }))) })))
# Across-dataset comparison (Demsar 2006): one number per data set.
across_dataset <- function(skew, diffs) {
  do.call(rbind, lapply(names(diffs), function(m) {
    x <- diffs[[m]]
    wt <- suppressWarnings(wilcox.test(x, mu = 0, exact = TRUE))
    bt <- binom.test(sum(x < 0), sum(x != 0))
    ct <- suppressWarnings(cor.test(skew, x, method = "spearman", exact = TRUE))
    data.frame(method = m, better = sum(x < 0), n = length(x), median = median(x),
               wil_p = wt$p.value, sign_p = bt$p.value, rho = unname(ct$estimate), rho_p = ct$p.value)
  }))
}
across_table_tex <- function(res) {
  f2 <- function(x, k) formatC(x, format = "f", digits = k)
  sgn <- function(x, k) paste0(ifelse(x < 0, "$-$", "$+$"), f2(abs(x), k))
  body <- vapply(seq_len(nrow(res)), function(i) { r <- res[i, ]
    paste0(c(Quantile = "", kmeans = "$k$-means", pam = "PAM")[[r$method]], " & ", r$better, " of ", r$n, " & ",
           sgn(r$median, 1), "\\% & ", f2(r$wil_p, 2), " & ", f2(r$sign_p, 2), " & ", sgn(r$rho, 2), " & ",
           f2(r$rho_p, 2), " \\\\") }, character(1))
  c("\\begin{table}[H]", "\\centering",
    paste0("\\caption{Real data results across datasets: difference in cross-validated MSE from Quantile. ",
           "Columns are described in Section~\\ref{sec:real-summary}.}"),
    "\\label{tab:real-across}",
    "\\begin{tabular}{lcccccc}", "\\toprule",
    " & Data sets with & Median & Wilcoxon & Sign test & \\multicolumn{2}{c}{Spearman (skewness, diff.)} \\\\",
    "\\cmidrule(lr){6-7}",
    "Method & lower MSE & difference & $p$ & $p$ & $\\rho$ & $p$ \\\\", "\\midrule", body,
    "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

acr <- across_dataset(chars$avg_skew, lapply(c(kmeans = "kmeans", pam = "pam"), function(m)
  sapply(chars$dataset, function(ds) mean(rel$rel[rel$dataset == ds & rel$method == m]))))
writeLines(across_table_tex(acr), file.path(out, "tables", "real_across.tex"))
rel$label <- factor(dlab[rel$dataset], levels = rev(dlab[chars$dataset]))
rel$method_lab <- factor(rel$method, c("kmeans", "pam"), c("k-means vs Quantile", "PAM vs Quantile"))
g2 <- ggplot(rel, aes(rel, label)) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(aes(colour = knot_start), size = 1.4, position = position_jitter(height = 0.15, seed = 1)) +
  stat_summary(fun = mean, geom = "point", shape = 124, size = 6) + facet_wrap(~ method_lab) +
  scale_colour_viridis_c(name = "First-order\nknots") +
  labs(x = "CV-MSE difference vs Quantile (% of Quantile MSE)", y = NULL) + theme_bw(base_size = 10)
ggsave(file.path(out, "figures", "real_rel_mse.pdf"), g2, width = 8, height = 3.8)

# ================================== main.tex ====================================
sec <- function(title, inputs) c(title, paste0("\\input{", inputs, "}"), "\\clearpage", "")
sub_skew <- unlist(lapply(c(0.5, 1, 2, 3), function(s) {
  gs <- simulation_groups$group[!is.na(simulation_groups$skewness) & simulation_groups$skewness == s]
  c(sprintf("\\subsubsection{Skewness %s}", s), paste0("\\input{", sim_input[gs], "}"), "\\clearpage")
}))
main <- c(
  "\\documentclass{article}", "\\usepackage{graphicx} % Required for inserting images",
  "\\usepackage{booktabs}", "\\usepackage{float}", "\\usepackage[margin=2cm]{geometry}", "",
  "\\title{hal results}", "\\date{November 2026}", "", "\\begin{document}", "", "\\maketitle", "",
  "\\section{Notes on the tables}",
  paste0("Simulations: every replicate draws new training ($n=1000$) and test data. ",
         if (style == "ci") "Brackets give 95\\% $t$-intervals " else "Parentheses give SDs ",
         "over the 10 replicates. Because all methods are fitted to the same data, a $k$-means or PAM ",
         "result is marked $^{+}$ ($^{-}$) when its paired difference from Quantile has a 95\\% CI above ",
         "(below) 0; overlapping intervals of two methods do not imply that they perform the same. ",
         "Real data: parentheses give the SD over 10 repeats of 10-fold CV (sensitivity to the fold split); ",
         "$^{+}$/$^{-}$ use the corrected repeated-CV $t$-test (Nadeau and Bengio, 2003), which accounts ",
         "for overlapping training sets."),
  "", "\\section{Simulations}", "\\subsection{Overview}", "\\input{tables/sim_overview}",
  "\\begin{figure}[ht]", "\\centering", "\\includegraphics[width=\\textwidth]{figures/sim_delta_r2.pdf}",
  "\\caption{Average paired difference in test $R^2$ vs Quantile, by covariate skewness and latent correlation.}",
  "\\label{fig:sim-delta}", "\\end{figure}", "\\clearpage", "",
  "\\subsection{Regular}", paste0("\\input{", sim_input["regular"], "}"), "",
  "\\subsection{Correlated variables}", paste0("\\input{", sim_input[c("corr_0.1", "corr_0.5", "corr_0.9")], "}"), "\\clearpage", "",
  "\\subsection{Skewed variables}", paste0("\\input{", sim_input["skewed"], "}"), "",
  "\\subsection{Skewed and correlated variables}", sub_skew, "",
  "\\section{Real World Data}", "\\subsection{Overview}", paste0("\\input{", real_overview, "}"),
  "\\input{tables/real_across}",
  "\\begin{figure}[ht]", "\\centering", "\\includegraphics[width=\\textwidth]{figures/real_rel_mse.pdf}",
  "\\caption{CV-MSE difference vs Quantile for each knot setting (dots) and their mean (bar).}",
  "\\label{fig:real-rel}", "\\end{figure}", "\\clearpage", "",
  unlist(lapply(datasets, function(ds) c(sprintf("\\subsection{%s}", dlab[[ds]]),
                                         paste0("\\input{", real_input[ds], "}"), ""))),
  "\\end{document}")
writeLines(main, file.path(out, "main.tex"))
message("Overleaf project written to ", out, " (", length(sim_input) + length(real_input) + 2, " tables)")
