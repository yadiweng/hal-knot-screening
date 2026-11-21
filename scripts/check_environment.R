#!/usr/bin/env Rscript
source("R/common.R")
root <- project_root()
message("Project root: ", root)
message("R version: ", R.version.string)
message(".libPaths(): ", paste(.libPaths(), collapse = ", "))
for (pkg in c("hal9001", "glmnet", "MASS", "rlang", "remotes")) {
  message(pkg, ": ", if (requireNamespace(pkg, quietly = TRUE))
          as.character(packageVersion(pkg)) else "MISSING")
}
if (requireNamespace("hal9001", quietly = TRUE)) {
  f <- names(formals(hal9001::enumerate_basis))
  message("enumerate_basis arguments: ", paste(f, collapse = ", "))
  if (!"algorithm" %in% f) {
    stop("The installed hal9001 does not expose algorithm=; use the custom HAL fork")
  }
}
message("Input dataset count: ", length(list.files("data/raw", recursive = TRUE)))
writeLines(capture.output(sessionInfo()), "results/tables/sessionInfo.txt")
message("Saved results/tables/sessionInfo.txt")
message("Environment inspection complete; execute smoke tests to validate fitting")
