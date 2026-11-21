#!/usr/bin/env Rscript
# Run this ONCE in an interactive compute session with internet access.
# Set HAL9001_REF to install a different commit or branch of the fork.
repo_lib <- file.path(getwd(), ".R", "library")
dir.create(repo_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(repo_lib, .libPaths()))
mirror <- "https://cloud.r-project.org"
for (pkg in c("remotes", "glmnet", "MASS", "rlang")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = mirror, lib = repo_lib)
  }
}
# Default: the commit of tq21/hal9001@kmeans-knot-screen used for the paper.
ref <- Sys.getenv("HAL9001_REF", "a81af69b28639221b4146ba8cfdcfcc83bb68aa2")
remote <- paste0("tq21/hal9001@", ref)
message("Installing custom HAL fork: ", remote)
remotes::install_github(remote, lib = repo_lib, upgrade = "never",
                        dependencies = TRUE)
message("Installed HAL fork; run scripts/check_environment.R")
