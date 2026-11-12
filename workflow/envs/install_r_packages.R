#!/usr/bin/env Rscript
# Install the R/Bioconductor packages the workflow needs.
#
#   Rscript workflow/envs/install_r_packages.R
#
# Packages go into the first entry of .libPaths(). To keep them out of your
# main R library, point R_LIBS_USER at a project folder first:
#
#   mkdir -p .Rlib && export R_LIBS_USER="$PWD/.Rlib"

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)

lib <- .libPaths()[1]
message("Installing into: ", lib)

cran <- c("BiocManager", "data.table", "R.utils", "ggplot2", "ggrepel", "patchwork", "ashr", "jsonlite")
bioc <- c("Rsubread", "Rfastp", "DESeq2", "fgsea")

need <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(need(cran)) > 0) install.packages(need(cran), lib = lib)
if (length(need(bioc)) > 0) BiocManager::install(need(bioc), lib = lib, update = FALSE, ask = FALSE)

missing <- need(c(cran, bioc))
if (length(missing) > 0) stop("Not installed: ", paste(missing, collapse = ", "))

message("Bioconductor ", BiocManager::version(), " on ", R.version.string)
for (p in c(cran, bioc)) message(sprintf("  %-12s %s", p, as.character(packageVersion(p))))
