# Build the Subread index for GRCh38 + ERCC.
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(Rsubread)

buildindex(
  basename = snakemake@params[["basename"]],
  reference = snakemake@input[["fasta"]],
  gappedIndex = isTRUE(snakemake@params[["gapped"]]),
  indexSplit = FALSE
)
