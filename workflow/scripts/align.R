# Align one paired-end library to GRCh38 + ERCC with the Subread aligner.
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages(library(Rsubread))

# unique = TRUE: only read pairs with a single best mapping location are reported, so reads
# from repeats and multi-copy genes never reach the counting step.
stats <- align(
  index = snakemake@params[["index"]],
  readfile1 = snakemake@input[["r1"]],
  readfile2 = snakemake@input[["r2"]],
  type = "rna",
  unique = TRUE,
  output_file = snakemake@output[["bam"]],
  nthreads = snakemake@threads
)

write.table(
  data.frame(metric = rownames(stats), value = stats[[1]]),
  snakemake@output[["stats"]],
  sep = "\t", quote = FALSE, row.names = FALSE
)
