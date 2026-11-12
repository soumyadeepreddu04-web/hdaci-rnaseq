# Adapter and quality trimming of one paired-end library with fastp.
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")

suppressPackageStartupMessages(library(Rfastp))

# rfastp() appends _R1.fastq.gz / _R2.fastq.gz / .json / .html to this prefix
prefix <- sub("_R1\\.fastq\\.gz$", "", snakemake@output[["r1"]])

rfastp(
  read1 = snakemake@input[["r1"]],
  read2 = snakemake@input[["r2"]],
  outputFastq = prefix,
  thread = snakemake@threads,
  verbose = FALSE
)

stopifnot(
  file.rename(paste0(prefix, ".json"), snakemake@output[["json"]]),
  file.rename(paste0(prefix, ".html"), snakemake@output[["html"]])
)
