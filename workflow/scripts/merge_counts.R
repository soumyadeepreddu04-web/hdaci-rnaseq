# Combine per-sample counts into one matrix and collect per-sample QC metrics.
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

samples <- unlist(snakemake@params[["samples"]])
named <- function(files) setNames(unlist(files), samples)  # inputs are in sample-sheet order

# ---- count matrix: gene_id, length, one column per sample
counts <- lapply(named(snakemake@input[["counts"]]), fread)
matrix <- counts[[1]][, .(gene_id, length)]
for (s in samples) {
  stopifnot(identical(counts[[s]]$gene_id, matrix$gene_id))
  set(matrix, j = s, value = counts[[s]]$count)
}
fwrite(matrix, snakemake@output[["counts"]], sep = "\t")

# ---- QC table
fastp <- lapply(named(snakemake@input[["fastp"]]), fromJSON)
align <- lapply(named(snakemake@input[["stats"]]), function(f) {
  x <- fread(f)
  setNames(x$value, x$metric)
})
assigned <- lapply(named(snakemake@input[["summaries"]]), function(f) {
  x <- fread(f)
  setNames(x[[2]], x$Status)
})
is_ercc <- grepl("^ERCC-", matrix$gene_id)

qc <- rbindlist(lapply(samples, function(s) {
  before <- fastp[[s]]$summary$before_filtering
  after <- fastp[[s]]$summary$after_filtering
  a <- align[[s]]
  n_assigned <- assigned[[s]][["Assigned"]]
  data.table(
    sample = s,
    read_pairs = as.integer(before$total_reads / 2),
    q30_rate = before$q30_rate,
    pairs_after_trimming = as.integer(after$total_reads / 2),
    duplication_rate = fastp[[s]]$duplication$rate,
    uniquely_mapped_fraction = a[["Uniquely_mapped_fragments"]] / a[["Total_fragments"]],
    properly_paired_fraction = a[["Properly_paired_fragments"]] / a[["Total_fragments"]],
    assigned_pairs = as.integer(n_assigned),
    assigned_fraction = n_assigned / a[["Total_fragments"]],
    ercc_fraction_of_assigned = sum(matrix[[s]][is_ercc]) / n_assigned,
    genes_detected = sum(matrix[[s]][!is_ercc] >= 10),
    no_feature_fraction = assigned[[s]][["Unassigned_NoFeatures"]] / a[["Total_fragments"]]
  )
}))
fwrite(qc, snakemake@output[["qc"]], sep = "\t", scipen = 50)

# ---- alignment table for MultiQC (custom content)
header <- c(
  "# id: 'subread_align'",
  "# section_name: 'Subread alignment'",
  "# description: 'Read pairs aligned to GRCh38 + ERCC with the Subread aligner (Rsubread::align).'",
  "# plot_type: 'table'"
)
mqc <- qc[, .(
  Sample = sample,
  `Read pairs aligned` = pairs_after_trimming,
  `% uniquely mapped` = round(100 * uniquely_mapped_fraction, 2),
  `% properly paired` = round(100 * properly_paired_fraction, 2)
)]
writeLines(header, snakemake@output[["mqc"]])
fwrite(mqc, snakemake@output[["mqc"]], sep = "\t", append = TRUE, col.names = TRUE, scipen = 50)
