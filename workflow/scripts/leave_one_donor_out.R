# Does the result depend on any one donor?
#
# The model is refitted four times, each time without one donor, and the drug-vs-DMSO fold
# changes and significant genes are compared with the full four-donor fit.
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

alpha <- snakemake@params[["de"]][["alpha"]]
drugs <- c(SAHA_vs_DMSO = "SAHA", RMD_vs_DMSO = "RMD")

full <- readRDS(snakemake@input[["dds"]])

rows <- list()
for (left_out in levels(full$donor)) {
  subset <- full[, full$donor != left_out]
  subset$donor <- droplevels(subset$donor)
  subset <- DESeq(subset, quiet = TRUE)

  for (name in names(drugs)) {
    contrast <- c("treatment", drugs[[name]], "DMSO")
    a <- results(full, contrast = contrast, alpha = alpha)
    b <- results(subset, contrast = contrast, alpha = alpha)
    sig_full <- which(a$padj < alpha)
    sig_subset <- which(b$padj < alpha)
    kept <- intersect(sig_full, sig_subset)
    rows[[length(rows) + 1]] <- data.table(
      contrast = name,
      donor_left_out = left_out,
      lfc_pearson_r = cor(a$log2FoldChange, b$log2FoldChange, use = "complete.obs"),
      significant_full = length(sig_full),
      significant_without_donor = length(sig_subset),
      still_significant_fraction = length(kept) / length(sig_full),
      direction_changes = sum(sign(a$log2FoldChange[kept]) != sign(b$log2FoldChange[kept]))
    )
  }
}

fwrite(rbindlist(rows), snakemake@output[[1]])
