# Gene set enrichment analysis (fgsea) on the MSigDB Hallmark collection.
#
# Genes are ranked by the DESeq2 Wald statistic, which combines the direction and the
# confidence of each change, so no significance cut-off is needed to define a gene list.
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
})
source(file.path(snakemake@scriptdir, "theme.R"))

pathways <- gmtPathways(snakemake@input[["gmt"]])
files <- unlist(snakemake@input[["tables"]])
names(files) <- sub("\\.csv\\.gz$", "", basename(files))

run_gsea <- function(name) {
  res <- fread(files[[name]])[!is.na(stat) & gene_name != ""]
  # one statistic per gene symbol: keep the entry with the strongest signal
  res <- res[order(-abs(stat))][!duplicated(gene_name)]
  ranks <- sort(setNames(res$stat, res$gene_name), decreasing = TRUE)
  set.seed(1)
  out <- fgsea(pathways, ranks, minSize = 15, maxSize = 500, eps = 0,
               BPPARAM = BiocParallel::SerialParam())
  out[, .(contrast = name, pathway, size, NES, pval, padj,
          leading_edge = vapply(leadingEdge, function(g) paste(head(g, 15), collapse = ";"), ""))]
}

gsea <- rbindlist(lapply(names(files), run_gsea))
setorder(gsea, contrast, padj)
fwrite(gsea, snakemake@output[["table"]])

# ---- figure: Hallmark sets significant for at least one drug
drugs <- c("SAHA_vs_DMSO", "RMD_vs_DMSO")
shown <- gsea[contrast %in% drugs]
keep <- shown[padj < 0.05, unique(pathway)]
shown <- shown[pathway %in% keep]
shown[, label := gsub("_", " ", sub("^HALLMARK_", "", pathway))]
order_by <- shown[, .(mean_nes = mean(NES)), by = label][order(mean_nes), label]
shown[, `:=`(
  label = factor(label, levels = order_by),
  contrast = factor(contrast, levels = drugs),
  direction = fifelse(padj >= 0.05, "ns", fifelse(NES > 0, "up", "down"))
)]

plot <- ggplot(shown, aes(NES, label, fill = direction)) +
  geom_col(width = 0.62) +
  geom_vline(xintercept = 0, colour = col$axis, linewidth = 0.3) +
  facet_wrap(~contrast, labeller = as_labeller(contrast_label)) +
  scale_fill_manual(
    values = c(up = col$up, down = col$down, ns = col$neutral),
    breaks = c("up", "down", "ns"),
    labels = c("enriched among genes that go up", "enriched among genes that go down", "not significant (FDR ≥ 0.05)"),
    name = NULL
  ) +
  labs(x = "Normalised enrichment score", y = NULL) +
  theme_report(base_size = 9) +
  theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 6.5), legend.key.size = unit(0.3, "cm"))
save_png(plot, snakemake@output[["figure"]], 8, 0.9 + 0.2 * length(order_by))
