# Differential expression with DESeq2.
#
# Each donor contributes one DMSO, one SAHA and one RMD library, so donor is a blocking factor:
# the model ~ donor + treatment compares treatments within donors.
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(jsonlite)
})

cfg <- snakemake@params[["de"]]
alpha <- cfg[["alpha"]]
lfc_strong <- cfg[["lfc_strong"]]

samples <- fread(snakemake@input[["samples"]])[, .(sample, donor, treatment)]
counts <- fread(snakemake@input[["counts"]])
genes <- fread(snakemake@input[["genes"]])

mat <- as.matrix(counts[, samples$sample, with = FALSE])
rownames(mat) <- counts$gene_id
mat <- mat[!grepl("^ERCC-", rownames(mat)), ]  # spike-ins are analysed separately (ercc.R)

treatments <- c(cfg[["reference_level"]], setdiff(unique(samples$treatment), cfg[["reference_level"]]))
coldata <- data.frame(
  row.names = samples$sample,
  donor = factor(samples$donor),
  treatment = factor(samples$treatment, levels = treatments)
)

dds <- DESeqDataSetFromMatrix(mat, coldata, design = ~ donor + treatment)
keep <- rowSums(counts(dds) >= cfg[["min_count"]]) >= cfg[["min_samples"]]
dds <- DESeq(dds[keep, ])

# Contrasts are named "<numerator>_vs_<denominator>" in the Snakefile (CONTRASTS)
contrast_names <- sub("\\.csv\\.gz$", "", basename(unlist(snakemake@output[["tables"]])))
contrasts <- lapply(strsplit(contrast_names, "_vs_"), function(x) c("treatment", x[1], x[2]))
names(contrasts) <- contrast_names

summary <- list(
  genes_in_annotation = nrow(mat),
  genes_tested = sum(keep),
  alpha = alpha,
  lfc_strong = lfc_strong,
  size_factors = as.list(round(sizeFactors(dds), 3)),
  contrasts = list()
)

for (i in seq_along(contrasts)) {
  name <- names(contrasts)[i]
  res <- results(dds, contrast = contrasts[[i]], alpha = alpha)
  # Shrinkage pulls the fold changes of low-count, noisy genes towards zero. The shrunken
  # values are for ranking and plotting; p-values come from the unshrunken Wald test.
  shrunk <- lfcShrink(dds, contrast = contrasts[[i]], res = res, type = "ashr", quiet = TRUE)

  table <- data.table(
    gene_id = rownames(res),
    baseMean = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcSE = res$lfcSE,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj,
    log2FoldChange_shrunk = shrunk$log2FoldChange
  )
  table <- merge(genes[, .(gene_id, gene_name, gene_type)], table, by = "gene_id")
  setorder(table, padj, -baseMean, na.last = TRUE)
  fwrite(table, snakemake@output[["tables"]][i])

  sig <- table[!is.na(padj) & padj < alpha]
  summary$contrasts[[name]] <- list(
    up = sig[log2FoldChange > 0, .N],
    down = sig[log2FoldChange < 0, .N],
    up_strong = sig[log2FoldChange_shrunk >= lfc_strong, .N],
    down_strong = sig[log2FoldChange_shrunk <= -lfc_strong, .N],
    top_up = sig[log2FoldChange > 0][seq_len(min(.N, 10)), gene_name],
    top_down = sig[log2FoldChange < 0][seq_len(min(.N, 10)), gene_name],
    highlight = table[gene_name %in% cfg[["highlight_genes"]],
                      .(gene_name, log2FoldChange = round(log2FoldChange, 3), padj = signif(padj, 3))]
  )
}

vsd <- vst(dds, blind = FALSE)
fwrite(as.data.table(assay(vsd), keep.rownames = "gene_id"), snakemake@output[["vst"]], sep = "\t")
fwrite(
  as.data.table(counts(dds, normalized = TRUE), keep.rownames = "gene_id"),
  snakemake@output[["normalised"]], sep = "\t"
)
saveRDS(dds, snakemake@output[["dds"]])
write_json(summary, snakemake@output[["summary"]], auto_unbox = TRUE, pretty = TRUE)
