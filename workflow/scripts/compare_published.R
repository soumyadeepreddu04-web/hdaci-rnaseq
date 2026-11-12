# External check against the count table deposited by the original authors.
#
# Their counts come from a different toolchain (TopHat + HTSeq, GENCODE v21) and from every
# read in each library. Fitting the same DESeq2 model to their counts and comparing the fold
# changes shows how much of the result depends on the aligner, the annotation release and
# the sequencing depth used here.
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})
source(file.path(snakemake@scriptdir, "theme.R"))

cfg <- snakemake@params[["de"]]
alpha <- cfg[["alpha"]]
drugs <- c(SAHA_vs_DMSO = "SAHA", RMD_vs_DMSO = "RMD")
unversioned <- function(id) sub("\\.[0-9]+$", "", id)

samples <- fread(snakemake@input[["samples"]])[, .(sample, donor, treatment)]

# ---- authors' counts: columns are named like "8_D2_uninf_SAHA"
published <- fread(snakemake@input[["published"]])
setnames(published, 1, "gene_id")
published <- published[grepl("^ENSG[0-9]", gene_id)]
columns <- grep("_uninf_", names(published), value = TRUE)
pub <- as.matrix(published[, columns, with = FALSE])
colnames(pub) <- sub("^[0-9]+_(D[0-9]+)_uninf_(.+)$", "\\1_\\2", columns)
rownames(pub) <- unversioned(published$gene_id)
pub <- pub[, samples$sample]

# ---- this workflow's counts
counts <- fread(snakemake@input[["counts"]])[grepl("^ENSG[0-9]", gene_id) & !grepl("_PAR_Y$", gene_id)]
mine <- as.matrix(counts[, samples$sample, with = FALSE])
rownames(mine) <- unversioned(counts$gene_id)

# ---- 1. do the two pipelines agree on expression levels, library by library?
shared <- intersect(rownames(mine), rownames(pub))
expressed <- shared[rowSums(mine[shared, ] >= 10) >= 4 | rowSums(pub[shared, ] >= 10) >= 4]
log_cpm <- function(m) log2(t(t(m) / colSums(m)) * 1e6 + 1)
a <- log_cpm(mine[shared, ])[expressed, ]
b <- log_cpm(pub[shared, ])[expressed, ]
per_sample <- data.table(
  sample = samples$sample,
  pearson = vapply(samples$sample, function(s) cor(a[, s], b[, s]), 0),
  spearman = vapply(samples$sample, function(s) cor(a[, s], b[, s], method = "spearman"), 0)
)

# ---- 2. same model, authors' counts
coldata <- data.frame(
  row.names = samples$sample,
  donor = factor(samples$donor),
  treatment = factor(samples$treatment, levels = c("DMSO", "SAHA", "RMD"))
)
dds <- DESeqDataSetFromMatrix(pub, coldata, design = ~ donor + treatment)
dds <- DESeq(dds[rowSums(counts(dds) >= cfg[["min_count"]]) >= cfg[["min_samples"]], ])

files <- unlist(snakemake@input[["tables"]])
names(files) <- sub("\\.csv\\.gz$", "", basename(files))

comparison <- list()
stats <- list()
for (name in names(drugs)) {
  res <- results(dds, contrast = c("treatment", drugs[[name]], "DMSO"), alpha = alpha)
  theirs <- data.table(ensembl_id = rownames(res), lfc_published = res$log2FoldChange, padj_published = res$padj)
  ours <- fread(files[[name]])[, .(ensembl_id = unversioned(gene_id), gene_name,
                                   lfc_this_workflow = log2FoldChange, padj_this_workflow = padj)]
  both <- merge(ours, theirs, by = "ensembl_id")[!is.na(lfc_this_workflow) & !is.na(lfc_published)]
  both[, contrast := name]
  comparison[[name]] <- both

  sig_ours <- both[!is.na(padj_this_workflow) & padj_this_workflow < alpha]
  sig_theirs <- both[!is.na(padj_published) & padj_published < alpha]
  confirmed <- sig_ours[!is.na(padj_published) & padj_published < alpha &
                          sign(lfc_this_workflow) == sign(lfc_published)]
  stats[[name]] <- list(
    genes_compared = nrow(both),
    lfc_pearson_r = both[, cor(lfc_this_workflow, lfc_published)],
    significant_this_workflow = nrow(sig_ours),
    significant_published_counts = nrow(sig_theirs),
    published_counts_up = sum(res$padj < alpha & res$log2FoldChange > 0, na.rm = TRUE),
    published_counts_down = sum(res$padj < alpha & res$log2FoldChange < 0, na.rm = TRUE),
    confirmed = nrow(confirmed),
    confirmed_fraction = nrow(confirmed) / nrow(sig_ours),
    opposite_direction = sig_ours[!is.na(padj_published) & padj_published < alpha &
                                    sign(lfc_this_workflow) != sign(lfc_published), .N]
  )
}
comparison <- rbindlist(comparison)
fwrite(comparison, snakemake@output[["table"]])

write_json(list(
  genes_shared_between_annotations = length(shared),
  genes_expressed = length(expressed),
  per_sample_log_cpm_pearson = list(min = min(per_sample$pearson), median = median(per_sample$pearson)),
  per_sample_log_cpm_spearman = list(min = min(per_sample$spearman), median = median(per_sample$spearman)),
  published_counts_genes_tested = nrow(dds),
  contrasts = stats,
  # edgeR + RUVSeq at FDR < 0.05 in mock-infected cells, as stated in the GEO record of GSE114883
  reported_by_authors = list(
    genes_after_filtering = 16058,
    SAHA_vs_DMSO = list(up = 3971, down = 2940),
    RMD_vs_DMSO = list(up = 5068, down = 4050)
  )
), snakemake@output[["summary"]], auto_unbox = TRUE, pretty = TRUE, digits = 4)

# ---- figure
comparison[, contrast := factor(contrast, levels = names(drugs))]
limit <- max(abs(quantile(c(comparison$lfc_this_workflow, comparison$lfc_published), c(0.001, 0.999))))
r_labels <- comparison[, .(text = sprintf("r = %.3f\n%s genes", cor(lfc_this_workflow, lfc_published),
                                          format(.N, big.mark = ","))), by = contrast]

plot <- ggplot(comparison, aes(lfc_published, lfc_this_workflow)) +
  geom_bin_2d(bins = 110) +
  geom_abline(slope = 1, intercept = 0, colour = col$muted, linewidth = 0.4) +
  geom_text(data = r_labels, aes(x = -limit, y = limit, label = text), hjust = 0, vjust = 1,
            colour = col$ink2, size = 3) +
  facet_wrap(~contrast, labeller = as_labeller(contrast_label)) +
  scale_fill_gradientn(colours = col$sequential, trans = "log10", name = "genes per bin") +
  coord_equal(xlim = c(-limit, limit), ylim = c(-limit, limit)) +
  labs(x = "log2 fold change from the authors' counts (TopHat + HTSeq, all reads)",
       y = "log2 fold change, this workflow") +
  theme_report() +
  theme(legend.position = "right", legend.justification = "center")
save_png(plot, snakemake@output[["figure"]], 8.4, 4.2)
