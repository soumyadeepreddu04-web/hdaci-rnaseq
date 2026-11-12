# QC and differential-expression figures.
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(jsonlite)
})
source(file.path(snakemake@scriptdir, "theme.R"))

cfg <- snakemake@params[["de"]]
alpha <- cfg[["alpha"]]

samples <- fread(snakemake@input[["samples"]])[, .(sample, donor, treatment)]
samples[, treatment := factor(treatment, levels = names(col$treatment))]
tables <- lapply(snakemake@input[["tables"]], fread)
names(tables) <- sub("\\.csv\\.gz$", "", basename(unlist(snakemake@input[["tables"]])))

vst <- fread(snakemake@input[["vst"]])
expr <- as.matrix(vst[, samples$sample, with = FALSE])
rownames(expr) <- vst$gene_id

# Expression with each donor's own mean removed, gene by gene. This is what the paired model
# compares, and it is used for display only.
within_donor <- expr
for (d in unique(samples$donor)) {
  idx <- samples[donor == d, sample]
  within_donor[, idx] <- expr[, idx] - rowMeans(expr[, idx])
}

# ---------------------------------------------------------------- PCA

pca_scores <- function(m, n_top = 500) {
  variable <- order(apply(m, 1, var), decreasing = TRUE)[seq_len(n_top)]
  p <- prcomp(t(m[variable, ]))
  explained <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  list(scores = cbind(samples, PC1 = p$x[, 1], PC2 = p$x[, 2]), explained = explained)
}

pca_panel <- function(m, title) {
  p <- pca_scores(m)
  ggplot(p$scores, aes(PC1, PC2, colour = treatment, shape = treatment)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = donor), colour = col$ink2, size = 2.8, show.legend = FALSE,
                    box.padding = 0.35, point.padding = 0.4, min.segment.length = Inf, seed = 1) +
    scale_colour_manual(values = col$treatment, name = NULL) +
    scale_shape_manual(values = shape_treatment, name = NULL) +
    labs(title = title,
         x = sprintf("PC1 (%s%% of variance)", p$explained[1]),
         y = sprintf("PC2 (%s%% of variance)", p$explained[2])) +
    theme_report()
}

pca <- pca_panel(expr, "All variation") +
  pca_panel(within_donor, "Donor differences removed") +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
save_png(pca, snakemake@output[["pca"]], 8, 4.2)

# ---------------------------------------------------------------- sample distances

distance <- as.matrix(dist(t(expr)))
ordering <- hclust(as.dist(distance), method = "average")$order
ordered <- rownames(distance)[ordering]
long <- as.data.table(as.table(distance))
setnames(long, c("a", "b", "distance"))
long[, `:=`(a = factor(a, levels = ordered), b = factor(b, levels = rev(ordered)))]

heat_distance <- ggplot(long, aes(a, b, fill = distance)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_gradientn(colours = col$sequential, name = "Euclidean distance\n(variance-stabilised counts)") +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_report() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        legend.position = "right", legend.justification = "center")
save_png(heat_distance, snakemake@output[["distance"]], 6.2, 4.8)

# ---------------------------------------------------------------- volcano plots

volcano_cap <- 150  # -log10(adjusted p) above this is drawn at the cap

volcano_data <- rbindlist(lapply(c("SAHA_vs_DMSO", "RMD_vs_DMSO"), function(name) {
  x <- tables[[name]][!is.na(padj)]
  x[, direction := fifelse(padj >= alpha, "ns", fifelse(log2FoldChange > 0, "up", "down"))]
  x[, y := pmin(-log10(padj), volcano_cap)]
  x[, contrast := name]
  x[, label := ""]
  # label the largest effects among confidently changed, named protein-coding genes
  named <- x[gene_type == "protein_coding" & padj < 1e-20]
  for (dir in c("up", "down")) {
    top <- head(named[direction == dir][order(-abs(log2FoldChange_shrunk))], 6)$gene_id
    x[gene_id %in% top, label := gene_name]
  }
  x
}))
volcano_data[, contrast := factor(contrast, levels = c("SAHA_vs_DMSO", "RMD_vs_DMSO"))]
count_genes <- function(name, dir) format(volcano_data[contrast == name & direction == dir, .N], big.mark = ",")
volcano_titles <- vapply(levels(volcano_data$contrast), function(name) {
  sprintf("%s\n%s genes higher, %s lower", contrast_label[[name]], count_genes(name, "up"), count_genes(name, "down"))
}, "")

volcano <- ggplot(volcano_data, aes(log2FoldChange_shrunk, y)) +
  geom_point(aes(colour = direction), size = 0.5, alpha = 0.5, stroke = 0) +
  geom_text_repel(aes(label = label), colour = col$ink, size = 2.5, max.overlaps = Inf,
                  segment.colour = col$muted, segment.size = 0.2, min.segment.length = 0,
                  box.padding = 0.5, force = 3, seed = 1) +
  facet_wrap(~contrast, labeller = as_labeller(volcano_titles)) +
  scale_y_sqrt(breaks = c(0, 5, 20, 50, 100, 150), expand = expansion(mult = c(0.01, 0.06))) +
  scale_colour_manual(
    values = c(up = col$up, down = col$down, ns = col$neutral),
    breaks = c("up", "down", "ns"),
    labels = c("higher with drug", "lower with drug", "not significant"),
    name = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "log2 fold change (shrunken)", y = "-log10 adjusted p-value") +
  theme_report()
save_png(volcano, snakemake@output[["volcano"]], 8, 4.6)

# ---------------------------------------------------------------- do the two drugs agree?

both <- merge(
  tables[["SAHA_vs_DMSO"]][, .(gene_id, gene_name, gene_type, saha = log2FoldChange, saha_padj = padj, saha_stat = stat)],
  tables[["RMD_vs_DMSO"]][, .(gene_id, rmd = log2FoldChange, rmd_padj = padj, rmd_stat = stat)],
  by = "gene_id"
)
both <- both[!is.na(saha) & !is.na(rmd)]
sig_saha <- both[!is.na(saha_padj) & saha_padj < alpha]
sig_rmd <- both[!is.na(rmd_padj) & rmd_padj < alpha]
sig_both <- both[!is.na(saha_padj) & !is.na(rmd_padj) & saha_padj < alpha & rmd_padj < alpha]

# Major-axis slope: treats the two fold changes symmetrically, unlike ordinary regression.
major_axis <- function(x, y) {
  v <- eigen(cov(cbind(x, y)))$vectors[, 1]
  v[2] / v[1]
}
concordance <- list(
  genes = nrow(both),
  pearson_r = cor(both$saha, both$rmd),
  spearman_rho = cor(both$saha, both$rmd, method = "spearman"),
  major_axis_slope = major_axis(both$saha, both$rmd),
  significant_saha = nrow(sig_saha),
  significant_rmd = nrow(sig_rmd),
  significant_both = nrow(sig_both),
  significant_both_same_direction = sig_both[sign(saha) == sign(rmd), .N],
  saha_genes_also_significant_with_rmd = nrow(sig_both) / nrow(sig_saha),
  median_abs_lfc_saha_in_shared = median(abs(sig_both$saha)),
  median_abs_lfc_rmd_in_shared = median(abs(sig_both$rmd))
)
write_json(concordance, snakemake@output[["concordance_stats"]], auto_unbox = TRUE, pretty = TRUE, digits = 4)

limit <- max(abs(quantile(c(both$saha, both$rmd), c(0.001, 0.999))))
scatter <- ggplot(both, aes(saha, rmd)) +
  geom_bin_2d(bins = 110) +
  geom_abline(slope = 1, intercept = 0, colour = col$muted, linewidth = 0.4) +
  geom_abline(slope = concordance$major_axis_slope,
              intercept = mean(both$rmd) - concordance$major_axis_slope * mean(both$saha),
              colour = col$up, linewidth = 0.6) +
  annotate("text", x = -limit, y = limit, hjust = 0, vjust = 1, colour = col$ink2, size = 3,
           label = sprintf("Pearson r = %.2f\nslope = %.2f (red line)\ngrey line: equal effect",
                           concordance$pearson_r, concordance$major_axis_slope)) +
  scale_fill_gradientn(colours = col$sequential, trans = "log10", name = "genes per bin") +
  coord_equal(xlim = c(-limit, limit), ylim = c(-limit, limit)) +
  labs(x = "log2 fold change, vorinostat vs DMSO", y = "log2 fold change, romidepsin vs DMSO") +
  theme_report() +
  theme(legend.position = "right", legend.justification = "center")
save_png(scatter, snakemake@output[["concordance"]], 5.6, 4.6)

# ---------------------------------------------------------------- heatmap of the most consistent genes

# Protein-coding genes changed in the same direction by both drugs, ranked by the weaker of
# the two Wald statistics: a gene ranks high only if both drugs move it clearly.
n_each <- 20
ranked <- sig_both[gene_type == "protein_coding" & sign(saha) == sign(rmd)]
ranked[, `:=`(score = pmin(abs(saha_stat), abs(rmd_stat)), mean_lfc = (saha + rmd) / 2)]
setorder(ranked, -score)
top <- rbind(head(ranked[mean_lfc > 0], n_each), head(ranked[mean_lfc < 0], n_each))

z <- within_donor[top$gene_id, ]
z <- z / apply(z, 1, sd)
rownames(z) <- make.unique(top$gene_name)
gene_order <- rownames(z)[order(top$mean_lfc)]
heat <- as.data.table(as.table(z))
setnames(heat, c("gene", "sample", "z"))
heat <- merge(heat, samples, by = "sample")
heat[, `:=`(gene = factor(gene, levels = gene_order), donor = factor(donor))]

heatmap <- ggplot(heat, aes(donor, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(~treatment, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = col$down, mid = col$mid, high = col$up, midpoint = 0,
                       name = "Expression relative to\nthe donor's mean (z-score)") +
  labs(x = NULL, y = NULL) +
  theme_report(base_size = 9) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(),
        strip.text = element_text(hjust = 0.5), axis.text.y = element_text(size = 6.5),
        legend.position = "right", legend.justification = "center")
save_png(heatmap, snakemake@output[["heatmap"]], 6, 6.8)
