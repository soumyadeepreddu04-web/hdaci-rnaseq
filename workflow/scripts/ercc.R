# ERCC spike-in checks.
#
# The authors added ERCC Mix 1 to DMSO lysates and Mix 2 to drug-treated lysates, in proportion
# to cell number. The 92 spike-ins fall into four subgroups with known Mix 1 : Mix 2 ratios
# (A 4:1, B 1:1, C 0.67:1, D 0.5:1). That gives two things to test:
#
#   1. Fold-change accuracy. With subgroup B (equal in both mixes) as the reference, subgroups
#      A, C and D should show log2 fold changes of -2, +0.58 and +1 (drug vs DMSO).
#   2. The normalisation assumption. DESeq2's default size factors assume that most genes do
#      not change. HDAC inhibitors act on chromatin genome-wide, so that is worth checking:
#      size factors from subgroup B measure RNA per cell instead, and the difference between
#      the two normalisations is a global shift applied to every host fold change.
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(jsonlite)
})
source(file.path(snakemake@scriptdir, "theme.R"))

cfg <- snakemake@params[["de"]]
alpha <- cfg[["alpha"]]
drugs <- c(SAHA_vs_DMSO = "SAHA", RMD_vs_DMSO = "RMD")

conc <- fread(snakemake@input[["conc"]])
setnames(conc, c("sort_id", "gene_id", "subgroup", "mix1", "mix2", "ratio", "log2_ratio"))
conc[, expected_lfc := -log2_ratio]  # drug (Mix 2) relative to DMSO (Mix 1): -2, 0, 0.58, 1

samples <- fread(snakemake@input[["samples"]])[, .(sample, donor, treatment)]
counts <- fread(snakemake@input[["counts"]])
ercc <- as.matrix(counts[grepl("^ERCC-", gene_id), samples$sample, with = FALSE])
rownames(ercc) <- counts[grepl("^ERCC-", gene_id), gene_id]
ercc <- ercc[rowSums(ercc >= cfg[["min_count"]]) >= cfg[["min_samples"]], ]

main <- readRDS(snakemake@input[["dds"]])
subgroup_b <- intersect(conc[subgroup == "B", gene_id], rownames(ercc))

# ---- refit the same model with spike-in size factors
dds <- DESeqDataSetFromMatrix(
  rbind(counts(main), ercc), as.data.frame(colData(main)[, c("donor", "treatment")]),
  design = ~ donor + treatment
)
dds <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% subgroup_b)
dds <- DESeq(dds)

observed <- list()
shift <- list()
for (name in names(drugs)) {
  contrast <- c("treatment", drugs[[name]], "DMSO")
  spike <- as.data.table(as.data.frame(results(dds, contrast = contrast, alpha = alpha)), keep.rownames = "gene_id")
  default <- as.data.table(as.data.frame(results(main, contrast = contrast, alpha = alpha)), keep.rownames = "gene_id")

  observed[[name]] <- merge(conc[, .(gene_id, subgroup, expected_lfc)],
                            spike[, .(gene_id, baseMean, observed_lfc = log2FoldChange, lfcSE)], by = "gene_id")
  observed[[name]][, contrast := name]

  host <- merge(default[, .(gene_id, lfc_default = log2FoldChange)],
                spike[, .(gene_id, lfc_spike = log2FoldChange, padj)], by = "gene_id")
  shift[[name]] <- list(
    up_spike_in_normalised = host[!is.na(padj) & padj < alpha & lfc_spike > 0, .N],
    down_spike_in_normalised = host[!is.na(padj) & padj < alpha & lfc_spike < 0, .N]
  )
}
observed <- rbindlist(observed)
fwrite(observed, snakemake@output[["table"]])

# ---- the global shift, donor by donor
# Changing size factors moves every log2 normalised count in a library by the same amount
# (delta), so within a donor the shift in a drug-vs-DMSO fold change is delta(drug) - delta(DMSO).
# Each donor gives one estimate per drug; with four donors all four are reported, and the
# median is used as the summary because a single odd library can dominate the mean.
median_of_ratios <- function(m) {
  m <- m[rowSums(m == 0) == 0, , drop = FALSE]
  log_m <- log(m)
  exp(apply(log_m - rowMeans(log_m), 2, median))
}
delta <- log2(sizeFactors(main)) - log2(median_of_ratios(ercc[subgroup_b, , drop = FALSE]))
delta <- setNames(delta[samples$sample], samples$sample)

shifts <- rbindlist(lapply(names(drugs), function(name) {
  rbindlist(lapply(unique(samples$donor), function(d) {
    data.table(
      contrast = name,
      donor = d,
      shift_log2 = delta[[paste(d, drugs[[name]], sep = "_")]] - delta[[paste(d, "DMSO", sep = "_")]]
    )
  }))
}))
fwrite(shifts, snakemake@output[["shifts"]])

shift_summary <- shifts[, .(
  median_shift_log2 = median(shift_log2),
  min_shift_log2 = min(shift_log2),
  max_shift_log2 = max(shift_log2),
  donors_negative = sum(shift_log2 < 0)
), by = contrast]
for (name in names(drugs)) {
  shift_summary[contrast == name, `:=`(
    significant_default = sum(results(main, contrast = c("treatment", drugs[[name]], "DMSO"), alpha = alpha)$padj < alpha, na.rm = TRUE),
    significant_spike_in_normalised = shift[[name]]$up_spike_in_normalised + shift[[name]]$down_spike_in_normalised
  )]
}

# ---- accuracy of the fold changes in subgroups A, C and D
test <- observed[subgroup != "B"]
fit <- lm(observed_lfc ~ expected_lfc, data = test)
by_subgroup <- observed[, .(spike_ins = .N, expected = expected_lfc[1], observed_median = median(observed_lfc)),
                        by = .(contrast, subgroup)][order(contrast, subgroup)]

write_json(list(
  spike_ins_detected = nrow(ercc),
  subgroup_b_used = length(subgroup_b),
  accuracy_slope = unname(coef(fit)[2]),
  accuracy_r = test[, cor(observed_lfc, expected_lfc)],
  by_subgroup = by_subgroup,
  global_shift = shift_summary,
  global_shift_by_donor = shifts
), snakemake@output[["summary"]], auto_unbox = TRUE, pretty = TRUE, digits = 4)

# ---- figure
observed[, contrast := factor(contrast, levels = names(drugs))]
medians <- observed[, .(observed_lfc = median(observed_lfc)), by = .(contrast, subgroup, expected_lfc)]
labels <- unique(observed[, .(subgroup, expected_lfc)])
labels[, text := fifelse(subgroup == "B", "B (reference)", subgroup)]

accuracy <- ggplot(observed, aes(expected_lfc, observed_lfc)) +
  geom_abline(slope = 1, intercept = 0, colour = col$muted, linewidth = 0.4) +
  geom_point(colour = col$series, alpha = 0.6, size = 1.8, stroke = 0,
             position = position_jitter(width = 0.05, height = 0, seed = 1)) +
  geom_segment(data = medians, aes(x = expected_lfc - 0.16, xend = expected_lfc + 0.16, yend = observed_lfc),
               colour = col$ink, linewidth = 0.7) +
  geom_text(data = labels, aes(x = expected_lfc, y = Inf, label = text), vjust = 1.5, colour = col$ink2, size = 2.8) +
  facet_wrap(~contrast, labeller = as_labeller(c(SAHA_vs_DMSO = "Vorinostat vs DMSO", RMD_vs_DMSO = "Romidepsin vs DMSO"))) +
  scale_x_continuous(breaks = sort(unique(observed$expected_lfc))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(title = "Spike-in fold changes", x = "Expected log2 fold change (Mix 2 / Mix 1)",
       y = "Observed log2 fold change (drug / DMSO)") +
  theme_report()

drug_names <- c(SAHA_vs_DMSO = "Vorinostat", RMD_vs_DMSO = "Romidepsin")
shifts[, x := match(contrast, names(drug_names))]
shift_medians <- shifts[, .(shift_log2 = median(shift_log2)), by = x]
by_donor <- ggplot(shifts, aes(x, shift_log2)) +
  geom_hline(yintercept = 0, colour = col$muted, linewidth = 0.4) +
  geom_segment(data = shift_medians, aes(x = x - 0.3, xend = x + 0.3, yend = shift_log2),
               colour = col$ink, linewidth = 0.7) +
  geom_point(colour = col$series, alpha = 0.8, size = 2.2, stroke = 0) +
  geom_text_repel(aes(label = donor), colour = col$ink2, size = 2.6, direction = "y", nudge_x = 0.2,
                  hjust = 0, min.segment.length = Inf, seed = 1) +
  scale_x_continuous(breaks = seq_along(drug_names), labels = unname(drug_names), limits = c(0.5, 2.6)) +
  labs(title = "Global shift, by donor", x = NULL, y = "Shift in host log2 fold changes\nunder spike-in normalisation") +
  theme_report() +
  theme(panel.grid.major.x = element_blank())

plot <- accuracy + by_donor + plot_layout(widths = c(2.7, 1))
save_png(plot, snakemake@output[["figure"]], 9, 3.9)
