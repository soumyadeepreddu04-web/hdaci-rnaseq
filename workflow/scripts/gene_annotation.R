# Gene-level annotation table from the combined GENCODE + ERCC GTF.
suppressPackageStartupMessages(library(data.table))

gtf <- snakemake@input[["gtf"]]

attribute <- function(x, key) {
  hit <- regmatches(x, regexpr(sprintf('%s "[^"]+"', key), x))
  sub(sprintf('%s "([^"]+)"', key), "\\1", hit)
}

genes <- fread(
  cmd = sprintf("awk -F'\\t' '$3 == \"gene\"' %s", shQuote(gtf)),
  sep = "\t", header = FALSE, quote = "", select = c(1, 9), col.names = c("chr", "attributes")
)
genes <- genes[, .(
  gene_id = attribute(attributes, "gene_id"),
  gene_name = attribute(attributes, "gene_name"),
  gene_type = attribute(attributes, "gene_type"),
  chr
)]

# The ERCC annotation has exon lines only, one per spike-in.
ercc <- fread(
  cmd = sprintf("awk -F'\\t' '$2 == \"ERCC\"' %s", shQuote(gtf)),
  sep = "\t", header = FALSE, quote = "", select = 1, col.names = "chr"
)
ercc <- ercc[, .(gene_id = chr, gene_name = chr, gene_type = "ERCC_spike_in", chr)]

annotation <- rbind(genes, ercc)
# Ensembl id without the version suffix, for joining to other annotation releases
annotation[, ensembl_id := sub("\\.[0-9]+(_PAR_Y)?$", "", gene_id)]
stopifnot(!anyDuplicated(annotation$gene_id))

fwrite(annotation, snakemake@output[[1]], sep = "\t")
