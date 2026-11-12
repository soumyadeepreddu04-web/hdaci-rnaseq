# Gene-level read-pair counts for one library with featureCounts.
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages(library(Rsubread))

bam <- snakemake@input[["bam"]]

# A pair is counted when it overlaps the exons of exactly one gene on the expected strand.
# Pairs overlapping more than one gene are not counted (multi-mapping pairs were already
# left out at the alignment step).
fc <- featureCounts(
  files = bam,
  annot.ext = snakemake@input[["gtf"]],
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  isPairedEnd = TRUE,
  countReadPairs = TRUE,
  countMultiMappingReads = FALSE,
  allowMultiOverlap = FALSE,
  strandSpecific = snakemake@params[["strand"]],
  nthreads = snakemake@threads,
  tmpDir = tempdir()
)

write.table(
  data.frame(gene_id = fc$annotation$GeneID, length = fc$annotation$Length, count = fc$counts[, 1]),
  snakemake@output[["counts"]],
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Same layout as the command-line featureCounts summary, so MultiQC can read it
summary <- fc$stat
colnames(summary) <- c("Status", bam)
write.table(summary, snakemake@output[["summary"]], sep = "\t", quote = FALSE, row.names = FALSE)
