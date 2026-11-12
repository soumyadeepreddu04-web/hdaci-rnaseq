# Count one library under each strandedness setting. For a reverse-stranded (dUTP) library
# almost all assignable pairs are assigned with strandSpecific = 2 and almost none with 1.
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages(library(Rsubread))

modes <- c(unstranded = 0, forward = 1, reverse = 2)

rows <- lapply(names(modes), function(mode) {
  fc <- featureCounts(
    files = snakemake@input[["bam"]],
    annot.ext = snakemake@input[["gtf"]],
    isGTFAnnotationFile = TRUE,
    isPairedEnd = TRUE,
    countMultiMappingReads = FALSE,
    strandSpecific = modes[[mode]],
    nthreads = snakemake@threads,
    tmpDir = tempdir()
  )
  stat <- setNames(fc$stat[[2]], fc$stat$Status)
  data.frame(
    setting = mode,
    strandSpecific = modes[[mode]],
    assigned = stat[["Assigned"]],
    no_feature = stat[["Unassigned_NoFeatures"]],
    ambiguous = stat[["Unassigned_Ambiguity"]],
    assigned_fraction = round(stat[["Assigned"]] / sum(stat), 4)
  )
})

write.table(do.call(rbind, rows), snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE)
