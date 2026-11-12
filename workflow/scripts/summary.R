# Collect the headline numbers of a run into one JSON file, plus the software versions.
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

qc <- fread(snakemake@input[["qc"]])
gsea <- fread(snakemake@input[["gsea"]])
donors <- fread(snakemake@input[["donors"]])

# A library is flagged when its mapping or gene-assignment rate is more than five percentage
# points below the median of all libraries.
flag <- function(x) qc$sample[x < median(x) - 0.05]
flagged <- union(flag(qc$uniquely_mapped_fraction), flag(qc$assigned_fraction))

range_of <- function(x) list(min = min(x), median = median(x), max = max(x))
top_sets <- function(name, direction) {
  x <- gsea[contrast == name & padj < 0.05 & sign(NES) == direction][order(-abs(NES))]
  sub("^HALLMARK_", "", head(x$pathway, 5))
}

drug_sets <- dcast(gsea[contrast %in% c("SAHA_vs_DMSO", "RMD_vs_DMSO")], pathway ~ contrast, value.var = c("NES", "padj"))
sets_both <- drug_sets[padj_SAHA_vs_DMSO < 0.05 & padj_RMD_vs_DMSO < 0.05]

packages <- c("Rsubread", "Rfastp", "DESeq2", "ashr", "fgsea", "ggplot2", "data.table")
tool_version <- function(command) {
  tryCatch(trimws(system2(command, "--version", stdout = TRUE, stderr = TRUE)[1]), error = function(e) NA)
}

summary <- list(
  dataset = list(geo = "GSE114883", bioproject = "PRJNA472973", libraries = nrow(qc)),
  max_read_pairs = snakemake@params[["max_read_pairs"]],
  gencode_release = snakemake@params[["gencode"]],
  qc = list(
    read_pairs = range_of(qc$read_pairs),
    q30_rate = range_of(qc$q30_rate),
    pairs_after_trimming_fraction = range_of(qc$pairs_after_trimming / qc$read_pairs),
    uniquely_mapped_fraction = range_of(qc$uniquely_mapped_fraction),
    assigned_fraction = range_of(qc$assigned_fraction),
    assigned_pairs = range_of(qc$assigned_pairs),
    ercc_fraction_of_assigned = range_of(qc$ercc_fraction_of_assigned),
    genes_detected = range_of(qc$genes_detected),
    flagged_libraries = as.list(flagged)
  ),
  strandedness = fread(snakemake@input[["strandedness"]]),
  differential_expression = fromJSON(snakemake@input[["de"]]),
  leave_one_donor_out = donors,
  drug_concordance = fromJSON(snakemake@input[["concordance"]]),
  gsea = list(
    sets_tested = nrow(drug_sets),
    significant_saha = drug_sets[padj_SAHA_vs_DMSO < 0.05, .N],
    significant_rmd = drug_sets[padj_RMD_vs_DMSO < 0.05, .N],
    significant_both = nrow(sets_both),
    significant_both_same_direction = sets_both[sign(NES_SAHA_vs_DMSO) == sign(NES_RMD_vs_DMSO), .N],
    opposite_direction = as.list(sub("^HALLMARK_", "", sets_both[sign(NES_SAHA_vs_DMSO) != sign(NES_RMD_vs_DMSO), pathway])),
    SAHA_vs_DMSO = list(up = top_sets("SAHA_vs_DMSO", 1), down = top_sets("SAHA_vs_DMSO", -1)),
    RMD_vs_DMSO = list(up = top_sets("RMD_vs_DMSO", 1), down = top_sets("RMD_vs_DMSO", -1)),
    RMD_vs_SAHA = list(up = top_sets("RMD_vs_SAHA", 1), down = top_sets("RMD_vs_SAHA", -1))
  ),
  ercc = fromJSON(snakemake@input[["ercc"]]),
  published_comparison = fromJSON(snakemake@input[["published"]]),
  software = c(
    list(R = R.version.string, Bioconductor = as.character(BiocManager::version())),
    setNames(lapply(packages, function(p) as.character(packageVersion(p))), packages),
    list(snakemake = tool_version("snakemake"), multiqc = tool_version("multiqc"))
  )
)
write_json(summary, snakemake@output[["summary"]], auto_unbox = TRUE, pretty = TRUE, digits = 4)

for (p in packages) suppressPackageStartupMessages(library(p, character.only = TRUE))
writeLines(capture.output(sessionInfo()), snakemake@output[["session"]])
