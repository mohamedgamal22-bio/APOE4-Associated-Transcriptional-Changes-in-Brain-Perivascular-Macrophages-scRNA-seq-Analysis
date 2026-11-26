#!/usr/bin/env Rscript

# 1. Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Error: Missing arguments.\nUsage: Rscript QC_validation.R <seurat_object.rds> <mt_cutoff>", call. = FALSE)
}

rds_path <- args[1]
cutoff <- as.numeric(args[2])

cat(paste0("--> Loading RDS file: ", rds_path, "\n"))
cat(paste0("--> Mitochondrial Cutoff: ", cutoff, "%\n\n"))

# 2. Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(scuttle)
  library(edgeR)
})

# 3. Read Seurat object
sample <- readRDS(rds_path)

# 4. Calculate percent.mt if not present (Case-insensitive pattern for Human & Mouse)
mt_pattern <- "(?i)^mt[-_]"

if (!"percent.mt" %in% colnames(sample@meta.data)) {
  cat("--> Calculating 'percent.mt' for the Seurat object...\n")
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = mt_pattern)
}

# 5. Define discarded and kept cells
discard <- sample$percent.mt > cutoff
cat("Cell filtering summary (FALSE = Kept, TRUE = Discarded):\n")
print(table(discard))
cat("\n")

if (sum(discard, na.rm = TRUE) == 0) {
  stop("Warning: Zero cells were discarded with this cutoff! No comparison plot can be generated.", call. = FALSE)
}

# 6. Extract raw counts matrix (Compatible with Seurat v4 and v5)
if (packageVersion("Seurat") >= "5.0.0") {
  counts_mat <- GetAssayData(sample, assay = "RNA", layer = "counts")
} else {
  counts_mat <- GetAssayData(sample, assay = "RNA", slot = "counts")
}

# 7. Calculate average gene expression for lost and kept groups
lost_avg <- calculateAverage(counts_mat[, discard, drop = FALSE])
kept_avg <- calculateAverage(counts_mat[, !discard, drop = FALSE])

# 8. Compute Log-Fold Change via edgeR
logged <- cpm(cbind(lost_avg, kept_avg), log = TRUE, prior.count = 2)
logFC <- logged[, 1] - logged[, 2]
abundance <- rowMeans(logged)

# 9. Identify Mitochondrial Genes (Case-insensitive matching)
is_mito <- grepl(mt_pattern, names(logFC), perl = TRUE)

# 10. Generate and save the MA plot to PNG
output_plot <- "QC_mitochondrial_plot.png"
png(filename = output_plot, width = 800, height = 600, res = 120)

plot(
  abundance,
  logFC,
  xlab = "Average count",
  ylab = "Log-FC (lost/kept)",
  pch = 16,
  col = "black",
  cex = 0.6
)
abline(h = 0, lty = 2, col = "gray")

# Highlight mitochondrial genes in dodgerblue
points(
  abundance[is_mito],
  logFC[is_mito],
  col = "dodgerblue",
  pch = 16,
  cex = 0.8
)

dev.off()
cat(paste0("--> Plot successfully saved to: ", output_plot, "\n\n"))

# 11. Extract top N non-MT genes and calculate cell detection counts/percentages
non_mt_genes <- !is_mito
top_non_mt_genes <- sort(logFC[non_mt_genes], decreasing = TRUE)

top_n <- 15
genes_to_check <- head(top_non_mt_genes, top_n)

results <- data.frame(
  Gene = names(genes_to_check),
  LogFC = round(as.numeric(genes_to_check), 3),
  Cells_In_Lost = sapply(names(genes_to_check), function(g) sum(counts_mat[g, discard] > 0)),
  Cells_In_Kept = sapply(names(genes_to_check), function(g) sum(counts_mat[g, !discard] > 0)),
  Pct_In_Lost = sapply(names(genes_to_check), function(g) round((sum(counts_mat[g, discard] > 0) / sum(discard)) * 100, 2)),
  Pct_In_Kept = sapply(names(genes_to_check), function(g) round((sum(counts_mat[g, !discard] > 0) / sum(!discard)) * 100, 2))
)

# 12. Output results matrix
cat("Top Non-Mitochondrial Genes Enriched in Discarded Cells:\n")
print(results, row.names = FALSE)
