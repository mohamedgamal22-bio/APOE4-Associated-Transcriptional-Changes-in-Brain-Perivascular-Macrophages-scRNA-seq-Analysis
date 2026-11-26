
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(DESeq2)
  library(Matrix)
  library(dplyr)
  library(tibble)
})

# ============================================================
# 1. CLI Arguments
# ============================================================

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Input Seurat RDS file"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "pseudobulk_deseq2",
    help = "Output directory"
  ),
  make_option(
    c("-p", "--padj_cutoff"),
    type = "numeric",
    default = 0.05,
    help = "Adjusted p-value cutoff"
  )
)

args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$input)) {
  stop("Please provide an input RDS file using -i")
}

if (args$padj_cutoff <= 0 || args$padj_cutoff >= 1) {
  stop("padj_cutoff must be between 0 and 1")
}

# ============================================================
# 2. Load Seurat Object
# ============================================================

cat("--> Loading Seurat object...\n")

obj <- readRDS(args$input)

DefaultAssay(obj) <- "RNA"

expected <- c(
  "APOE3_88",
  "APOE3_89",
  "APOE3_90",
  "APOE4_91",
  "APOE4_92",
  "APOE4_93"
)

samples <- as.character(obj$orig.ident)

if (
  anyNA(samples) ||
  !setequal(unique(samples), expected)
) {
  stop(
    "Unexpected sample IDs: ",
    paste(unique(samples), collapse = ", ")
  )
}

# ============================================================
# 3. Define Experimental Conditions
# ============================================================

obj$condition <- ifelse(
  grepl("^APOE3_", samples),
  "APOE3",
  "APOE4"
)

obj$condition <- factor(
  obj$condition,
  levels = c("APOE3", "APOE4")
)

cat("\n--> Sample distribution:\n")

print(table(
  obj$orig.ident,
  obj$condition
))

# ============================================================
# 4. Check Raw RNA Counts
# ============================================================

rna_layers <- Layers(obj[["RNA"]])

count_layers <- grep(
  "^counts(\\.|$)",
  rna_layers,
  value = TRUE
)

if (length(count_layers) == 0) {
  stop("No raw RNA count layers found.")
}

cat("\n--> Raw count layers:\n")
print(count_layers)

# Join multiple count layers if necessary
if (length(count_layers) > 1) {

  cat("--> Joining RNA count layers...\n")

  obj[["RNA"]] <- JoinLayers(
    obj[["RNA"]],
    layers = count_layers,
    new = "counts"
  )
}

# ============================================================
# 5. Pseudobulk Aggregation
# ============================================================

cat("\n--> Aggregating raw counts by sample...\n")

pb <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = "orig.ident",
  return.seurat = FALSE,
  verbose = FALSE
)$RNA

# Seurat replaces underscores with dashes.
# Restore the original sample names.

colnames(pb) <- gsub(
  "-",
  "_",
  colnames(pb),
  fixed = TRUE
)

cat("\n--> Aggregated sample names:\n")
print(colnames(pb))

# Validate sample names
if (!setequal(colnames(pb), expected)) {
  stop(
    "Unexpected pseudobulk columns: ",
    paste(colnames(pb), collapse = ", ")
  )
}

# Ensure consistent column ordering
pb <- pb[, expected, drop = FALSE]

cat("\n--> Pseudobulk dimensions:\n")
print(dim(pb))

cat("\n--> Total counts per sample:\n")
print(colSums(pb))

# ============================================================
# 6. Prepare DESeq2 Metadata
# ============================================================

col_data <- data.frame(
  sample = expected,
  condition = factor(
    c(
      rep("APOE3", 3),
      rep("APOE4", 3)
    ),
    levels = c("APOE3", "APOE4")
  ),
  row.names = expected
)

cat("\n--> DESeq2 sample metadata:\n")
print(col_data)

# ============================================================
# 7. Validate Raw Counts
# ============================================================

if (anyNA(pb) || any(pb < 0)) {
  stop("Counts contain missing or negative values.")
}

if (inherits(pb, "sparseMatrix")) {
  values <- pb@x
} else {
  values <- as.vector(pb)
}

if (any(abs(values - round(values)) > 1e-6)) {
  stop("Counts are not integers. Check raw RNA counts.")
}

# ============================================================
# 8. Create Output Directory
# ============================================================

dir.create(
  args$outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 9. Create DESeq2 Dataset
# ============================================================

cat("\n--> Creating DESeq2 dataset...\n")

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(pb),
  colData = col_data,
  design = ~ condition
)

# Filter low-expression genes
keep <- rowSums(counts(dds) >= 10) >= 3

dds <- dds[keep, ]

cat("--> Genes retained:", nrow(dds), "\n")

if (nrow(dds) == 0) {
  stop("No genes passed the expression filter.")
}

# ============================================================
# 10. Run DESeq2
# ============================================================

cat("\n--> Running DESeq2...\n")

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(
    "condition",
    "APOE4",
    "APOE3"
  ),
  alpha = args$padj_cutoff
)

# ============================================================
# 11. Prepare Results
# ============================================================

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

# Save all genes
write.csv(
  res_df,
  file.path(
    args$outdir,
    "APOE4_vs_APOE3_All_DEGs.csv"
  ),
  row.names = FALSE
)

# Save significant genes
sig_df <- res_df %>%
  filter(
    !is.na(padj),
    padj < args$padj_cutoff
  )

write.csv(
  sig_df,
  file.path(
    args$outdir,
    "APOE4_vs_APOE3_Significant_DEGs.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 12. Save Pseudobulk Counts
# ============================================================

write.csv(
  as.data.frame(pb) %>%
    rownames_to_column("gene"),
  file.path(
    args$outdir,
    "Pseudobulk_Counts.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 13. Extract Maf and Igf1
# ============================================================

genes_of_interest <- c(
  "Maf",
  "Igf1",
  "MAF",
  "IGF1"
)

target_results <- res_df %>%
  filter(gene %in% genes_of_interest)

write.csv(
  target_results,
  file.path(
    args$outdir,
    "Maf_Igf1_Results.csv"
  ),
  row.names = FALSE
)

cat("\n--> Maf / Igf1 results:\n")
print(target_results)

# ============================================================
# 14. Final Summary
# ============================================================

cat("\n========================================\n")
cat("       PSEUDOBULK ANALYSIS COMPLETE\n")
cat("========================================\n")

cat("Total genes tested:", nrow(res_df), "\n")

cat(
  "Significant DEGs (padj <",
  args$padj_cutoff,
  "):",
  nrow(sig_df),
  "\n"
)

cat(
  "Upregulated in APOE4:",
  sum(sig_df$log2FoldChange > 0, na.rm = TRUE),
  "\n"
)

cat(
  "Downregulated in APOE4:",
  sum(sig_df$log2FoldChange < 0, na.rm = TRUE),
  "\n"
)

cat("\nResults saved to:", args$outdir, "\n")
