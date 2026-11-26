#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

# 1. Parsing Arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to PCA Seurat RDS file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "seurat_harmony.rds",
              help = "Output RDS file path [default= %default]", metavar = "character"),
  make_option(c("-d", "--dims"), type = "integer", default = 20,
              help = "Number of PCs for Harmony and UMAP [default= %default]", metavar = "number")
)

parser <- OptionParser(usage = "Usage: %prog -i <pca_object.rds> -o <seurat_harmony.rds>", option_list = option_list)
args <- parse_args(parser)

if (is.null(args$input)) {
  print_help(parser)
  stop("Error: Input RDS file is required (-i/--input).", call. = FALSE)
}

# 2. Load Object
cat(paste0("--> Loading Seurat object: ", args$input, "\n"))
seurat_obj <- readRDS(args$input)

# 3. Run Harmony Batch Correction
cat(paste0("--> Running Harmony Integration across first ", args$dims, " PCs...\n"))
seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "orig.ident",
  dims.use = 1:args$dims,
  verbose = FALSE
)

# 4. Generate UMAP Embeddings (Uncorrected PCA vs Corrected Harmony)
cat("--> Running UMAP on original PCA embeddings...\n")
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:args$dims, reduction.name = "umap_uncorrected", verbose = FALSE)

cat("--> Running UMAP on corrected Harmony embeddings...\n")
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:args$dims, reduction.name = "umap_harmony", verbose = FALSE)

# 5. Plot Comparison
cat("--> Generating comparison plots (Before vs After Harmony)...\n")
p1 <- DimPlot(seurat_obj, reduction = "umap_uncorrected", group.by = "orig.ident") +
  ggtitle("Before Integration (PCA)") + theme_minimal()

p2 <- DimPlot(seurat_obj, reduction = "umap_harmony", group.by = "orig.ident") +
  ggtitle("After Integration (Harmony)") + theme_minimal()

p_combined <- p1 + p2
ggsave(filename = "harmony_comparison.png", plot = p_combined, width = 12, height = 5, dpi = 300)
cat("--> Comparison saved to 'harmony_comparison.png'\n")

# 6. Save Object
cat(paste0("\n--> Saving object with Harmony embeddings to: ", args$output, "\n"))
saveRDS(seurat_obj, file = args$output)
cat("--> Step completed successfully!\n")
