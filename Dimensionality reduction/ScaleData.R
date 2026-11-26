#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

# 1. Parsing arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to normalized/HVG Seurat RDS file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "seurat_pca.rds",
              help = "Output RDS file path [default= %default]", metavar = "character"),
  make_option(c("-p", "--plot"), type = "character", default = "elbow_plot.png",
              help = "Path to save Elbow plot [default= %default]", metavar = "character")
)

parser <- OptionParser(usage = "Usage: %prog -i <featuer.object.rds> -o <seurat_pca.rds>", option_list = option_list)
args <- parse_args(parser)

if (is.null(args$input)) {
  print_help(parser)
  stop("Error: Input RDS file is required (-i/--input).", call. = FALSE)
}

# 2. Load Seurat Object
cat(paste0("--> Loading Seurat object: ", args$input, "\n"))
seurat_obj <- readRDS(args$input)

# 3. Scale ALL Genes in the dataset
cat("--> Scaling ALL genes in the matrix (this may take a few moments)...\n")
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes, verbose = FALSE)

# 4. Run PCA (Uses the HVGs set by modelGeneVar by default)
cat("--> Running PCA across 50 dimensions...\n")
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

# 5. Generate and save Elbow Plot
cat(paste0("--> Generating Elbow Plot: ", args$plot, "\n"))
p <- ElbowPlot(seurat_obj, ndims = 50) +
  ggtitle("Elbow Plot for PC Selection") +
  theme_minimal()

ggsave(filename = args$plot, plot = p, width = 8, height = 6, dpi = 300)

# 6. Save Object
cat(paste0("\n--> Saving scaled object with PCA to: ", args$output, "\n"))
saveRDS(seurat_obj, file = args$output)
cat("--> Process completed successfully!\n")
