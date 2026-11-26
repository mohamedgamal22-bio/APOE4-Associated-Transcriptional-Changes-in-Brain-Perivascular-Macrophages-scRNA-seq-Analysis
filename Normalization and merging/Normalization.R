#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(scran)
  library(SingleCellExperiment)
})

# 1. Parsing command-line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to Seurat RDS file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "seurat_norm_hvg.rds",
              help = "Output RDS file path [default= %default]", metavar = "character"),
  make_option(c("-v", "--var_output"), type = "character", default = NULL,
              help = "Optional output path to save gene_var object [default= %default]", metavar = "character"),
  make_option(c("-n", "--nhvgs"), type = "integer", default = 2000,
              help = "Number of HVGs to select [default= %default]", metavar = "number")
)

parser <- OptionParser(usage = "Usage: %prog -i <input.rds> -o <output.rds>", option_list = option_list)
args <- parse_args(parser)

if (is.null(args$input)) {
  print_help(parser)
  stop("Error: Input RDS file is required (-i/--input).", call. = FALSE)
}

# 2. Load Seurat Object
cat(paste0("--> Loading Seurat object: ", args$input, "\n"))
seurat_obj <- readRDS(args$input)

# 3. LogNormalize
cat("--> Performing Standard LogNormalization...\n")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

# 4. Join Layers (Required for Seurat v5 assay compatibility)
cat("--> Joining layers for Seurat v5 assay...\n")
seurat_obj <- JoinLayers(seurat_obj)

# 5. Feature Selection using scran::modelGeneVar
cat("--> Extracting normalized matrix for scran::modelGeneVar...\n")
norm_counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

# Check number of unique samples in orig.ident
num_samples <- length(unique(seurat_obj$orig.ident))

if (num_samples > 1) {
  cat(paste0("--> Detected ", num_samples, " samples in orig.ident. Building Design Matrix...\n"))
  design_mat <- model.matrix(~ orig.ident, data = seurat_obj@meta.data)
  cat("--> Modeling gene variance with design matrix...\n")
  gene_var <- modelGeneVar(norm_counts, design = design_mat)
} else {
  cat("--> Single sample detected in orig.ident. Modeling gene variance without design matrix...\n")
  gene_var <- modelGeneVar(norm_counts)
}

# Save gene_var object if requested
if (!is.null(args$var_output)) {
  cat(paste0("--> Saving gene_var statistics to: ", args$var_output, "\n"))
  saveRDS(gene_var, file = args$var_output)
}

# Extract top HVGs
cat(paste0("--> Selecting top ", args$nhvgs, " Highly Variable Genes (HVGs)...\n"))
top_hvgs <- getTopHVGs(gene_var, n = args$nhvgs)

# 6. Assign HVGs back to Seurat object
VariableFeatures(seurat_obj) <- top_hvgs

cat(paste0("--> Successfully set ", length(top_hvgs), " HVGs in Seurat object.\n"))

# 7. Save Object
cat(paste0("\n--> Saving processed Seurat object to: ", args$output, "\n"))
saveRDS(seurat_obj, file = args$output)
cat("--> Normalization and Feature Selection completed successfully!\n")
