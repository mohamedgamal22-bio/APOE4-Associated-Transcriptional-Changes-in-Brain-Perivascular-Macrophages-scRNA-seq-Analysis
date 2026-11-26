#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

# 1. إعداد الـ Arguments الشاملة والجديدة
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to input Seurat RDS file [Required]", metavar="character"),
  make_option(c("-m", "--mito_cutoff"), type="numeric", default=15,
              help="Mitochondrial percentage cutoff [default= %default]", metavar="numeric"),
  make_option(c("--min_features"), type="numeric", default=200,
              help="Minimum number of features/genes per cell [default= %default]", metavar="numeric"),
  make_option(c("--max_features"), type="numeric", default=6000,
              help="Maximum number of features/genes per cell [default= %default]", metavar="numeric"),
  make_option(c("--min_counts"), type="numeric", default=500,
              help="Minimum UMI count per cell [default= %default]", metavar="numeric"),
  make_option(c("--max_counts"), type="numeric", default=50000,
              help="Maximum UMI count per cell [default= %default]", metavar="numeric"),
  make_option(c("-p", "--prefix"), type="character", default="(?i)^mt[-_]",
              help="Regex pattern for Mitochondrial genes [default= %default]", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output RDS file path [default= auto-generated]", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) {
  stop("Error: Input RDS file is required (-i/--input).", call.=FALSE)
}

rds_path <- opt$input
mito_cutoff <- opt$mito_cutoff

# تحديد مسار ملف المخرجات تلقائياً إذا لم يُحدد
if (!is.null(opt$output)) {
  output_path <- opt$output
} else {
  base_name <- sub("\\.rds$", "", rds_path, ignore.case = TRUE)
  output_path <- paste0(base_name, "_filtered.rds")
}

cat(paste0("--> Loading RDS file: ", rds_path, "\n"))
cat(paste0("--> Applied Thresholds:\n"))
cat(sprintf("    - Mitochondrial Cutoff : <= %.2f%%\n", mito_cutoff))
cat(sprintf("    - nFeature_RNA Range   : %d - %d\n", opt$min_features, opt$max_features))
cat(sprintf("    - nCount_RNA Range     : %d - %d\n\n", opt$min_counts, opt$max_counts))

# 2. Read Seurat object
sample <- readRDS(rds_path)
original_cells <- ncol(sample)

# 3. Calculate percent.mt (Case-insensitive pattern for Human & Mouse)
if (!"percent.mt" %in% colnames(sample@meta.data)) {
  cat("--> 'percent.mt' not found in metadata. Calculating now...\n")
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = opt$prefix)
}

# Ensure basic metadata metrics exist
if (!"nFeature_RNA" %in% colnames(sample@meta.data)) {
  cat("--> 'nFeature_RNA' not found in metadata. Calculating now...\n")
  sample[["nFeature_RNA"]] <- Matrix::colSums(GetAssayData(sample, assay = "RNA", slot = "counts") > 0)
}

if (!"nCount_RNA" %in% colnames(sample@meta.data)) {
  cat("--> 'nCount_RNA' not found in metadata. Calculating now...\n")
  sample[["nCount_RNA"]] <- Matrix::colSums(GetAssayData(sample, assay = "RNA", slot = "counts"))
}

# 4. Apply MT cutoff filtering
sample_filtered <- subset(sample, subset = percent.mt <= mito_cutoff)
cat(paste0("--> Cells after MT filtering: ", ncol(sample_filtered), "\n"))

# 5. Apply dynamic nFeature_RNA and nCount_RNA filtering
sample_filtered <- subset(sample_filtered,
                          subset = nFeature_RNA > opt$min_features &
                                   nFeature_RNA < opt$max_features &
                                   nCount_RNA > opt$min_counts &
                                   nCount_RNA < opt$max_counts)

cat(paste0("--> Cells after nFeature & nCount filtering: ", ncol(sample_filtered), "\n"))

# 6. Print final summary
filtered_cells <- ncol(sample_filtered)
removed_cells <- original_cells - filtered_cells
pct_removed <- round((removed_cells / original_cells) * 100, 2)

cat("\n==========================================\n")
cat("          QC Filtering Summary            \n")
cat("==========================================\n")
cat(sprintf("Original cells count : %d\n", original_cells))
cat(sprintf("Filtered cells count : %d\n", filtered_cells))
cat(sprintf("Total cells removed  : %d (%.2f%%)\n", removed_cells, pct_removed))
cat("==========================================\n\n")

# 7. Save the filtered Seurat object
saveRDS(sample_filtered, file = output_path)
cat(paste0("--> Filtered object successfully saved to: ", output_path, "\n\n"))
