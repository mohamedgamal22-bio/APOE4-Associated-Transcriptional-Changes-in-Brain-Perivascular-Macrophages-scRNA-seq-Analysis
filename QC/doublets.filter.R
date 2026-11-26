#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
})

# 1. إعداد الـ Arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to input Seurat RDS file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="seurat_filtered_doublets.rds",
              help="Path to output RDS file after removing doublets", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Error: Input RDS file (-i) is required.", call.value=FALSE)
}

cat("--> Loading Seurat object from:", opt$input, "\n")
seurat_obj <- readRDS(opt$input)

cat("--> Converting Seurat object to SingleCellExperiment for scDblFinder...\n")
# تحويل كائن سيورات ليتوافق مع حزمة التحليل
sce <- as.SingleCellExperiment(seurat_obj)

cat("--> Running scDblFinder to detect doublets...\n")
# تشغيل خوارزمية الكشف عن الـ Doublets
sce <- scDblFinder(sce)

# إضافة نتائج الفحص إلى الـ Seurat metadata
seurat_obj$scDblFinder_class <- colData(sce)$scDblFinder.class
seurat_obj$scDblFinder_score <- colData(sce)$scDblFinder.score

# طباعة ملخص النتائج في التيرمنال
cat("\nDoublet Filtering Summary:\n")
print(table(seurat_obj$scDblFinder_class))

# فصل الخلايا السليمة (Singlets) عن المزدوجة (Doublets)
singlets_obj <- subset(seurat_obj, scDblFinder_class == "singlet")

cat(sprintf("\n--> Removed %d doublets. Kept %d singlets.\n",
            sum(seurat_obj$scDblFinder_class == "doublet"),
            ncol(singlets_obj)))

# حفظ الكائن الناتج بعد التصفية
saveRDS(singlets_obj, file = opt$output)
cat(sprintf("--> Filtered Seurat object successfully saved to: %s\n", opt$output))
