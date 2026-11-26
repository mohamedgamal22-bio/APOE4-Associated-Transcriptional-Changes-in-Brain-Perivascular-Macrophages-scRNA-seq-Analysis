#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

# 1. إعداد الـ Arguments مع تحسين نمط النمط الافتراضي ليشمل الإنسان والفئران
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to input Seurat RDS file", metavar="character"),
  make_option(c("-m", "--mads"), type="numeric", default=3,
              help="Number of MADs for cutoff calculation [default= %default]", metavar="numeric"),
  make_option(c("-p", "--prefix"), type="character", default="(?i)^mt[-_]",
              help="Regex pattern for Mitochondrial genes [default= %default]", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default="mito_cutoff.txt",
              help="Text file to save suggested cutoff [default= %default]", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) {
  stop("Error: Input RDS file is required (-i/--input).", call.=FALSE)
}

# 2. تحميل Seurat Object
seurat_obj <- readRDS(opt$input)

# 3. حساب نسبة الميتوكندريا في Seurat
# ملاحظة: تم إدراج الـ pattern المعدل
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = opt$prefix)

mito_percentages <- seurat_obj$percent.mt

# التأكد من العثور على جينات ميتوكندريا
if (all(mito_percentages == 0, na.rm = TRUE)) {
  warning("Warning: Percent MT is 0 for all cells! Check gene names or Ensembl ID formats.")
}

# 4. حساب الـ MAD والـ Cutoff المقترح
data_median <- median(mito_percentages, na.rm = TRUE)
data_mad <- mad(mito_percentages, na.rm = TRUE)

calculated_cutoff <- data_median + (opt$mads * data_mad)
mito_outliers <- mito_percentages > calculated_cutoff

# 5. طباعة التقرير الشامل
cat("\n==========================================\n")
cat("  Mitochondrial QC Threshold (Seurat Engine) \n")
cat("==========================================\n")
cat(sprintf("Total Cells Evaluated  : %d\n", length(mito_percentages)))
cat(sprintf("Median %% Mitochondrial : %.2f%%\n", data_median))
cat(sprintf("MAD Value              : %.2f%%\n", data_mad))
cat(sprintf("Applied Multiplier     : %.2f MADs\n", opt$mads))
cat("------------------------------------------\n")
cat(sprintf("--> RECOMMENDED CUTOFF : %.2f%%\n", calculated_cutoff))
cat(sprintf("--> EXPECTED REMOVAL   : %d cells (%.2f%%)\n",
            sum(mito_outliers, na.rm = TRUE),
            (sum(mito_outliers, na.rm = TRUE) / length(mito_percentages)) * 100))
cat("==========================================\n\n")

# 6. حفظ الرقم في ملف نصي
writeLines(as.character(round(calculated_cutoff, 2)), opt$out_file)
cat(sprintf("[INFO] Suggested cutoff saved to: %s\n", opt$out_file))
