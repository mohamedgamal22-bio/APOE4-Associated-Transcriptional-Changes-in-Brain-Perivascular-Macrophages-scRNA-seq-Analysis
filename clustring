library(Seurat)
library(ggplot2)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input RDS file"),
  make_option(c("-o", "--output"), type="character", default="cluster_object.rds", help="Output RDS file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

seurat_obj <- readRDS(opt$input)

# 1. Neighbors & Clustering
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:17)
res_range <- seq(0.3, 1.1, by = 0.1)
seurat_obj <- FindClusters(seurat_obj, resolution = res_range)

# 2. UMAP calculation
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:17)

# 3. رسم وحفظ الـ UMAP لكل resolution
for (res in res_range) {
  res_col <- paste0("RNA_snn_res.", res)
  Idents(seurat_obj) <- res_col

  p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
       ggtitle(paste0("UMAP - Resolution ", res))

  ggsave(filename = paste0("umap_res_", res, ".png"), plot = p, width = 8, height = 6)
}

# 4. حفظ الـ Object النهائي
saveRDS(seurat_obj, opt$output)
cat("--> Clustering, UMAP calculation, and plots saved successfully!\n")
