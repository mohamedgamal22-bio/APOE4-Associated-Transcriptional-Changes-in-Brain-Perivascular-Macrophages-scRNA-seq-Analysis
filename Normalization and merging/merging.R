#!/usr/bin/env Rscript

# ============================================================
# CPU configuration - use up to 8 threads
# ============================================================

Sys.setenv(
  OMP_NUM_THREADS = "8",
  OPENBLAS_NUM_THREADS = "8",
  MKL_NUM_THREADS = "8",
  VECLIB_MAXIMUM_THREADS = "8",
  NUMEXPR_NUM_THREADS = "8"
)

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
})

cat("CPU cores available: ", parallel::detectCores(), "\n")
cat("Requested threads: 8\n")
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
})

# ============================================================
# Arguments
# ============================================================

option_list <- list(
  make_option(
    c("-f", "--files"),
    type = "character",
    default = NULL,
    help = "Comma-separated list of normalized Seurat RDS files"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "seurat_merged_hvg.rds",
    help = "Output RDS file"
  ),
  make_option(
    c("-n", "--nhvgs"),
    type = "integer",
    default = 2000,
    help = "Number of consensus HVGs"
  )
)

args <- parse_args(
  OptionParser(option_list = option_list)
)

if (is.null(args$files)) {
  stop("You must provide files using -f")
}

files <- trimws(
  strsplit(args$files, ",")[[1]]
)

if (length(files) < 2) {
  stop("At least 2 input files are required.")
}

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {
  stop(
    "Missing files: ",
    paste(missing_files, collapse = ", ")
  )
}

# ============================================================
# Memory helper
# ============================================================

print_memory <- function(stage) {

  gc(verbose = FALSE)

  cat("\n--- MEMORY: ", stage, " ---\n", sep = "")

  if (file.exists("/proc/meminfo")) {

    x <- readLines("/proc/meminfo")

    wanted <- grep(
      "^(MemAvailable|SwapTotal|SwapFree):",
      x,
      value = TRUE
    )

    cat(paste(wanted, collapse = "\n"), "\n")
  }

  cat("----------------------------------------\n")
}

# ============================================================
# Start
# ============================================================

cat("\n========================================\n")
cat("SEURAT MERGE\n")
cat("========================================\n")

cat("Samples: ", length(files), "\n", sep = "")
cat("HVG target: ", args$nhvgs, "\n\n", sep = "")

cat(
  "Seurat version: ",
  as.character(packageVersion("Seurat")),
  "\n",
  sep = ""
)

cat(
  "SeuratObject version: ",
  as.character(packageVersion("SeuratObject")),
  "\n\n",
  sep = ""
)

# ============================================================
# 1. Load objects
# ============================================================

objects <- vector(
  mode = "list",
  length = length(files)
)

sample_hvgs <- vector(
  mode = "list",
  length = length(files)
)

sample_names <- character(length(files))

for (i in seq_along(files)) {

  s_name <- gsub(
    "\\.norm\\.rds$",
    "",
    basename(files[i])
  )

  s_name <- gsub(
    "\\.rds$",
    "",
    s_name
  )

  sample_names[i] <- s_name

  cat(
    sprintf(
      "--> Loading [%d/%d]: %s\n",
      i,
      length(files),
      s_name
    )
  )

  objects[[i]] <- readRDS(files[i])

  # Make sample identity explicit
  objects[[i]]$orig.ident <- s_name

  # Get existing HVGs
  sample_hvgs[[i]] <- VariableFeatures(
    objects[[i]]
  )

  cat(
    "    Genes: ",
    nrow(objects[[i]]),
    " | Cells: ",
    ncol(objects[[i]]),
    " | HVGs: ",
    length(sample_hvgs[[i]]),
    "\n",
    sep = ""
  )

  if (length(sample_hvgs[[i]]) == 0) {
    stop(
      "Sample ",
      s_name,
      " has no VariableFeatures."
    )
  }

  gc(verbose = FALSE)
}

print_memory("All samples loaded")

# ============================================================
# 2. Consensus HVGs
# ============================================================

cat("\n--> Calculating consensus HVGs...\n")

hvg_frequency <- sort(
  table(unlist(sample_hvgs)),
  decreasing = TRUE
)

n_hvg <- min(
  args$nhvgs,
  length(hvg_frequency)
)

consensus_hvgs <- names(
  hvg_frequency
)[seq_len(n_hvg)]

cat(
  "--> Selected ",
  length(consensus_hvgs),
  " consensus HVGs.\n",
  sep = ""
)

rm(sample_hvgs, hvg_frequency)

gc(verbose = FALSE)

# ============================================================
# 3. ONE merge
# ============================================================

cat("\n========================================\n")
cat("MERGING ALL SAMPLES\n")
cat("========================================\n")

print_memory("Before merge")

merge_start <- Sys.time()

merged_obj <- merge(
  x = objects[[1]],
  y = objects[-1],
  add.cell.ids = sample_names,
  merge.data = TRUE
)

merge_end <- Sys.time()

cat(
  "\n--> Merge completed in ",
  round(
    as.numeric(
      difftime(
        merge_end,
        merge_start,
        units = "mins"
      )
    ),
    2
  ),
  " minutes.\n",
  sep = ""
)

cat(
  "--> Genes: ",
  nrow(merged_obj),
  "\n",
  sep = ""
)

cat(
  "--> Cells: ",
  ncol(merged_obj),
  "\n",
  sep = ""
)

# ============================================================
# IMPORTANT: release original 6 objects
# ============================================================

rm(objects)

gc(verbose = FALSE)

print_memory("After merge + removing original objects")

# ============================================================
# 4. Verify sample identities
# ============================================================

merged_obj$orig.ident <- factor(
  merged_obj$orig.ident,
  levels = sample_names
)

cat("\n--> Cells per sample:\n")
print(table(merged_obj$orig.ident))

# ============================================================
# 5. Join Seurat v5 layers
# ============================================================

cat("\n========================================\n")
cat("JOINING RNA LAYERS\n")
cat("========================================\n")

cat("\n--> Layers BEFORE JoinLayers:\n")
print(Layers(merged_obj, assay = "RNA"))

print_memory("Before JoinLayers")

join_start <- Sys.time()

# Use SeuratObject namespace explicitly.
# This avoids "could not find function JoinLayers".
merged_obj <- SeuratObject::JoinLayers(
  object = merged_obj,
  assay = "RNA"
)

join_end <- Sys.time()

gc(verbose = FALSE)

cat(
  "\n--> JoinLayers completed in ",
  round(
    as.numeric(
      difftime(
        join_end,
        join_start,
        units = "mins"
      )
    ),
    2
  ),
  " minutes.\n",
  sep = ""
)

cat("\n--> Layers AFTER JoinLayers:\n")
print(Layers(merged_obj, assay = "RNA"))

print_memory("After JoinLayers")

# ============================================================
# 6. Assign consensus HVGs
# ============================================================

consensus_hvgs <- intersect(
  consensus_hvgs,
  rownames(merged_obj)
)

VariableFeatures(merged_obj) <- consensus_hvgs

cat(
  "\n--> Assigned ",
  length(VariableFeatures(merged_obj)),
  " consensus HVGs.\n",
  sep = ""
)

# ============================================================
# 7. Final validation
# ============================================================

cat("\n========================================\n")
cat("FINAL OBJECT\n")
cat("========================================\n")

cat(
  "Genes:   ",
  nrow(merged_obj),
  "\n",
  sep = ""
)

cat(
  "Cells:   ",
  ncol(merged_obj),
  "\n",
  sep = ""
)

cat(
  "HVGs:    ",
  length(VariableFeatures(merged_obj)),
  "\n",
  sep = ""
)

cat(
  "Samples: ",
  length(unique(merged_obj$orig.ident)),
  "\n",
  sep = ""
)

cat("\nSample counts:\n")
print(table(merged_obj$orig.ident))

cat("\nFinal RNA layers:\n")
print(Layers(merged_obj, assay = "RNA"))

# ============================================================
# 8. Save
# ============================================================

print_memory("Before save")

cat(
  "\n--> Saving: ",
  args$output,
  "\n",
  sep = ""
)

save_start <- Sys.time()

saveRDS(
  merged_obj,
  file = args$output,
  compress = FALSE
)

save_end <- Sys.time()

cat(
  "--> Save completed in ",
  round(
    as.numeric(
      difftime(
        save_end,
        save_start,
        units = "mins"
      )
    ),
    2
  ),
  " minutes.\n",
  sep = ""
)

cat("\n========================================\n")
cat("DONE\n")
cat("========================================\n")
