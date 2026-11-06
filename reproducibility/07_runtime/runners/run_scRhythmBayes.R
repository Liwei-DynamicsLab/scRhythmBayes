# Isolated scRhythmBayes CellPhase runner; original source is unchanged.
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
})
input_dir <- "outputs/simulation/output/scRITMO_extended_simulation/sigma_gradient_1p5_2p5_3p5_4p5_5p5_6p5_empirical_fixed_v1/R_input"
output_dir <- "07_runtime/runs/scRhythmBayes"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
count <- as(Matrix::readMM(file.path(input_dir, "count_gene_by_cell.mtx")), "CsparseMatrix")
genes <- read.csv(file.path(input_dir, "genes.csv"), stringsAsFactors = FALSE)
metadata <- read.csv(file.path(input_dir, "metadata.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(count) == nrow(genes), ncol(count) == nrow(metadata), !anyDuplicated(metadata$cell))
rownames(count) <- genes$gene
colnames(count) <- metadata$cell
metadata <- metadata[match(colnames(count), metadata$cell), , drop = FALSE]
source("../R/infer_cell_phase.R")
started <- Sys.time()
result <- infer_cell_phase(
  count = count,
  metadata = metadata,
  ZT_by = "ZT",
  strata_by = "celltype",
  Species = "mouse",
  min_seed = 0.01,
  family = "NB",
  setseed = 1,
  Verbose = TRUE
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
stopifnot(nrow(result) == nrow(metadata), !anyDuplicated(result$cell), setequal(result$cell, metadata$cell))
saveRDS(result, file.path(output_dir, "scRhythmBayes_cell_phase.rds"))
write.csv(result, file.path(output_dir, "scRhythmBayes_cell_phase.csv"), row.names = FALSE)
write.csv(data.frame(method_only_seconds = elapsed), file.path(output_dir, "method_timing.csv"), row.names = FALSE)
