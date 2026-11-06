# Use the R environment prepared for this analysis.
stopifnot(utils::packageVersion("SeuratObject") >= "5.0.0")
root <- "outputs/heart"
objects_dir <- "data/scHeart_by_celltype/objects"
out <- file.path(root, "figures_complete_noif_20261023")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
phase <- readRDS(file.path(root, "heart_infer_cell_phase.rds"))
rho <- readRDS(file.path(root, "heart_run_phase_rho.rds"))
seed_env <- new.env(parent = emptyenv())
load("data/merge_27seedgenes_bmal1.Rdata", envir = seed_env)
seeds <- unique(toupper(as.character(unlist(as.list(seed_env), use.names = FALSE))))
seeds[seeds == "BMAL1"] <- "ARNTL"
celltypes <- c("Cardiomyocytes", "Endothelial", "Lymphocytes", "Myeloid", "Adipocytes",
               "Neuron", "Pericytes", "SMC", "Fibroblasts", "Myofibroblasts")
meta_list <- expression_list <- source_list <- vector("list", length(celltypes))
names(meta_list) <- names(expression_list) <- names(source_list) <- celltypes

for (celltype in celltypes) {
  message("Read original subset: ", celltype)
  source_file <- file.path(objects_dir, paste0(celltype, "_object.Rdata"))
  e <- new.env(parent = emptyenv())
  load(source_file, envir = e)
  obj <- e$obj_use
  stopifnot(inherits(obj, "Seurat"), "RNA" %in% SeuratObject::Assays(obj))
  z <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")
  stopifnot(nrow(z) == 24386L)
  genes <- intersect(seeds, rownames(z))
  stopifnot(length(genes) >= 20L)
  m <- obj@meta.data
  m$cell <- rownames(m)
  stopifnot(identical(colnames(z), m$cell), all(as.character(m$celltype) == celltype))
  ix <- match(m$cell, rho$data$cell)
  stopifnot(!anyNA(ix))
  stopifnot(all(as.character(m$Group) == as.character(rho$data$Group[ix])))
  stopifnot(all(as.character(m$Sample) == as.character(rho$data$Sample[ix])))
  expression_list[[celltype]] <- z[genes, , drop = FALSE]
  meta_list[[celltype]] <- m
  source_list[[celltype]] <- data.frame(
    celltype = celltype, source = source_file, n_cell = ncol(z), n_gene_RNA = nrow(z),
    n_clock_gene = length(genes), source_bytes = file.info(source_file)$size
  )
  rm(obj, z, m, e)
  invisible(gc())
}

heart_meta <- do.call(rbind, unname(meta_list))
rownames(heart_meta) <- heart_meta$cell
stopifnot(all(vapply(expression_list, function(x) identical(rownames(x), rownames(expression_list[[1]])), logical(1))))
heart_data <- do.call(cbind, expression_list)
stopifnot(identical(heart_meta$cell, colnames(heart_data)))
stopifnot(setequal(heart_meta$cell, phase$cell), ncol(heart_data) == 127942L)
stopifnot(!anyDuplicated(heart_meta$cell))
stopifnot(all(is.finite(heart_data@x)), all(heart_data@x >= 0))
heart_meta$ZT_num <- as.numeric(sub("^ZT", "", as.character(heart_meta$ZT)))
heart_meta$Group <- factor(as.character(heart_meta$Group), levels = c("NCR", "ACR"))
heart_meta$celltype <- factor(as.character(heart_meta$celltype), levels = celltypes)
design <- readr::read_tsv(file.path(root, "verified_sample_design.tsv"), show_col_types = FALSE)
stopifnot(all(heart_meta$ZT_num == design$ZT_num[match(heart_meta$Sample, design$Sample)]))

saveRDS(list(
  metadata = heart_meta, clock_expression = heart_data, celltype_levels = celltypes,
  RNA_features = 24386L, phase_md5 = unname(tools::md5sum(file.path(root, "heart_infer_cell_phase.rds"))),
  source_manifest = do.call(rbind, source_list),
  provenance = "Read-only original per-celltype subsets; all cells retained; RNA data layer; no normalization/refitting"
), file.path(out, "heart_plot_inputs_20261023.rds"))
readr::write_tsv(do.call(rbind, source_list), file.path(out, "plot_input_provenance.tsv"))
cat("Prepared", ncol(heart_data), "cells and", nrow(heart_data), "clock genes; no model fit.\n")
