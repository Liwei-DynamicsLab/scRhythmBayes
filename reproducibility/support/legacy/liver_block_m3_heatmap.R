options(future.globals.maxSize = Inf)

build_scRhythm_OVR_fullmatrix_heatmap <- function(
    expr_mat,
    cell_meta,
    ovr_df,
    out_dir = "data/scLiver_global_phase 2/",
    celltype_col = "celltype",
    theta_col = "theta_pred",
    gene_col = "gene",
    target_col = "TargetCellType",
    phase_col = "Phase_group1",
    flag_col = "Module3_flag",
    n_grid = 60,
    k = 8,
    normalize_method = "LogNormalize",
    scale_factor = 10000,
    row_anno_width = 0.8,
    top_anno_height = 0.8
) {

  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(dplyr)
    library(tibble)
    library(mgcv)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
    library(Matrix)
    library(future)
    library(furrr)
    library(RColorBrewer)
  })

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  stopifnot(!is.null(rownames(expr_mat)))
  stopifnot(!is.null(colnames(expr_mat)))
  stopifnot(all(colnames(expr_mat) %in% rownames(cell_meta)))
  stopifnot(all(c(celltype_col, theta_col) %in% colnames(cell_meta)))
  stopifnot(all(c(gene_col, target_col, phase_col, flag_col) %in% colnames(ovr_df)))

  wrap_phase <- function(x) {
    ((as.numeric(x) %% (2 * pi)) + 2 * pi) %% (2 * pi)
  }

  row_zscore <- function(mat) {
    z <- t(apply(mat, 1, function(x) {
      sdx <- stats::sd(x, na.rm = TRUE)
      out <- (x - mean(x, na.rm = TRUE)) / sdx
      out[!is.finite(out)] <- 0
      out
    }))
    z <- as.matrix(z)
    rownames(z) <- rownames(mat)
    colnames(z) <- colnames(mat)
    z[!is.finite(z)] <- 0
    z
  }

  make_phase_cat_from_theta <- function(theta_grid) {
    zt_grid <- (theta_grid / (2 * pi) * 24) %% 24
    phase_cat <- floor((zt_grid - 2) / 4) * 4 + 2
    phase_cat <- phase_cat %% 24
    factor(phase_cat, levels = c(2, 6, 10, 14, 18, 22))
  }

  cyclic_smooth_matrix <- function(expr_sub, theta_sub, k = 8, n_grid = 60) {
    expr_sub <- as.matrix(expr_sub)
    theta_use <- wrap_phase(theta_sub)
    theta_grid <- seq(0, 2 * pi, length.out = n_grid + 1)[-1]

    pred_list <- furrr::future_map(
      seq_len(nrow(expr_sub)),
      function(i) {
        dat <- data.frame(
          expr = as.numeric(expr_sub[i, ]),
          theta = theta_use
        )

        fit <- try(
          mgcv::gam(
            expr ~ s(theta, bs = "cc", k = k),
            data = dat,
            method = "REML",
            knots = list(theta = c(0, 2 * pi))
          ),
          silent = TRUE
        )

        pred <- tryCatch(
          stats::predict(fit, newdata = data.frame(theta = theta_grid)),
          error = function(e) rep(NA_real_, length(theta_grid))
        )

        pred
      },
      .options = furrr::furrr_options(seed = TRUE)
    )

    smoothed_mat <- do.call(rbind, pred_list)
    smoothed_mat <- as.matrix(smoothed_mat)
    rownames(smoothed_mat) <- rownames(expr_sub)
    colnames(smoothed_mat) <- paste0("Theta_", seq_len(ncol(smoothed_mat)))

    keep <- stats::complete.cases(smoothed_mat)
    smoothed_mat <- smoothed_mat[keep, , drop = FALSE]

    list(
      mat_raw = smoothed_mat,
      mat_z = row_zscore(smoothed_mat),
      theta_grid = theta_grid,
      keep_genes = rownames(smoothed_mat)
    )
  }

  cell_meta_use <- cell_meta[colnames(expr_mat), , drop = FALSE]
  cell_meta_use[[celltype_col]] <- as.character(cell_meta_use[[celltype_col]])
  cell_meta_use[[theta_col]] <- wrap_phase(cell_meta_use[[theta_col]])

  ovr_use <- ovr_df %>%
    dplyr::filter(.data[[flag_col]]) %>%
    dplyr::mutate(
      TargetCellType = as.character(.data[[target_col]]),
      phase_wrap = wrap_phase(.data[[phase_col]])
    ) %>%
    dplyr::filter(.data[[gene_col]] %in% rownames(expr_mat))

  cell_count_df <- cell_meta_use %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::count(.data[[celltype_col]], name = "n_cells") %>%
    dplyr::rename(TargetCellType = 1) %>%
    dplyr::arrange(dplyr::desc(n_cells), TargetCellType)

  ct_levels <- cell_count_df$TargetCellType

  row_block_info <- ovr_use %>%
    dplyr::count(TargetCellType, name = "n_genes_m3") %>%
    dplyr::right_join(cell_count_df, by = "TargetCellType") %>%
    dplyr::mutate(
      n_genes_m3 = dplyr::coalesce(n_genes_m3, 0L)
    ) %>%
    dplyr::filter(n_genes_m3 > 0) %>%
    dplyr::arrange(match(TargetCellType, ct_levels))

  stopifnot(nrow(row_block_info) > 0)

  obj_all <- SeuratObject::CreateSeuratObject(
    counts = expr_mat
  )

  obj_all <- Seurat::NormalizeData(
    object = obj_all,
    normalization.method = normalize_method,
    scale.factor = scale_factor,
    verbose = FALSE
  )

  obj_all <- Seurat::ScaleData(
    object = obj_all,
    features = rownames(expr_mat),
    verbose = FALSE
  )

  mat_scale_all <- obj_all[["RNA"]]@scale.data
  mat_scale_all <- mat_scale_all[, colnames(expr_mat), drop = FALSE]

  # column block meta: width proportional to cell number
  col_block_meta <- vector("list", length(ct_levels))
  names(col_block_meta) <- ct_levels

  for (j in seq_along(ct_levels)) {
    ct_col <- ct_levels[j]

    meta_sub <- cell_meta_use %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::filter(.data[[celltype_col]] == ct_col) %>%
      dplyr::arrange(.data[[theta_col]])

    n_grid_use <- round(nrow(meta_sub) / max(cell_count_df$n_cells) * n_grid)
    n_grid_use <- max(n_grid_use, 12)

    theta_grid_use <- seq(0, 2 * pi, length.out = n_grid_use + 1)[-1]

    col_block_meta[[j]] <- list(
      celltype = ct_col,
      cells = meta_sub$cell,
      theta = meta_sub[[theta_col]],
      n_cells = nrow(meta_sub),
      n_grid_use = n_grid_use,
      theta_grid = theta_grid_use,
      phase_cat = make_phase_cat_from_theta(theta_grid_use)
    )
  }

  row_blocks <- vector("list", nrow(row_block_info))
  names(row_blocks) <- row_block_info$TargetCellType

  for (i in seq_len(nrow(row_block_info))) {

    ct_row <- row_block_info$TargetCellType[i]
    message("Building full-matrix row block for: ", ct_row)

    gene_order_df <- ovr_use %>%
      dplyr::filter(TargetCellType == ct_row) %>%
      dplyr::select(
        gene = all_of(gene_col),
        phase_wrap
      ) %>%
      dplyr::distinct(gene, .keep_all = TRUE) %>%
      dplyr::arrange(phase_wrap, gene)

    gene_order <- gene_order_df$gene

    block_piece_list <- vector("list", length(ct_levels))
    names(block_piece_list) <- ct_levels

    gene_keep_intersection <- gene_order

    for (j in seq_along(ct_levels)) {
      cells_col <- col_block_meta[[j]]$cells
      theta_vec_col <- col_block_meta[[j]]$theta
      n_grid_use <- col_block_meta[[j]]$n_grid_use

      expr_sub <- mat_scale_all[gene_order, cells_col, drop = FALSE]

      smooth_res <- cyclic_smooth_matrix(
        expr_sub = expr_sub,
        theta_sub = theta_vec_col,
        k = k,
        n_grid = n_grid_use
      )

      gene_keep_intersection <- intersect(gene_keep_intersection, smooth_res$keep_genes)
      block_piece_list[[j]] <- smooth_res
    }

    gene_keep <- gene_order[gene_order %in% gene_keep_intersection]

    mat_row_block <- do.call(
      cbind,
      lapply(block_piece_list, function(x) {
        x$mat_z[gene_keep, , drop = FALSE]
      })
    )

    row_blocks[[i]] <- list(
      celltype = ct_row,
      gene_order = gene_keep,
      mat_z = mat_row_block,
      n_genes = length(gene_keep)
    )

    gc()
  }

  big_mat <- do.call(
    rbind,
    lapply(row_blocks, function(x) x$mat_z)
  )

  row_ct <- unlist(lapply(row_blocks, function(x) rep(x$celltype, x$n_genes)))
  row_gene <- unlist(lapply(row_blocks, function(x) x$gene_order))

  col_ct <- unlist(lapply(col_block_meta, function(x) rep(x$celltype, length(x$phase_cat))))
  col_phase <- unlist(lapply(col_block_meta, function(x) as.character(x$phase_cat)))

  rownames(big_mat) <- paste0(row_ct, "__", row_gene)
  colnames(big_mat) <- paste0(col_ct, "__", seq_len(ncol(big_mat)))

  npg_palette_internal <- c(
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#91D1C2",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#DC0000",
    "#7E6148",
    "#B09C85"
  )

  ct_colors <- structure(
    npg_palette_internal[seq_along(ct_levels)],
    names = ct_levels
  )

  row_anno <- ComplexHeatmap::rowAnnotation(
    CellType = factor(row_ct, levels = ct_levels),
    col = list(CellType = ct_colors),
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(row_anno_width, "cm"),
    gap = grid::unit(0, "mm"),
    border = FALSE,
    gp = grid::gpar(col = NA, lwd = 0)
  )

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    CellType = factor(col_ct, levels = ct_levels),
    CellPhase = factor(col_phase, levels = c("2", "6", "10", "14", "18", "22")),
    col = list(
      CellType = ct_colors,
      CellPhase = c(
        "2"  = "#000000",
        "6"  = "#000000",
        "10" = "#F5F5F5",
        "14" = "#F5F5F5",
        "18" = "#F5F5F5",
        "22" = "#000000"
      )
    ),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      CellType = list(
        title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 12)
      ),
      CellPhase = list(
        title = "Cell phase",
        title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 12)
      )
    ),
    simple_anno_size = grid::unit(top_anno_height, "cm"),
    gap = grid::unit(c(1, 1), "mm"),
    border = FALSE,
    gp = grid::gpar(col = NA, lwd = 0)
  )

  col_fun <- circlize::colorRamp2(
    seq(-2, 2, length.out = 11),
    rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  )

  ht <- ComplexHeatmap::Heatmap(
    matrix = big_mat,
    name = "z-score",
    col = col_fun,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    left_annotation = row_anno,
    top_annotation = top_anno,
    row_split = factor(row_ct, levels = ct_levels),
    column_split = factor(col_ct, levels = ct_levels),
    row_gap = grid::unit(1.5, "mm"),
    column_gap = grid::unit(1.5, "mm"),
    use_raster = TRUE,
    raster_quality = 4,
    border = FALSE,
    rect_gp = grid::gpar(col = NA, lwd = 0),
    heatmap_legend_param = list(
      title = "z-score",
      title_gp = grid::gpar(fontsize = 15, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 12)
    ),
    column_title = NULL
  )

  out <- list(
    heatmap = ht,
    big_mat = big_mat,
    row_block_info = row_block_info,
    cell_count_df = cell_count_df,
    row_ct = row_ct,
    row_gene = row_gene,
    col_ct = col_ct,
    col_phase = col_phase,
    out_dir = out_dir,
    celltype_colors = ct_colors
  )

  class(out) <- "scRhythmOVRFullHeatmap"
  out
}

plot.scRhythmOVRFullHeatmap <- function(x, ...) {
  ComplexHeatmap::draw(
    x$heatmap,
    merge_legends = TRUE,
    ...
  )
  invisible(x)
}
