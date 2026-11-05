#' Plot one-vs-rest rhythmic marker heatmap
#'
#' This function builds a full-matrix heatmap for one-vs-rest target rhythmic
#' markers across multiple cell types. Genes are grouped by their target cell
#' type and ordered by target-group peak phase. Columns are grouped by cell type
#' and smoothed along inferred cell phase.
#'
#' This function accepts a gene-by-cell count matrix and uses Seurat internally
#' for log-normalization and scaling. Users do not need to provide a Seurat
#' object.
#'
#' @param count A gene-by-cell count matrix. Row names must be gene names and
#' column names must be cell names.
#' @param res_phase A data.frame returned by \code{infer_cell_phase()}. It must
#' contain \code{"cell"} and the column specified by \code{theta_by}.
#' @param ovr_result A data.frame returned by \code{infer_ovr_rhythm()}, usually
#' \code{res_ovr$ovr_result}. Current results use the column named by
#' \code{celltype_by}, \code{"Phase_target"}, and \code{"ovr_marker"}; legacy
#' column names are also accepted.
#' @param metadata A data.frame containing cell-level metadata. Cells are matched
#' to \code{count} by the \code{cell} column if available, otherwise by row
#' names.
#' @param celltype_by A character vector specifying the column name in
#' \code{metadata} used for cell-type groups. Defaults to \code{"celltype"}.
#' @param theta_by A character vector specifying the column name in
#' \code{res_phase} containing inferred cell phase in radians. Defaults to
#' \code{"theta_pred"}.
#' @param scale_factor A numeric value passed to \code{Seurat::NormalizeData()}.
#' Defaults to \code{10000}.
#' @param palette A named list specifying optional palettes for \code{celltype},
#' \code{phase}, and \code{heatmap}. If \code{NULL}, default scRhythmBayes
#' palettes are used. Defaults to \code{NULL}.
#' @param Verbose Logical; if \code{TRUE}, progress messages are printed.
#' Defaults to \code{FALSE}.
#'
#' @return An object of class \code{scRhythmOVRHeatmap}.
#'
#' @examples
#' \dontrun{
#' data("srb_example_data")
#'
#' res_phase <- infer_cell_phase(
#'   count = srb_example_data$Count,
#'   metadata = srb_example_data$meta,
#'   ZT_by = "ZT",
#'   strata_by = "celltype",
#'   Species = "mouse",
#'   setseed = 1,
#'   Verbose = FALSE
#' )
#'
#' res_ovr <- infer_ovr_rhythm(
#'   count = srb_example_data$Count,
#'   metadata = srb_example_data$meta,
#'   theta_pred = res_phase$theta_pred,
#'   r_cell = res_phase$r_cell,
#'   group_by = "celltype",
#'   ZT_by = "ZT",
#'   B = 5,
#'   n_worker = 1,
#'   setseed = 1,
#'   Verbose = FALSE
#' )
#'
#' ht_ovr <- plot_ovr_heatmap(
#'   count = srb_example_data$Count,
#'   res_phase = res_phase,
#'   ovr_result = res_ovr$ovr_result,
#'   metadata = srb_example_data$meta,
#'   celltype_by = "celltype"
#' )
#'
#' ht_ovr
#' plot(ht_ovr)
#' }
#'
#' @importFrom stats predict sd setNames complete.cases
#' @importFrom grDevices colorRampPalette
#' @importFrom grid unit
#' @export
plot_ovr_heatmap <- function(
    count,
    res_phase,
    ovr_result,
    metadata,
    celltype_by = "celltype",
    theta_by = "theta_pred",
    scale_factor = 10000,
    palette = NULL,
    Verbose = FALSE
) {
  req_pkgs <- c(
    "Seurat",
    "SeuratObject",
    "mgcv",
    "ComplexHeatmap",
    "circlize",
    "grid",
    "RColorBrewer"
  )

  for (pkg in req_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for plot_ovr_heatmap().")
    }
  }

  if (!inherits(count, c("matrix", "Matrix"))) {
    stop("count must be a matrix or Matrix object.")
  }

  if (is.null(rownames(count))) {
    stop("count must have row names as gene names.")
  }

  if (is.null(colnames(count))) {
    stop("count must have column names as cell names.")
  }

  if (!is.data.frame(res_phase)) {
    stop("res_phase must be a data.frame.")
  }

  if (!is.data.frame(ovr_result)) {
    stop("ovr_result must be a data.frame.")
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame.")
  }

  if (!"cell" %in% colnames(res_phase)) {
    stop("res_phase must contain a 'cell' column.")
  }

  if (!theta_by %in% colnames(res_phase)) {
    stop("theta_by is not found in res_phase.")
  }

  if (!is.character(celltype_by) || length(celltype_by) != 1) {
    stop("celltype_by must be a single character value.")
  }

  modern_cols <- c("gene", celltype_by, "Phase_target", "ovr_marker")
  if (all(modern_cols %in% names(ovr_result))) {
    ovr_result$TargetCellType <- ovr_result[[celltype_by]]
    ovr_result$Phase_group1 <- ovr_result$Phase_target
    ovr_result$target_specific_flag <- ovr_result$ovr_marker
  }

  required_ovr_cols <- c(
    "gene",
    "TargetCellType",
    "Phase_group1",
    "target_specific_flag"
  )

  miss_ovr_cols <- setdiff(required_ovr_cols, colnames(ovr_result))

  if (length(miss_ovr_cols) > 0) {
    stop(
      "ovr_result is missing required columns: ",
      paste(miss_ovr_cols, collapse = ", ")
    )
  }

  if (!is.numeric(scale_factor) || length(scale_factor) != 1 || scale_factor <= 0) {
    stop("scale_factor must be a single positive numeric value.")
  }

  n_grid <- 60
  k <- 8
  min_grid <- 12
  row_anno_width <- 0.8
  top_anno_height <- 0.8

  wrap_phase <- function(x) {
    ((as.numeric(x) %% (2 * pi)) + 2 * pi) %% (2 * pi)
  }

  row_zscore <- function(mat) {
    mat <- as.matrix(mat)

    if (is.null(rownames(mat))) {
      stop("Input matrix must have row names.")
    }

    if (is.null(colnames(mat))) {
      stop("Input matrix must have column names.")
    }

    if (nrow(mat) == 0 || ncol(mat) == 0) {
      return(mat)
    }

    z <- t(apply(mat, 1, function(x) {
      sdx <- stats::sd(x, na.rm = TRUE)

      if (!is.finite(sdx) || sdx <= 0) {
        out <- rep(0, length(x))
      } else {
        out <- (x - mean(x, na.rm = TRUE)) / sdx
      }

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

  get_scaled_data <- function(obj) {
    out <- tryCatch(
      SeuratObject::GetAssayData(
        object = obj,
        assay = "RNA",
        layer = "scale.data"
      ),
      error = function(e) {
        SeuratObject::GetAssayData(
          object = obj,
          assay = "RNA",
          slot = "scale.data"
        )
      }
    )

    as.matrix(out)
  }

  metadata_to_cell_df <- function(meta_use) {
    meta_use <- as.data.frame(meta_use, stringsAsFactors = FALSE)

    if (!"cell" %in% colnames(meta_use)) {
      if (!is.null(rownames(meta_use))) {
        meta_use$cell <- rownames(meta_use)
      } else {
        stop("metadata must contain a 'cell' column or rownames.")
      }
    }

    meta_use$cell <- as.character(meta_use$cell)

    if (anyDuplicated(meta_use$cell)) {
      stop("metadata contains duplicated cell identifiers.")
    }

    meta_use
  }

  cyclic_smooth_matrix <- function(expr_sub, theta_sub, n_grid_use) {
    expr_sub <- as.matrix(expr_sub)

    if (is.null(rownames(expr_sub))) {
      stop("expr_sub must have row names.")
    }

    if (is.null(colnames(expr_sub))) {
      stop("expr_sub must have column names.")
    }

    theta_use <- wrap_phase(theta_sub)
    theta_grid <- seq(0, 2 * pi, length.out = n_grid_use + 1)[-1]

    ok_cell <- is.finite(theta_use)
    expr_sub <- expr_sub[, ok_cell, drop = FALSE]
    theta_use <- theta_use[ok_cell]

    theta_names <- paste0("Theta_", seq_along(theta_grid))

    if (nrow(expr_sub) == 0 || ncol(expr_sub) < k) {
      empty_mat <- matrix(
        numeric(0),
        nrow = 0,
        ncol = length(theta_grid),
        dimnames = list(character(0), theta_names)
      )

      return(list(
        mat_raw = empty_mat,
        mat_z = empty_mat,
        theta_grid = theta_grid,
        keep_genes = character(0)
      ))
    }

    s <- mgcv::s

    pred_list <- lapply(seq_len(nrow(expr_sub)), function(i) {
      dat <- data.frame(
        expr = as.numeric(expr_sub[i, ]),
        theta = theta_use
      )

      dat <- dat[
        is.finite(dat$expr) & is.finite(dat$theta),
        ,
        drop = FALSE
      ]

      if (nrow(dat) < k) {
        return(rep(NA_real_, length(theta_grid)))
      }

      fit <- try(
        mgcv::gam(
          expr ~ s(theta, bs = "cc", k = k),
          data = dat,
          method = "REML",
          knots = list(theta = c(0, 2 * pi))
        ),
        silent = TRUE
      )

      if (inherits(fit, "try-error")) {
        return(rep(NA_real_, length(theta_grid)))
      }

      pred <- tryCatch(
        stats::predict(
          fit,
          newdata = data.frame(theta = theta_grid)
        ),
        error = function(e) rep(NA_real_, length(theta_grid))
      )

      as.numeric(pred)
    })

    smoothed_mat <- matrix(
      unlist(pred_list, use.names = FALSE),
      nrow = nrow(expr_sub),
      ncol = length(theta_grid),
      byrow = TRUE
    )

    rownames(smoothed_mat) <- rownames(expr_sub)
    colnames(smoothed_mat) <- theta_names

    keep <- complete.cases(smoothed_mat)
    smoothed_mat <- smoothed_mat[keep, , drop = FALSE]

    list(
      mat_raw = smoothed_mat,
      mat_z = row_zscore(smoothed_mat),
      theta_grid = theta_grid,
      keep_genes = rownames(smoothed_mat)
    )
  }

  meta <- metadata_to_cell_df(metadata)

  if (!celltype_by %in% colnames(meta)) {
    stop("celltype_by is not found in metadata.")
  }

  phase_df <- as.data.frame(res_phase, stringsAsFactors = FALSE)
  phase_df <- phase_df[, c("cell", theta_by), drop = FALSE]
  phase_df$cell <- as.character(phase_df$cell)

  if (anyDuplicated(phase_df$cell)) {
    stop("res_phase contains duplicated cell identifiers.")
  }

  idx_meta <- match(colnames(count), meta$cell)

  if (any(is.na(idx_meta))) {
    stop("Some cells in count are not found in metadata.")
  }

  meta <- meta[idx_meta, , drop = FALSE]

  idx_phase <- match(meta$cell, phase_df$cell)

  if (any(is.na(idx_phase))) {
    stop("Some cells in metadata/count are not found in res_phase.")
  }

  meta[[theta_by]] <- phase_df[[theta_by]][idx_phase]
  meta[[theta_by]] <- wrap_phase(meta[[theta_by]])
  meta[[celltype_by]] <- as.character(meta[[celltype_by]])

  rownames(meta) <- meta$cell
  meta <- meta[colnames(count), , drop = FALSE]

  if (any(!is.finite(meta[[theta_by]]))) {
    stop("Matched theta values contain NA/Inf.")
  }

  ovr_use <- as.data.frame(ovr_result, stringsAsFactors = FALSE)
  ovr_use <- ovr_use[
    as.logical(ovr_use$target_specific_flag) &
      as.character(ovr_use$gene) %in% rownames(count),
    ,
    drop = FALSE
  ]

  if (nrow(ovr_use) == 0) {
    stop("No target-specific genes in ovr_result are found in count.")
  }

  ovr_use$gene <- as.character(ovr_use$gene)
  ovr_use$TargetCellType <- as.character(ovr_use$TargetCellType)
  ovr_use$phase_wrap <- wrap_phase(ovr_use$Phase_group1)

  cell_count_df <- meta |>
    dplyr::count(.data[[celltype_by]], name = "n_cells") |>
    dplyr::rename(TargetCellType = 1) |>
    dplyr::arrange(dplyr::desc(.data$n_cells), .data$TargetCellType)

  ct_levels <- as.character(cell_count_df$TargetCellType)

  row_block_info <- ovr_use |>
    dplyr::count(.data$TargetCellType, name = "n_target_gene") |>
    dplyr::right_join(cell_count_df, by = "TargetCellType") |>
    dplyr::mutate(
      n_target_gene = dplyr::coalesce(.data$n_target_gene, 0L)
    ) |>
    dplyr::filter(.data$n_target_gene > 0) |>
    dplyr::arrange(match(.data$TargetCellType, ct_levels))

  if (nrow(row_block_info) == 0) {
    stop("No cell type has target-specific genes available for heatmap.")
  }

  if (isTRUE(Verbose)) {
    message(
      sprintf(
        "[plot_ovr_heatmap] %d cell types, %d target-specific genes",
        length(ct_levels),
        length(unique(ovr_use$gene))
      )
    )
  }

  obj_all <- Seurat::CreateSeuratObject(
    counts = count
  )

  obj_all <- Seurat::NormalizeData(
    object = obj_all,
    normalization.method = "LogNormalize",
    scale.factor = scale_factor,
    verbose = FALSE
  )

  obj_all <- Seurat::ScaleData(
    object = obj_all,
    features = rownames(count),
    verbose = FALSE
  )

  mat_scale_all <- get_scaled_data(obj_all)
  mat_scale_all <- mat_scale_all[, colnames(count), drop = FALSE]

  max_n_cell <- max(cell_count_df$n_cells)

  col_block_meta <- vector("list", length(ct_levels))
  names(col_block_meta) <- ct_levels

  for (j in seq_along(ct_levels)) {
    ct_col <- ct_levels[j]

    meta_sub <- meta |>
      dplyr::filter(.data[[celltype_by]] == ct_col) |>
      dplyr::arrange(.data[[theta_by]])

    n_grid_use <- round(nrow(meta_sub) / max_n_cell * n_grid)
    n_grid_use <- max(n_grid_use, min_grid)

    theta_grid_use <- seq(0, 2 * pi, length.out = n_grid_use + 1)[-1]

    col_block_meta[[j]] <- list(
      celltype = ct_col,
      cells = meta_sub$cell,
      theta = meta_sub[[theta_by]],
      n_cells = nrow(meta_sub),
      n_grid = n_grid_use,
      theta_grid = theta_grid_use,
      phase_cat = make_phase_cat_from_theta(theta_grid_use)
    )
  }

  row_blocks <- vector("list", nrow(row_block_info))
  names(row_blocks) <- row_block_info$TargetCellType

  for (i in seq_len(nrow(row_block_info))) {
    ct_row <- row_block_info$TargetCellType[i]

    if (isTRUE(Verbose)) {
      message("Building OVR row block: ", ct_row)
    }

    gene_order_df <- ovr_use |>
      dplyr::filter(.data$TargetCellType == ct_row) |>
      dplyr::select(
        gene = .data$gene,
        phase_wrap = .data$phase_wrap
      ) |>
      dplyr::distinct(.data$gene, .keep_all = TRUE) |>
      dplyr::arrange(.data$phase_wrap, .data$gene)

    gene_order <- gene_order_df$gene

    block_piece_list <- vector("list", length(ct_levels))
    names(block_piece_list) <- ct_levels

    gene_keep_intersection <- gene_order

    for (j in seq_along(ct_levels)) {
      cells_col <- col_block_meta[[j]]$cells
      theta_vec_col <- col_block_meta[[j]]$theta
      n_grid_use <- col_block_meta[[j]]$n_grid

      expr_sub <- mat_scale_all[gene_order, cells_col, drop = FALSE]

      smooth_res <- cyclic_smooth_matrix(
        expr_sub = expr_sub,
        theta_sub = theta_vec_col,
        n_grid_use = n_grid_use
      )

      gene_keep_intersection <- intersect(
        gene_keep_intersection,
        smooth_res$keep_genes
      )

      block_piece_list[[j]] <- smooth_res
    }

    gene_keep <- gene_order[gene_order %in% gene_keep_intersection]

    if (length(gene_keep) == 0) {
      row_blocks[[i]] <- list(
        celltype = ct_row,
        gene_order = character(0),
        mat_z = NULL,
        n_genes = 0
      )
      next
    }

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
  }

  row_blocks <- row_blocks[
    vapply(row_blocks, function(x) x$n_genes > 0, logical(1))
  ]

  if (length(row_blocks) == 0) {
    stop("No genes remain after smoothing.")
  }

  big_mat <- do.call(
    rbind,
    lapply(row_blocks, function(x) x$mat_z)
  )

  row_ct <- unlist(
    lapply(row_blocks, function(x) rep(x$celltype, x$n_genes)),
    use.names = FALSE
  )

  row_gene <- unlist(
    lapply(row_blocks, function(x) x$gene_order),
    use.names = FALSE
  )

  col_ct <- unlist(
    lapply(col_block_meta, function(x) rep(x$celltype, length(x$phase_cat))),
    use.names = FALSE
  )

  col_phase <- unlist(
    lapply(col_block_meta, function(x) as.character(x$phase_cat)),
    use.names = FALSE
  )

  rownames(big_mat) <- paste0(row_ct, "__", row_gene)
  colnames(big_mat) <- paste0(col_ct, "__", seq_len(ncol(big_mat)))

  celltype_pal_default <- c(
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

  if (length(ct_levels) > length(celltype_pal_default)) {
    celltype_pal_default <- grDevices::colorRampPalette(celltype_pal_default)(
      length(ct_levels)
    )
  } else {
    celltype_pal_default <- celltype_pal_default[seq_along(ct_levels)]
  }

  palette_default <- list(
    celltype = stats::setNames(celltype_pal_default, ct_levels),
    phase = c(
      "2"  = "#000000",
      "6"  = "#000000",
      "10" = "#F5F5F5",
      "14" = "#F5F5F5",
      "18" = "#F5F5F5",
      "22" = "#000000"
    ),
    heatmap = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  )

  palette_use <- palette_default

  if (!is.null(palette)) {
    palette_use[names(palette)] <- palette
  }

  row_anno <- ComplexHeatmap::rowAnnotation(
    CellType = factor(row_ct, levels = ct_levels),
    col = list(CellType = palette_use$celltype),
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
      CellType = palette_use$celltype,
      CellPhase = palette_use$phase
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
    palette_use$heatmap
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
    cell_count = cell_count_df,
    row_ct = row_ct,
    row_gene = row_gene,
    col_ct = col_ct,
    col_phase = col_phase,
    palette = palette_use,
    input = list(
      count = count,
      res_phase = res_phase,
      ovr_result = ovr_result,
      metadata = meta,
      celltype_by = celltype_by,
      theta_by = theta_by,
      scale_factor = scale_factor
    )
  )

  class(out) <- "scRhythmOVRHeatmap"

  out
}


#' Print a scRhythmOVRHeatmap object
#'
#' This S3 method prints a concise summary of a
#' \code{scRhythmOVRHeatmap} object returned by
#' \code{plot_ovr_heatmap()}.
#'
#' @param x A \code{scRhythmOVRHeatmap} object.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
plot.scRhythmOVRHeatmap <- function(x, ...) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required to draw this object.")
  }

  if (!inherits(x, "scRhythmOVRHeatmap")) {
    stop("x must be a scRhythmOVRHeatmap object.")
  }

  ComplexHeatmap::draw(
    x$heatmap,
    merge_legends = TRUE,
    ...
  )

  invisible(x)
}


#' Plot a scRhythmOVRHeatmap object
#'
#' This S3 method draws the ComplexHeatmap object stored in a
#' \code{scRhythmOVRHeatmap} object returned by \code{plot_ovr_heatmap()}.
#'
#' @param x A \code{scRhythmOVRHeatmap} object.
#' @param ... Additional arguments passed to \code{ComplexHeatmap::draw()}.
#'
#' @return Invisibly returns the drawn ComplexHeatmap object.
#'
#' @export
print.scRhythmOVRHeatmap <- function(x, ...) {
  cat("scRhythmOVRHeatmap object\n")
  cat("Genes: ", nrow(x$big_mat), "\n", sep = "")
  cat("Smoothed phase columns: ", ncol(x$big_mat), "\n", sep = "")

  if (!is.null(x$row_block_info)) {
    print(x$row_block_info)
  }

  invisible(x)
}
