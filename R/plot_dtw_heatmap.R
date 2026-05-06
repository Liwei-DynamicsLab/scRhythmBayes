#' Plot DTW-based phase-distribution similarity heatmap
#'
#' This function plots a ComplexHeatmap heatmap from a
#' \code{scRhythmDTWSimilarity} object returned by
#' \code{run_dtw_similarity()}. Rows and columns represent celltype-by-group
#' combinations. Heatmap values represent DTW-based similarity or distance
#' between circular cell-phase distributions.
#'
#' @param dtw_res A \code{scRhythmDTWSimilarity} object returned by
#' \code{run_dtw_similarity()}.
#' @param matrix_type A character value specifying which matrix to plot. Must be
#' one of \code{"similarity"} or \code{"distance"}. Defaults to
#' \code{"similarity"}.
#' @param cluster_rows A logical value indicating whether rows should be
#' clustered. Defaults to \code{TRUE}.
#' @param cluster_columns A logical value indicating whether columns should be
#' clustered. Defaults to \code{TRUE}.
#' @param show_row_names A logical value indicating whether row names should be
#' shown. Defaults to \code{TRUE}.
#' @param show_column_names A logical value indicating whether column names
#' should be shown. Defaults to \code{FALSE}.
#' @param group_colors A named character vector specifying colors for groups. If
#' \code{NULL}, the default scRhythmBayes group palette is used. Defaults to
#' \code{NULL}.
#' @param celltype_colors A named character vector specifying colors for cell
#' types. If \code{NULL}, the default scRhythmBayes discrete palette is used.
#' Defaults to \code{NULL}.
#' @param heatmap_colors A character vector specifying heatmap colors. If
#' \code{NULL}, the default scRhythmBayes DTW heatmap palette is used. Defaults
#' to \code{NULL}.
#' @param title A character vector specifying the heatmap title. If \code{NULL},
#' a default title is used. Defaults to \code{NULL}.
#' @param width A \code{\link[grid]{unit}} object specifying heatmap body width.
#' Defaults to \code{grid::unit(9, "cm")}.
#' @param height A \code{\link[grid]{unit}} object specifying heatmap body
#' height. Defaults to \code{grid::unit(9, "cm")}.
#'
#' @return An object of class \code{scRhythmDTWHeatmap}, containing:
#' \describe{
#'   \item{\code{heatmap}}{A ComplexHeatmap object.}
#'   \item{\code{matrix}}{The plotted DTW similarity or distance matrix.}
#'   \item{\code{annotation}}{Annotation data used for row and column annotations.}
#'   \item{\code{matrix_type}}{The matrix type used for plotting.}
#'   \item{\code{palette}}{Color palettes used for groups, cell types, and heatmap values.}
#'   \item{\code{params}}{Plotting parameters.}
#' }
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
#' dtw_res <- run_dtw_similarity(
#'   res_phase = res_phase,
#'   metadata = srb_example_data$meta,
#'   celltype_by = "celltype",
#'   group_by = "condition",
#'   min_cell = 10,
#'   Verbose = TRUE
#' )
#'
#' ht_dtw <- plot_dtw_heatmap(dtw_res = dtw_res)
#'
#' ht_dtw
#' plot(ht_dtw)
#' }
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom grid gpar unit
#' @importFrom stats as.dendrogram dist hclust setNames
#' @export
plot_dtw_heatmap <- function(
    dtw_res,
    matrix_type = c("similarity", "distance"),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    group_colors = NULL,
    celltype_colors = NULL,
    heatmap_colors = NULL,
    title = NULL,
    width = grid::unit(9, "cm"),
    height = grid::unit(9, "cm")
) {
  req_pkgs <- c(
    "ComplexHeatmap",
    "circlize",
    "grid",
    "RColorBrewer"
  )

  for (pkg in req_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for plot_dtw_heatmap().")
    }
  }

  if (!inherits(dtw_res, "scRhythmDTWSimilarity")) {
    stop("dtw_res must be a scRhythmDTWSimilarity object.")
  }

  matrix_type <- match.arg(matrix_type)

  if (!is.logical(cluster_rows) || length(cluster_rows) != 1) {
    stop("cluster_rows must be TRUE or FALSE.")
  }

  if (!is.logical(cluster_columns) || length(cluster_columns) != 1) {
    stop("cluster_columns must be TRUE or FALSE.")
  }

  if (!is.logical(show_row_names) || length(show_row_names) != 1) {
    stop("show_row_names must be TRUE or FALSE.")
  }

  if (!is.logical(show_column_names) || length(show_column_names) != 1) {
    stop("show_column_names must be TRUE or FALSE.")
  }

  if (!is.null(title) && (!is.character(title) || length(title) != 1)) {
    stop("title must be NULL or a single character value.")
  }

  mat_use <- if (matrix_type == "similarity") {
    dtw_res$similarity
  } else {
    dtw_res$distance
  }

  mat_use <- as.matrix(mat_use)

  if (nrow(mat_use) == 0 || ncol(mat_use) == 0) {
    stop("The selected DTW matrix is empty.")
  }

  if (is.null(rownames(mat_use)) || is.null(colnames(mat_use))) {
    stop("The selected DTW matrix must have row names and column names.")
  }

  annotation_df <- as.data.frame(
    dtw_res$annotation,
    stringsAsFactors = FALSE
  )

  required_cols <- c("combo", "cell_type", "Group")
  miss_cols <- setdiff(required_cols, colnames(annotation_df))

  if (length(miss_cols) > 0) {
    stop(
      "dtw_res$annotation is missing required columns: ",
      paste(miss_cols, collapse = ", ")
    )
  }

  rownames(annotation_df) <- annotation_df$combo

  if (!all(rownames(mat_use) %in% rownames(annotation_df))) {
    stop("Some matrix row names are not found in dtw_res$annotation.")
  }

  if (!all(colnames(mat_use) %in% rownames(annotation_df))) {
    stop("Some matrix column names are not found in dtw_res$annotation.")
  }

  annotation_row <- annotation_df[rownames(mat_use), , drop = FALSE]
  annotation_col <- annotation_df[colnames(mat_use), , drop = FALSE]

  group_levels <- as.character(dtw_res$params$group_levels)
  celltype_levels <- as.character(dtw_res$params$celltype_levels)

  group_levels <- group_levels[group_levels %in% annotation_df$Group]
  celltype_levels <- celltype_levels[celltype_levels %in% annotation_df$cell_type]

  if (length(group_levels) == 0) {
    group_levels <- unique(annotation_df$Group)
  }

  if (length(celltype_levels) == 0) {
    celltype_levels <- unique(annotation_df$cell_type)
  }

  annotation_row$Group <- factor(annotation_row$Group, levels = group_levels)
  annotation_col$Group <- factor(annotation_col$Group, levels = group_levels)

  annotation_row$cell_type <- factor(
    annotation_row$cell_type,
    levels = celltype_levels
  )

  annotation_col$cell_type <- factor(
    annotation_col$cell_type,
    levels = celltype_levels
  )

  match_named_colors <- function(cols, levels_use, arg_name) {
    levels_use <- as.character(levels_use)

    if (is.null(cols)) {
      return(NULL)
    }

    if (is.null(names(cols))) {
      if (length(cols) < length(levels_use)) {
        stop(arg_name, " must contain enough colors for all levels.")
      }

      cols <- stats::setNames(cols[seq_along(levels_use)], levels_use)
    } else {
      miss <- setdiff(levels_use, names(cols))

      if (length(miss) > 0) {
        stop(
          arg_name,
          " is missing colors for levels: ",
          paste(miss, collapse = ", ")
        )
      }

      cols <- cols[levels_use]
    }

    cols
  }

  if (is.null(group_colors)) {
    group_colors <- .srb_group_colors(group_levels)
  } else {
    group_colors <- match_named_colors(
      cols = group_colors,
      levels_use = group_levels,
      arg_name = "group_colors"
    )
  }

  if (is.null(celltype_colors)) {
    celltype_colors <- .srb_discrete_colors(celltype_levels)
  } else {
    celltype_colors <- match_named_colors(
      cols = celltype_colors,
      levels_use = celltype_levels,
      arg_name = "celltype_colors"
    )
  }

  if (is.null(heatmap_colors)) {
    if (matrix_type == "similarity") {
      heatmap_colors <- rev(
        grDevices::colorRampPalette(
          RColorBrewer::brewer.pal(9, "YlGnBu")
        )(100)
      )
    } else {
      heatmap_colors <- grDevices::colorRampPalette(
        rev(RColorBrewer::brewer.pal(9, "YlGnBu"))
      )(100)
    }
  }

  if (!is.character(heatmap_colors) || length(heatmap_colors) < 2) {
    stop("heatmap_colors must be a character vector with at least two colors.")
  }

  finite_vals <- mat_use[is.finite(mat_use)]

  if (length(finite_vals) == 0) {
    stop("The selected DTW matrix contains no finite values.")
  }

  if (matrix_type == "similarity") {
    col_fun <- circlize::colorRamp2(
      seq(0, 1, length.out = length(heatmap_colors)),
      heatmap_colors
    )

    legend_title <- "Similarity"

    if (is.null(title)) {
      title <- "DTW-based Similarity"
    }
  } else {
    val_min <- min(finite_vals, na.rm = TRUE)
    val_max <- max(finite_vals, na.rm = TRUE)

    if (!is.finite(val_min) || !is.finite(val_max) || val_max <= val_min) {
      val_min <- 0
      val_max <- 1
    }

    col_fun <- circlize::colorRamp2(
      seq(val_min, val_max, length.out = length(heatmap_colors)),
      heatmap_colors
    )

    legend_title <- "Distance"

    if (is.null(title)) {
      title <- "DTW Distance"
    }
  }

  top_ha <- ComplexHeatmap::HeatmapAnnotation(
    df = data.frame(
      Group = annotation_col$Group,
      cell_type = annotation_col$cell_type,
      check.names = FALSE
    ),
    col = list(
      Group = group_colors,
      cell_type = celltype_colors
    ),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 13, fontface = "bold"),
    annotation_name_side = "right",
    simple_anno_size = grid::unit(7, "mm"),
    annotation_height = grid::unit(c(7, 7), "mm"),
    gap = grid::unit(1.2, "mm"),
    border = FALSE,
    gp = grid::gpar(col = NA),
    annotation_legend_param = list(
      Group = list(
        title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 12),
        grid_width = grid::unit(0.55, "cm"),
        grid_height = grid::unit(0.55, "cm")
      ),
      cell_type = list(
        title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 12),
        grid_width = grid::unit(0.55, "cm"),
        grid_height = grid::unit(0.55, "cm")
      )
    )
  )

  left_ha <- ComplexHeatmap::rowAnnotation(
    df = data.frame(
      Group = annotation_row$Group,
      cell_type = annotation_row$cell_type,
      check.names = FALSE
    ),
    col = list(
      Group = group_colors,
      cell_type = celltype_colors
    ),
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(7, "mm"),
    annotation_width = grid::unit(c(7, 7), "mm"),
    gap = grid::unit(1.2, "mm"),
    border = FALSE,
    gp = grid::gpar(col = NA)
  )

  row_cluster_use <- cluster_rows
  col_cluster_use <- cluster_columns

  if (isTRUE(cluster_rows) && nrow(mat_use) > 1) {
    row_cluster_use <- stats::as.dendrogram(
      stats::hclust(stats::dist(mat_use))
    )
  }

  if (isTRUE(cluster_columns) && ncol(mat_use) > 1) {
    col_cluster_use <- stats::as.dendrogram(
      stats::hclust(stats::dist(t(mat_use)))
    )
  }

  ht <- ComplexHeatmap::Heatmap(
    mat_use,
    name = legend_title,
    col = col_fun,
    cluster_rows = row_cluster_use,
    cluster_columns = col_cluster_use,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    column_names_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    row_names_side = "right",
    top_annotation = top_ha,
    left_annotation = left_ha,
    border = FALSE,
    rect_gp = grid::gpar(col = NA, lwd = 0),
    heatmap_legend_param = list(
      title = legend_title,
      title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 13),
      legend_height = grid::unit(5, "cm")
    ),
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 24, fontface = "bold"),
    width = width,
    height = height
  )

  out <- list(
    heatmap = ht,
    matrix = mat_use,
    annotation = annotation_df,
    matrix_type = matrix_type,
    palette = list(
      group = group_colors,
      celltype = celltype_colors,
      heatmap = heatmap_colors
    ),
    params = list(
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      title = title,
      width = width,
      height = height
    )
  )

  class(out) <- "scRhythmDTWHeatmap"

  out
}


#' Plot a scRhythmDTWHeatmap object
#'
#' This S3 method draws a \code{scRhythmDTWHeatmap} object returned by
#' \code{plot_dtw_heatmap()}.
#'
#' @param x A \code{scRhythmDTWHeatmap} object returned by
#' \code{plot_dtw_heatmap()}.
#' @param merge_legends A logical value indicating whether heatmap and
#' annotation legends should be merged. Defaults to \code{TRUE}.
#' @param heatmap_legend_side A character value specifying the heatmap legend
#' side. Defaults to \code{"right"}.
#' @param annotation_legend_side A character value specifying the annotation
#' legend side. Defaults to \code{"right"}.
#' @param padding A \code{\link[grid]{unit}} object specifying padding passed to
#' \code{ComplexHeatmap::draw()}. Defaults to
#' \code{grid::unit(c(6, 4, 4, 4), "mm")}.
#' @param ... Additional arguments passed to \code{ComplexHeatmap::draw()}.
#'
#' @return Invisibly returns the drawn ComplexHeatmap object.
#'
#' @importFrom grid unit
#' @export
plot.scRhythmDTWHeatmap <- function(
    x,
    merge_legends = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = grid::unit(c(6, 4, 4, 4), "mm"),
    ...
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required to draw this object.")
  }

  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required to draw this object.")
  }

  if (!inherits(x, "scRhythmDTWHeatmap")) {
    stop("x must be a scRhythmDTWHeatmap object.")
  }

  if (!is.logical(merge_legends) || length(merge_legends) != 1) {
    stop("merge_legends must be TRUE or FALSE.")
  }

  ht_drawn <- ComplexHeatmap::draw(
    x$heatmap,
    merge_legends = merge_legends,
    heatmap_legend_side = heatmap_legend_side,
    annotation_legend_side = annotation_legend_side,
    padding = padding,
    ...
  )

  invisible(ht_drawn)
}


#' Print a scRhythmDTWHeatmap object
#'
#' @param x A \code{scRhythmDTWHeatmap} object.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.scRhythmDTWHeatmap <- function(x, ...) {
  if (!inherits(x, "scRhythmDTWHeatmap")) {
    stop("x must be a scRhythmDTWHeatmap object.")
  }

  cat("scRhythmDTWHeatmap object\n")
  cat("Matrix type:", x$matrix_type, "\n")
  cat("Matrix dimension:", nrow(x$matrix), "x", ncol(x$matrix), "\n")

  invisible(x)
}
