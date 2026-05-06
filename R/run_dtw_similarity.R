#' Run DTW-based phase-distribution similarity analysis
#'
#' This function compares cell-phase distributions across celltype-by-group
#' combinations using dynamic time warping (DTW). For each combination, inferred
#' cell phases are converted into a circular histogram. Pairwise DTW distances
#' between histograms are then converted to a 0-to-1 similarity matrix.
#'
#' Cells are matched between \code{res_phase} and \code{metadata} by the
#' \code{cell} column.
#'
#' @param res_phase A data.frame returned by \code{infer_cell_phase()}.
#' It must contain \code{"cell"} and the column specified by \code{theta_by}.
#' @param metadata A data.frame containing cell-level metadata.
#' @param theta_by A character vector specifying the column name in
#' \code{res_phase} that contains inferred cell phase in radians. Defaults to
#' \code{"theta_pred"}.
#' @param celltype_by A character vector specifying the column name in
#' \code{metadata} that defines cell types. Defaults to \code{"celltype"}.
#' @param group_by A character vector specifying the column name in
#' \code{metadata} that defines groups or conditions. Defaults to
#' \code{"condition"}.
#' @param group_levels A character vector specifying group order. If \code{NULL},
#' groups are sorted automatically. Defaults to \code{NULL}.
#' @param celltype_levels A character vector specifying cell-type order. If
#' \code{NULL}, cell types are sorted automatically. Defaults to \code{NULL}.
#' @param n_bin A numeric value specifying the number of circular histogram bins.
#' Defaults to \code{24}.
#' @param min_cell A numeric value specifying the minimum number of cells
#' required for a celltype-by-group combination. Defaults to \code{10}.
#' @param normalize_distance A logical value indicating whether DTW distance
#' should be normalized by warping path length. Defaults to \code{TRUE}.
#' @param Verbose A logical value indicating whether progress messages should be
#' printed. Defaults to \code{FALSE}.
#'
#' @return An object of class \code{scRhythmDTWSimilarity}, containing:
#' \describe{
#'   \item{\code{data}}{Merged cell-level phase and metadata table.}
#'   \item{\code{histogram}}{Circular histogram table for each celltype-by-group combination.}
#'   \item{\code{combo}}{Cell-count summary for retained combinations.}
#'   \item{\code{annotation}}{Annotation table matched to rows and columns of the matrices.}
#'   \item{\code{distance}}{Pairwise DTW distance matrix.}
#'   \item{\code{similarity}}{Pairwise 0-to-1 DTW similarity matrix.}
#'   \item{\code{params}}{Parameters used for the analysis.}
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
#' dtw_res
#' dim(dtw_res$similarity)
#' dtw_res$similarity[1:3, 1:3]
#' head(dtw_res$histogram)
#' }
#'
#' @importFrom graphics hist
#' @export
run_dtw_similarity <- function(
    res_phase,
    metadata,
    theta_by = "theta_pred",
    celltype_by = "celltype",
    group_by = "condition",
    group_levels = NULL,
    celltype_levels = NULL,
    n_bin = 24,
    min_cell = 10,
    normalize_distance = TRUE,
    Verbose = FALSE
) {
  if (!requireNamespace("dtw", quietly = TRUE)) {
    stop("Package 'dtw' is required for run_dtw_similarity().")
  }

  if (!is.data.frame(res_phase)) {
    stop("res_phase must be a data.frame.")
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

  if (!is.character(theta_by) || length(theta_by) != 1) {
    stop("theta_by must be a single character value.")
  }

  if (!is.character(celltype_by) || length(celltype_by) != 1) {
    stop("celltype_by must be a single character value.")
  }

  if (!is.character(group_by) || length(group_by) != 1) {
    stop("group_by must be a single character value.")
  }

  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!"cell" %in% colnames(meta)) {
    if (!is.null(rownames(meta))) {
      meta$cell <- rownames(meta)
    } else {
      stop("metadata must contain a 'cell' column or rownames.")
    }
  }

  if (!celltype_by %in% colnames(meta)) {
    stop("celltype_by is not found in metadata.")
  }

  if (!group_by %in% colnames(meta)) {
    stop("group_by is not found in metadata.")
  }

  if (!is.numeric(n_bin) || length(n_bin) != 1 ||
      !is.finite(n_bin) || n_bin < 4) {
    stop("n_bin must be a numeric value greater than or equal to 4.")
  }

  n_bin <- as.integer(n_bin)

  if (!is.numeric(min_cell) || length(min_cell) != 1 ||
      !is.finite(min_cell) || min_cell < 1) {
    stop("min_cell must be a positive numeric value.")
  }

  min_cell <- as.integer(min_cell)

  if (!is.logical(normalize_distance) || length(normalize_distance) != 1) {
    stop("normalize_distance must be TRUE or FALSE.")
  }

  if (!is.logical(Verbose) || length(Verbose) != 1) {
    stop("Verbose must be TRUE or FALSE.")
  }

  wrap_phase <- function(x) {
    ((as.numeric(x) %% (2 * pi)) + 2 * pi) %% (2 * pi)
  }

  hist_breaks <- seq(0, 2 * pi, length.out = n_bin + 1)
  hist_start <- hist_breaks[seq_len(n_bin)]
  hist_end <- hist_breaks[seq_len(n_bin) + 1]
  hist_mid <- (hist_start + hist_end) / 2

  phase_to_hist <- function(phase_rad) {
    phase_rad <- as.numeric(phase_rad)
    phase_rad <- phase_rad[is.finite(phase_rad)]
    phase_rad <- wrap_phase(phase_rad)

    if (length(phase_rad) == 0) {
      return(rep(NA_real_, n_bin))
    }

    h <- graphics::hist(
      phase_rad,
      breaks = hist_breaks,
      plot = FALSE,
      include.lowest = TRUE,
      right = FALSE
    )

    total <- sum(h$counts)

    if (!is.finite(total) || total <= 0) {
      rep(NA_real_, n_bin)
    } else {
      as.numeric(h$counts) / total
    }
  }

  idx <- match(res_phase$cell, meta$cell)

  if (any(is.na(idx))) {
    stop("Some cells in res_phase are not found in metadata.")
  }

  df_phase <- data.frame(
    cell = as.character(res_phase$cell),
    phase_rad = wrap_phase(res_phase[[theta_by]]),
    cell_type = as.character(meta[[celltype_by]][idx]),
    Group = as.character(meta[[group_by]][idx]),
    stringsAsFactors = FALSE
  )

  df_phase$phase_hour <- ((df_phase$phase_rad / (2 * pi)) * 24) %% 24

  df_phase <- df_phase[
    is.finite(df_phase$phase_rad) &
      !is.na(df_phase$cell_type) &
      nzchar(df_phase$cell_type) &
      !is.na(df_phase$Group) &
      nzchar(df_phase$Group),
    ,
    drop = FALSE
  ]

  if (nrow(df_phase) == 0) {
    stop("No finite phase values are available for DTW similarity analysis.")
  }

  if (is.null(group_levels)) {
    group_levels <- sort(unique(df_phase$Group))
  } else {
    group_levels <- as.character(group_levels)
  }

  if (is.null(celltype_levels)) {
    celltype_levels <- sort(unique(df_phase$cell_type))
  } else {
    celltype_levels <- as.character(celltype_levels)
  }

  if (length(group_levels) == 0) {
    stop("group_levels must contain at least one group.")
  }

  if (length(celltype_levels) == 0) {
    stop("celltype_levels must contain at least one cell type.")
  }

  df_phase <- df_phase[
    df_phase$Group %in% group_levels &
      df_phase$cell_type %in% celltype_levels,
    ,
    drop = FALSE
  ]

  if (nrow(df_phase) == 0) {
    stop("No cells remain after applying group_levels and celltype_levels.")
  }

  df_phase$Group <- factor(df_phase$Group, levels = group_levels)
  df_phase$cell_type <- factor(df_phase$cell_type, levels = celltype_levels)
  df_phase$combo <- paste(df_phase$cell_type, df_phase$Group, sep = " | ")

  combo_df <- as.data.frame(
    table(df_phase$cell_type, df_phase$Group),
    stringsAsFactors = FALSE
  )

  colnames(combo_df) <- c("cell_type", "Group", "n_cell")

  combo_df <- combo_df[
    combo_df$n_cell >= min_cell,
    ,
    drop = FALSE
  ]

  combo_df$cell_type <- factor(combo_df$cell_type, levels = celltype_levels)
  combo_df$Group <- factor(combo_df$Group, levels = group_levels)

  combo_df <- combo_df[
    order(combo_df$cell_type, combo_df$Group),
    ,
    drop = FALSE
  ]

  combo_df$combo <- paste(combo_df$cell_type, combo_df$Group, sep = " | ")

  if (nrow(combo_df) < 2) {
    stop(
      "Fewer than two celltype-by-group combinations remain after min_cell filtering."
    )
  }

  combo_names <- as.character(combo_df$combo)

  df_phase <- df_phase[
    df_phase$combo %in% combo_names,
    ,
    drop = FALSE
  ]

  if (isTRUE(Verbose)) {
    message("Use ", nrow(df_phase), " cells.")
    message("Use ", length(combo_names), " celltype-by-group combinations.")
  }

  hist_list <- lapply(combo_names, function(combo_now) {
    phase_to_hist(
      phase_rad = df_phase$phase_rad[df_phase$combo == combo_now]
    )
  })

  names(hist_list) <- combo_names

  hist_df <- do.call(
    rbind,
    lapply(seq_along(hist_list), function(i) {
      combo_now <- names(hist_list)[i]
      combo_meta <- combo_df[combo_df$combo == combo_now, , drop = FALSE]

      data.frame(
        combo = combo_now,
        cell_type = as.character(combo_meta$cell_type),
        Group = as.character(combo_meta$Group),
        bin = seq_len(n_bin),
        bin_start = hist_start,
        bin_end = hist_end,
        bin_mid = hist_mid,
        prop = as.numeric(hist_list[[i]]),
        stringsAsFactors = FALSE
      )
    })
  )

  distance_mat <- matrix(
    NA_real_,
    nrow = length(combo_names),
    ncol = length(combo_names),
    dimnames = list(combo_names, combo_names)
  )

  for (i in seq_along(combo_names)) {
    distance_mat[i, i] <- 0

    if (i < length(combo_names)) {
      for (j in seq.int(i + 1L, length(combo_names))) {
        fit <- dtw::dtw(
          hist_list[[i]],
          hist_list[[j]],
          keep = TRUE
        )

        distance_use <- as.numeric(fit$distance)

        if (isTRUE(normalize_distance)) {
          path_len <- length(fit$index1)

          if (is.finite(path_len) && path_len > 0) {
            distance_use <- distance_use / path_len
          }
        }

        distance_mat[i, j] <- distance_use
        distance_mat[j, i] <- distance_use
      }
    }
  }

  off_diag <- distance_mat[row(distance_mat) != col(distance_mat)]
  off_diag <- off_diag[is.finite(off_diag)]

  if (length(off_diag) == 0) {
    stop("Cannot compute DTW similarity from the distance matrix.")
  }

  d_min <- min(off_diag, na.rm = TRUE)
  d_max <- max(off_diag, na.rm = TRUE)

  if (!is.finite(d_min) || !is.finite(d_max) || d_max <= d_min) {
    similarity_mat <- matrix(
      1,
      nrow = nrow(distance_mat),
      ncol = ncol(distance_mat),
      dimnames = dimnames(distance_mat)
    )
  } else {
    similarity_mat <- 1 - (distance_mat - d_min) / (d_max - d_min)
    similarity_mat[similarity_mat < 0] <- 0
    similarity_mat[similarity_mat > 1] <- 1
    diag(similarity_mat) <- 1
  }

  annotation_df <- combo_df[, c("combo", "cell_type", "Group", "n_cell")]
  annotation_df$cell_type <- as.character(annotation_df$cell_type)
  annotation_df$Group <- as.character(annotation_df$Group)

  rownames(annotation_df) <- annotation_df$combo

  out <- list(
    data = df_phase,
    histogram = hist_df,
    combo = combo_df,
    annotation = annotation_df,
    distance = distance_mat,
    similarity = similarity_mat,
    params = list(
      theta_by = theta_by,
      celltype_by = celltype_by,
      group_by = group_by,
      group_levels = group_levels,
      celltype_levels = celltype_levels,
      n_bin = n_bin,
      min_cell = min_cell,
      normalize_distance = normalize_distance
    )
  )

  class(out) <- "scRhythmDTWSimilarity"

  out
}


#' Print a scRhythmDTWSimilarity object
#'
#' @param x A \code{scRhythmDTWSimilarity} object returned by
#' \code{run_dtw_similarity()}.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.scRhythmDTWSimilarity <- function(x, ...) {
  if (!inherits(x, "scRhythmDTWSimilarity")) {
    stop("x must be a scRhythmDTWSimilarity object.")
  }

  cat("scRhythmDTWSimilarity object\n")
  cat("Combinations:", nrow(x$annotation), "\n")
  cat("Bins:", x$params$n_bin, "\n")
  cat("Groups:", paste(x$params$group_levels, collapse = ", "), "\n")
  cat("Cell types:", paste(x$params$celltype_levels, collapse = ", "), "\n\n")

  print(x$combo)

  invisible(x)
}
