#' Plot paired module heatmaps across cell groups
#'
#' This function builds paired module-level heatmaps showing smoothed rhythmic
#' expression patterns along inferred cell phase. For each selected module,
#' genes are ordered by module-specific reference phase and expression is
#' smoothed separately in two specified cell groups.
#'
#' This function accepts a gene-by-cell count matrix and uses Seurat internally
#' for log-normalization and scaling. Users do not need to provide a Seurat
#' object.
#'
#' @param count A gene-by-cell count matrix. Row names must be gene names and
#' column names must be cell names.
#' @param res_phase A data.frame returned by \code{infer_cell_phase()}. It must
#' contain \code{"cell"} and the column specified by \code{theta_by}.
#' @param gene_meta A data.frame returned by \code{infer_rhythm_module()}.
#' It must contain \code{"gene"} and \code{"module_map"} columns, as well as
#' group-specific rhythm parameter columns such as \code{"Phase_A"},
#' \code{"Amplitude_A"}, and \code{"FDR_A"}. The alternative column style
#' \code{"phi_A"}, \code{"A_A"}, and \code{"q_A"} is also supported.
#' @param metadata A data.frame containing cell-level metadata. Cells are matched
#' to \code{count} by the \code{cell} column if available, otherwise by rownames.
#' @param group_by A character vector specifying the column name in
#' \code{metadata} used to define the two paired cell groups.
#' @param group1 A character vector specifying the first group to plot. This
#' should usually match \code{group1} used in \code{infer_rhythm_module()}.
#' @param group2 A character vector specifying the second group to plot. This
#' should usually match \code{group2} used in \code{infer_rhythm_module()}.
#' @param theta_by A character vector specifying the column name in
#' \code{res_phase} containing inferred cell phase in radians. Defaults to
#' \code{"theta_pred"}.
#' @param module_ids An integer vector specifying modules to plot. Defaults to
#' \code{0:7}.
#' @param phase_levels A numeric vector specifying displayed phase groups in
#' hours. Defaults to \code{c(2, 6, 10, 14, 18, 22)}.
#' @param scale_factor A numeric value passed to \code{Seurat::NormalizeData()}.
#' Defaults to \code{10000}.
#' @param palette A named list specifying optional palettes for \code{phase},
#' \code{group}, and \code{heatmap}. If \code{NULL}, black/white phase labels,
#' default scRhythmBayes group colors, and a blue-red heatmap palette are used.
#' Defaults to \code{NULL}.
#' @param Verbose Logical; if \code{TRUE}, progress messages are printed.
#'
#' @return An object of class \code{scRhythmModuleHeatmap}.
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
#' cell_use <- as.character(srb_example_data$meta$celltype) == "1"
#'
#' res_param <- infer_gene_param(
#'   count = srb_example_data$Count[, cell_use, drop = FALSE],
#'   metadata = srb_example_data$meta[cell_use, , drop = FALSE],
#'   theta_pred = res_phase$theta_pred[cell_use],
#'   r_cell = res_phase$r_cell[cell_use],
#'   strata_by = "condition",
#'   B = 100,
#'   n_worker = 1,
#'   setseed = 1,
#'   Verbose = FALSE
#' )
#'
#' res_module <- infer_rhythm_module(
#'   param_wide = res_param$param_wide,
#'   compare_by = "condition",
#'   group1 = "A",
#'   group2 = "B",
#'   out_dir = NULL,
#'   Verbose = FALSE
#' )
#'
#' hm <- plot_module_heatmap(
#'   count = srb_example_data$Count[, cell_use, drop = FALSE],
#'   res_phase = res_phase[cell_use, , drop = FALSE],
#'   gene_meta = res_module$module_result,
#'   metadata = srb_example_data$meta[cell_use, , drop = FALSE],
#'   group_by = "condition",
#'   group1 = "A",
#'   group2 = "B"
#' )
#'
#' hm
#' plot(hm, module_id = 0)
#' }
#'
#' @importFrom stats predict sd setNames complete.cases
#' @importFrom grDevices colorRampPalette
#' @importFrom grid unit
#' @export
plot_module_heatmap <- function(
    count,
    res_phase,
    gene_meta,
    metadata,
    group_by,
    group1,
    group2,
    theta_by = "theta_pred",
    module_ids = 0:7,
    phase_levels = c(2, 6, 10, 14, 18, 22),
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
      stop("Package '", pkg, "' is required for plot_module_heatmap().")
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

  if (!is.data.frame(gene_meta)) {
    stop("gene_meta must be a data.frame.")
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

  if (!"gene" %in% colnames(gene_meta)) {
    stop("gene_meta must contain a 'gene' column.")
  }

  if (!"module_map" %in% colnames(gene_meta)) {
    stop("gene_meta must contain a 'module_map' column.")
  }

  if (!is.character(group_by) || length(group_by) != 1) {
    stop("group_by must be a single character value.")
  }

  if (!is.character(group1) || length(group1) != 1) {
    stop("group1 must be a single character value.")
  }

  if (!is.character(group2) || length(group2) != 1) {
    stop("group2 must be a single character value.")
  }

  if (identical(group1, group2)) {
    stop("group1 and group2 must be different.")
  }

  k <- 8
  n_grid <- 60

  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!"cell" %in% colnames(meta)) {
    if (!is.null(rownames(meta))) {
      meta$cell <- rownames(meta)
    } else {
      stop("metadata must contain a 'cell' column or rownames.")
    }
  }

  if (!group_by %in% colnames(meta)) {
    stop("group_by is not found in metadata.")
  }

  phase_df <- as.data.frame(res_phase, stringsAsFactors = FALSE)
  phase_df <- phase_df[, c("cell", theta_by), drop = FALSE]
  colnames(phase_df)[colnames(phase_df) == theta_by] <- "theta_pred_use_input"

  idx_meta <- match(colnames(count), meta$cell)

  if (any(is.na(idx_meta))) {
    stop("Some cells in count are not found in metadata.")
  }

  meta <- meta[idx_meta, , drop = FALSE]

  idx_phase <- match(meta$cell, phase_df$cell)

  if (any(is.na(idx_phase))) {
    stop("Some cells in metadata/count are not found in res_phase.")
  }

  meta[[theta_by]] <- phase_df$theta_pred_use_input[idx_phase]
  rownames(meta) <- meta$cell
  meta <- meta[colnames(count), , drop = FALSE]

  group_values <- unique(as.character(meta[[group_by]]))

  if (!group1 %in% group_values) {
    stop("group1 is not found in metadata[[group_by]].")
  }

  if (!group2 %in% group_values) {
    stop("group2 is not found in metadata[[group_by]].")
  }

  group_levels_use <- c(group1, group2)

  keep_group <- as.character(meta[[group_by]]) %in% group_levels_use
  meta <- meta[keep_group, , drop = FALSE]
  count <- count[, rownames(meta), drop = FALSE]

  if (ncol(count) == 0) {
    stop("No cells remain after filtering to group1 and group2.")
  }

  module_ids <- unique(as.integer(module_ids))
  module_ids <- module_ids[!is.na(module_ids)]

  if (length(module_ids) == 0) {
    stop("module_ids must contain at least one integer module ID.")
  }

  phase_levels <- as.numeric(phase_levels)
  phase_levels <- phase_levels[is.finite(phase_levels)]

  if (length(phase_levels) == 0) {
    stop("phase_levels must contain at least one finite value.")
  }

  wrap_phase <- function(x) {
    ((as.numeric(x) %% (2 * pi)) + 2 * pi) %% (2 * pi)
  }

  make_phase_cat_from_theta <- function(theta_grid) {
    zt_grid <- (theta_grid / (2 * pi) * 24) %% 24
    phase_cat <- floor((zt_grid - 2) / 4) * 4 + 2
    phase_cat <- phase_cat %% 24
    factor(phase_cat, levels = phase_levels)
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

  find_param_col <- function(prefix_candidates, group_name, all_cols) {
    candidate_vec <- unlist(lapply(prefix_candidates, function(px) {
      c(
        paste0(px, group_name),
        paste0(px, make.names(group_name)),
        paste0(px, gsub("\\.", "_", group_name)),
        paste0(px, gsub("_", ".", group_name))
      )
    }))

    candidate_vec <- unique(candidate_vec)
    hit <- candidate_vec[candidate_vec %in% all_cols]

    if (length(hit) == 0) {
      NA_character_
    } else {
      hit[1]
    }
  }

  infer_param_cols <- function(gene_meta, group_levels_use) {
    all_cols <- colnames(gene_meta)

    param_col_map <- lapply(group_levels_use, function(g) {
      list(
        phase = find_param_col(
          prefix_candidates = c("Phase_", "phi_"),
          group_name = g,
          all_cols = all_cols
        ),
        amplitude = find_param_col(
          prefix_candidates = c("Amplitude_", "A_"),
          group_name = g,
          all_cols = all_cols
        ),
        fdr = find_param_col(
          prefix_candidates = c("FDR_", "q_"),
          group_name = g,
          all_cols = all_cols
        )
      )
    })

    names(param_col_map) <- group_levels_use

    cols_need <- unlist(lapply(param_col_map[group_levels_use], function(x) {
      c(x$phase, x$amplitude, x$fdr)
    }))

    ok <- all(!is.na(cols_need)) && all(cols_need %in% colnames(gene_meta))

    if (!ok) {
      stop(
        paste0(
          "Cannot infer parameter columns from group1/group2.\n",
          "Expected columns similar to:\n",
          "  Phase_", group1, ", Amplitude_", group1, ", FDR_", group1, "\n",
          "  Phase_", group2, ", Amplitude_", group2, ", FDR_", group2, "\n",
          "or:\n",
          "  phi_", group1, ", A_", group1, ", q_", group1, "\n",
          "  phi_", group2, ", A_", group2, ", q_", group2, "\n",
          "Available columns in gene_meta:\n",
          paste(colnames(gene_meta), collapse = ", ")
        )
      )
    }

    param_col_map
  }

  get_module_title <- function(gene_df, module_use) {
    if ("module_name" %in% colnames(gene_df)) {
      title_use <- unique(as.character(gene_df$module_name))
      title_use <- title_use[!is.na(title_use) & nzchar(title_use)]

      if (length(title_use) > 0) {
        return(title_use[1])
      }
    }

    paste0("Module ", module_use)
  }

  param_col_map_use <- infer_param_cols(
    gene_meta = gene_meta,
    group_levels_use = group_levels_use
  )

  palette_default <- list(
    phase = stats::setNames(
      rep(c("#000000", "#F5F5F5"), length.out = length(phase_levels)),
      as.character(phase_levels)
    ),
    group = .srb_group_colors(group_levels_use),
    heatmap = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  )

  palette_use <- palette_default

  if (!is.null(palette)) {
    palette_use[names(palette)] <- palette
  }

  cyclic_smooth_matrix <- function(expr_sub, theta_sub) {
    expr_sub <- as.matrix(expr_sub)

    if (is.null(rownames(expr_sub))) {
      stop("expr_sub must have row names.")
    }

    if (is.null(colnames(expr_sub))) {
      stop("expr_sub must have column names.")
    }

    theta_use <- wrap_phase(theta_sub)
    theta_grid <- seq(0, 2 * pi, length.out = n_grid + 1)[-1]

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

  build_top_annotation <- function(theta_grid, group_name, side = c("left", "right")) {
    side <- match.arg(side)

    phase_cat <- make_phase_cat_from_theta(theta_grid)
    group_cat <- factor(
      rep(group_name, length(theta_grid)),
      levels = group_levels_use
    )

    anno_df <- data.frame(
      phase = phase_cat,
      group = group_cat,
      check.names = FALSE
    )

    if (side == "left") {
      colnames(anno_df) <- c("Cell phase_left", "Group_left")
      anno_names_use <- c("Cell phase_left", "Group_left")
      show_name_use <- c(FALSE, FALSE)
      show_legend_use <- c(FALSE, FALSE)
      legend_param_use <- list()
    } else {
      colnames(anno_df) <- c("Cell phase", "Group")
      anno_names_use <- c("Cell phase", "Group")
      show_name_use <- c(FALSE, FALSE)
      show_legend_use <- c(TRUE, TRUE)
      legend_param_use <- list(
        `Cell phase` = list(
          title_gp = grid::gpar(fontsize = 15, fontface = "bold"),
          labels_gp = grid::gpar(fontsize = 13),
          grid_width = grid::unit(0.7, "cm"),
          grid_height = grid::unit(0.7, "cm")
        ),
        Group = list(
          title_gp = grid::gpar(fontsize = 15, fontface = "bold"),
          labels_gp = grid::gpar(fontsize = 13),
          grid_width = grid::unit(0.7, "cm"),
          grid_height = grid::unit(0.7, "cm")
        )
      )
    }

    ComplexHeatmap::HeatmapAnnotation(
      df = anno_df,
      col = structure(
        list(
          palette_use$phase,
          palette_use$group
        ),
        names = anno_names_use
      ),
      show_annotation_name = show_name_use,
      show_legend = show_legend_use,
      annotation_name_gp = grid::gpar(fontsize = 16, fontface = "bold"),
      annotation_legend_param = legend_param_use,
      simple_anno_size = grid::unit(1.2, "cm"),
      gap = grid::unit(0, "mm"),
      border = FALSE,
      gp = grid::gpar(col = NA, lwd = 0)
    )
  }

  col_fun <- circlize::colorRamp2(
    seq(-2, 2, length.out = 11),
    palette_use$heatmap
  )

  order_ref_map <- c(
    "0" = group1,
    "1" = group1,
    "2" = group2,
    "3" = group1,
    "4" = group1,
    "5" = group1,
    "6" = group1,
    "7" = group1
  )

  module_list <- list()
  summary_list <- list()

  for (module_use in module_ids) {
    if (isTRUE(Verbose)) {
      message("Building module ", module_use)
    }

    gene_df <- gene_meta[
      as.integer(gene_meta$module_map) == module_use,
      ,
      drop = FALSE
    ]

    module_title_use <- get_module_title(gene_df, module_use)

    gene_pool <- unique(as.character(gene_df$gene))
    gene_pool <- intersect(gene_pool, rownames(count))

    if (length(gene_pool) == 0) {
      module_list[[paste0("module", module_use)]] <- list(
        module_id = module_use,
        status = "no_genes",
        gene_df = gene_df,
        title = module_title_use
      )

      summary_list[[paste0("module", module_use)]] <- data.frame(
        module = module_use,
        module_name = module_title_use,
        n_gene_input = 0,
        n_gene_final = 0,
        ref_group = NA_character_,
        status = "no_genes",
        stringsAsFactors = FALSE
      )

      next
    }

    ref_group_use <- unname(order_ref_map[as.character(module_use)])

    if (is.na(ref_group_use) || length(ref_group_use) == 0) {
      ref_group_use <- group1
    }

    ref_phase_col_use <- param_col_map_use[[ref_group_use]]$phase

    gene_order_df <- data.frame(
      gene = as.character(gene_df$gene),
      phase_ref = as.numeric(gene_df[[ref_phase_col_use]]),
      stringsAsFactors = FALSE
    )

    gene_order_df$phase_wrap <- wrap_phase(gene_order_df$phase_ref)

    gene_order_df <- gene_order_df[
      gene_order_df$gene %in% gene_pool,
      ,
      drop = FALSE
    ]

    gene_order_df <- gene_order_df[
      order(is.na(gene_order_df$phase_wrap), gene_order_df$phase_wrap),
      ,
      drop = FALSE
    ]

    gene_order_df <- gene_order_df[
      !duplicated(gene_order_df$gene),
      ,
      drop = FALSE
    ]

    gene_order <- gene_order_df$gene

    if (length(gene_order) == 0) {
      module_list[[paste0("module", module_use)]] <- list(
        module_id = module_use,
        status = "no_ordered_genes",
        gene_df = gene_df,
        ref_group = ref_group_use,
        title = module_title_use
      )

      summary_list[[paste0("module", module_use)]] <- data.frame(
        module = module_use,
        module_name = module_title_use,
        n_gene_input = length(gene_pool),
        n_gene_final = 0,
        ref_group = ref_group_use,
        status = "no_ordered_genes",
        stringsAsFactors = FALSE
      )

      next
    }

    group_res_list <- lapply(group_levels_use, function(group_use) {
      meta_sub <- meta
      meta_sub$theta_pred_use <- wrap_phase(meta_sub[[theta_by]])

      meta_sub <- meta_sub[
        as.character(meta_sub[[group_by]]) == group_use,
        ,
        drop = FALSE
      ]

      meta_sub <- meta_sub[
        order(meta_sub$theta_pred_use),
        ,
        drop = FALSE
      ]

      cell_order <- rownames(meta_sub)
      theta_vec <- meta_sub$theta_pred_use

      if (length(cell_order) == 0) {
        return(list(
          group_name = group_use,
          status = "no_cells",
          mat_raw = NULL,
          mat_z = NULL,
          theta_grid = NULL,
          gene_keep = character(0)
        ))
      }

      obj_sub <- Seurat::CreateSeuratObject(
        counts = count[gene_order, cell_order, drop = FALSE]
      )

      obj_sub <- Seurat::NormalizeData(
        object = obj_sub,
        normalization.method = "LogNormalize",
        scale.factor = scale_factor,
        verbose = FALSE
      )

      obj_sub <- Seurat::ScaleData(
        object = obj_sub,
        features = gene_order,
        verbose = FALSE
      )

      mat_z <- get_scaled_data(obj_sub)
      mat_plot <- mat_z[gene_order, cell_order, drop = FALSE]

      smooth_res <- cyclic_smooth_matrix(
        expr_sub = mat_plot,
        theta_sub = theta_vec
      )

      gene_keep <- gene_order[gene_order %in% smooth_res$keep_genes]

      list(
        group_name = group_use,
        status = "ok",
        mat_raw = smooth_res$mat_raw[gene_keep, , drop = FALSE],
        mat_z = smooth_res$mat_z[gene_keep, , drop = FALSE],
        theta_grid = smooth_res$theta_grid,
        gene_keep = gene_keep
      )
    })

    names(group_res_list) <- group_levels_use

    left_res <- group_res_list[[group1]]
    right_res <- group_res_list[[group2]]

    if (left_res$status != "ok" || right_res$status != "ok") {
      module_list[[paste0("module", module_use)]] <- list(
        module_id = module_use,
        status = "smooth_failed",
        gene_df = gene_df,
        ref_group = ref_group_use,
        title = module_title_use
      )

      summary_list[[paste0("module", module_use)]] <- data.frame(
        module = module_use,
        module_name = module_title_use,
        n_gene_input = length(gene_pool),
        n_gene_final = 0,
        ref_group = ref_group_use,
        status = "smooth_failed",
        stringsAsFactors = FALSE
      )

      next
    }

    gene_keep <- intersect(left_res$gene_keep, right_res$gene_keep)

    if (length(gene_keep) == 0) {
      module_list[[paste0("module", module_use)]] <- list(
        module_id = module_use,
        status = "no_common_genes_after_smooth",
        gene_df = gene_df,
        ref_group = ref_group_use,
        title = module_title_use
      )

      summary_list[[paste0("module", module_use)]] <- data.frame(
        module = module_use,
        module_name = module_title_use,
        n_gene_input = length(gene_pool),
        n_gene_final = 0,
        ref_group = ref_group_use,
        status = "no_common_genes_after_smooth",
        stringsAsFactors = FALSE
      )

      next
    }

    mat_left_raw <- left_res$mat_raw[gene_keep, , drop = FALSE]
    mat_right_raw <- right_res$mat_raw[gene_keep, , drop = FALSE]

    mat_joint <- cbind(mat_left_raw, mat_right_raw)
    mat_joint_z <- row_zscore(mat_joint)

    mat_left_z <- mat_joint_z[
      ,
      seq_len(ncol(mat_left_raw)),
      drop = FALSE
    ]

    mat_right_z <- mat_joint_z[
      ,
      ncol(mat_left_raw) + seq_len(ncol(mat_right_raw)),
      drop = FALSE
    ]

    rownames(mat_left_z) <- gene_keep
    rownames(mat_right_z) <- gene_keep

    ha_left <- build_top_annotation(
      theta_grid = left_res$theta_grid,
      group_name = group1,
      side = "left"
    )

    ha_right <- build_top_annotation(
      theta_grid = right_res$theta_grid,
      group_name = group2,
      side = "right"
    )

    ht_left <- ComplexHeatmap::Heatmap(
      mat_left_z,
      name = "z-score-left",
      col = col_fun,
      top_annotation = ha_left,
      show_row_names = FALSE,
      show_column_names = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      use_raster = TRUE,
      raster_quality = 4,
      border = FALSE,
      rect_gp = grid::gpar(col = NA, lwd = 0),
      heatmap_legend_param = list(
        title = "z-score",
        title_gp = grid::gpar(fontsize = 17, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 16)
      ),
      column_title = group1,
      column_title_gp = grid::gpar(fontsize = 23, fontface = "bold"),
      show_heatmap_legend = FALSE
    )

    ht_right <- ComplexHeatmap::Heatmap(
      mat_right_z,
      name = "z-score",
      col = col_fun,
      top_annotation = ha_right,
      show_row_names = FALSE,
      show_column_names = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      use_raster = TRUE,
      raster_quality = 4,
      border = FALSE,
      rect_gp = grid::gpar(col = NA, lwd = 0),
      heatmap_legend_param = list(
        title = "z-score",
        title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 14)
      ),
      column_title = group2,
      column_title_gp = grid::gpar(fontsize = 23, fontface = "bold"),
      show_heatmap_legend = TRUE
    )

    ht_pair <- ht_left + ht_right

    module_list[[paste0("module", module_use)]] <- list(
      module_id = module_use,
      status = "ok",
      module_type = "paired_heatmap",
      ref_group = ref_group_use,
      gene_df = gene_df,
      gene_pool = gene_pool,
      gene_keep = gene_keep,
      mat_left_z = mat_left_z,
      mat_right_z = mat_right_z,
      theta_grid_left = left_res$theta_grid,
      theta_grid_right = right_res$theta_grid,
      phase_annotation_left = make_phase_cat_from_theta(left_res$theta_grid),
      phase_annotation_right = make_phase_cat_from_theta(right_res$theta_grid),
      phase_annotation_name_left = "Cell phase_left",
      phase_annotation_name_right = "Cell phase",
      left_group = group1,
      right_group = group2,
      heatmap = ht_pair,
      title = module_title_use
    )

    summary_list[[paste0("module", module_use)]] <- data.frame(
      module = module_use,
      module_name = module_title_use,
      n_gene_input = length(gene_pool),
      n_gene_final = length(gene_keep),
      ref_group = ref_group_use,
      status = "ok",
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    input = list(
      count = count,
      res_phase = res_phase,
      gene_meta = gene_meta,
      metadata = meta,
      group_by = group_by,
      group1 = group1,
      group2 = group2,
      theta_by = theta_by,
      param_col_map = param_col_map_use,
      module_ids = module_ids,
      phase_levels = phase_levels,
      scale_factor = scale_factor
    ),
    palette = palette_use,
    group_levels = group_levels_use,
    module_list = module_list,
    summary = do.call(rbind, summary_list)
  )

  class(out) <- "scRhythmModuleHeatmap"

  out
}


#' Plot a scRhythmModuleHeatmap object
#'
#' This S3 method draws one module-specific heatmap from a
#' \code{scRhythmModuleHeatmap} object returned by
#' \code{plot_module_heatmap()}.
#'
#' @param x A \code{scRhythmModuleHeatmap} object returned by
#' \code{plot_module_heatmap()}.
#' @param module_id An integer specifying which module to draw.
#' @param title A character value specifying the plot title. If \code{NULL},
#' the stored module title is used. Defaults to \code{NULL}.
#' @param draw_separators Logical; if \code{TRUE}, separators are drawn between
#' phase blocks. Defaults to \code{TRUE}.
#' @param separator_lwd A numeric value specifying separator line width.
#' Defaults to \code{1.2}.
#' @param ... Additional arguments passed to \code{ComplexHeatmap::draw()}.
#'
#' @return Invisibly returns the drawn ComplexHeatmap object.
#'
#' @export
plot.scRhythmModuleHeatmap <- function(
    x,
    module_id,
    title = NULL,
    draw_separators = TRUE,
    separator_lwd = 1.2,
    ...
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required to draw this object.")
  }

  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required to draw this object.")
  }

  if (!inherits(x, "scRhythmModuleHeatmap")) {
    stop("x must be a scRhythmModuleHeatmap object.")
  }

  if (!is.numeric(module_id) && !is.character(module_id)) {
    stop("module_id must be an integer or character value.")
  }

  if (length(module_id) != 1) {
    stop("module_id must be a single value.")
  }

  module_name <- paste0("module", module_id)
  obj <- x$module_list[[module_name]]

  if (is.null(obj)) {
    stop("module ", module_id, " is not available in this object.")
  }

  if (!identical(obj$status, "ok")) {
    stop(
      "module ",
      module_id,
      " is not drawable. Current status: ",
      obj$status
    )
  }

  if (is.null(title)) {
    title_use <- obj$title

    if (is.null(title_use) || is.na(title_use) || !nzchar(title_use)) {
      title_use <- paste0("M", module_id)
    }
  } else {
    if (!is.character(title) || length(title) != 1) {
      stop("title must be NULL or a single character value.")
    }

    title_use <- title
  }

  if (!is.logical(draw_separators) || length(draw_separators) != 1) {
    stop("draw_separators must be a single logical value.")
  }

  if (
    !is.numeric(separator_lwd) ||
    length(separator_lwd) != 1 ||
    !is.finite(separator_lwd) ||
    separator_lwd < 0
  ) {
    stop("separator_lwd must be a single non-negative numeric value.")
  }

  ht_drawn <- ComplexHeatmap::draw(
    obj$heatmap,
    merge_legends = TRUE,
    column_title = title_use,
    column_title_gp = grid::gpar(fontsize = 28, fontface = "bold"),
    ...
  )

  if (isTRUE(draw_separators)) {
    phase_num_left <- as.character(obj$phase_annotation_left)
    phase_num_right <- as.character(obj$phase_annotation_right)

    boundary_idx_left <- which(
      phase_num_left[-1] != phase_num_left[-length(phase_num_left)]
    )

    boundary_idx_right <- which(
      phase_num_right[-1] != phase_num_right[-length(phase_num_right)]
    )

    x_pos_left <- boundary_idx_left / pmax(length(phase_num_left), 1)
    x_pos_right <- boundary_idx_right / pmax(length(phase_num_right), 1)

    ComplexHeatmap::decorate_annotation(obj$phase_annotation_name_left, {
      for (xline in x_pos_left) {
        grid::grid.lines(
          x = grid::unit(c(xline, xline), "npc"),
          y = grid::unit(c(0, 1), "npc"),
          gp = grid::gpar(col = "white", lwd = separator_lwd)
        )
      }
    })

    ComplexHeatmap::decorate_annotation(obj$phase_annotation_name_right, {
      for (xline in x_pos_right) {
        grid::grid.lines(
          x = grid::unit(c(xline, xline), "npc"),
          y = grid::unit(c(0, 1), "npc"),
          gp = grid::gpar(col = "white", lwd = separator_lwd)
        )
      }
    })
  }

  invisible(ht_drawn)
}


#' Print a scRhythmModuleHeatmap object
#'
#' @param x A \code{scRhythmModuleHeatmap} object.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.scRhythmModuleHeatmap <- function(x, ...) {
  cat("scRhythmModuleHeatmap object\n")
  print(x$summary)
  invisible(x)
}


#' Plot module-level rhythm remodeling in polar coordinates
#'
#' This function draws a polar scatter plot for one rhythm-state module. Each
#' gene is represented by two points, one for each compared group. The angular
#' coordinate represents peak phase, the radial coordinate represents rhythmic
#' amplitude, and line segments connect the same gene across the two groups.
#'
#' Genes are selected from a module result table returned by
#' \code{infer_rhythm_module()}. The function automatically detects
#' group-specific rhythm parameter columns such as \code{"Phase_A"},
#' \code{"Amplitude_A"}, and \code{"FDR_A"}. The alternative column style
#' \code{"phi_A"}, \code{"A_A"}, and \code{"q_A"} is also supported.
#'
#' @param gene_meta A data.frame returned by \code{infer_rhythm_module()}.
#' It must contain \code{"gene"} and \code{"module_map"} columns, as well as
#' group-specific rhythm parameter columns.
#' @param module_id An integer or character value specifying the module to plot.
#' @param group1 A character vector specifying the first group. This should
#' usually match \code{group1} used in \code{infer_rhythm_module()}.
#' @param group2 A character vector specifying the second group. This should
#' usually match \code{group2} used in \code{infer_rhythm_module()}.
#' @param q_cut A numeric value specifying the FDR cutoff used to define
#' rhythmic genes. Defaults to \code{0.05}.
#' @param amp_cut A numeric value specifying the minimum amplitude required for
#' selected genes. Defaults to \code{0.10}.
#' @param n_gene An integer specifying the number of top-ranked genes to plot.
#' Defaults to \code{30}.
#' @param n_label An integer specifying the number of top-ranked genes to label.
#' Defaults to \code{8}.
#' @param pt_size A numeric value specifying point size. Defaults to
#' \code{3.2}.
#' @param label_size A numeric value specifying gene-label size. Defaults to
#' \code{3.2}.
#' @param title A character value specifying the plot title. If \code{NULL},
#' the title is set to \code{"M"} followed by \code{module_id}. Defaults to
#' \code{NULL}.
#' @param palette A named list specifying optional palettes for \code{module}
#' and \code{group}. If \code{NULL}, default module colors and default
#' scRhythmBayes group colors are used. Defaults to \code{NULL}.
#' @param Verbose Logical; if \code{TRUE}, progress messages are printed.
#'
#' @return Invisibly returns a list containing the processed plotting data,
#' paired gene table, selected labels, top genes, parameter-column map, and a
#' ggplot object.
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
#' cell_use <- as.character(srb_example_data$meta$celltype) == "1"
#'
#' res_param <- infer_gene_param(
#'   count = srb_example_data$Count[, cell_use, drop = FALSE],
#'   metadata = srb_example_data$meta[cell_use, , drop = FALSE],
#'   theta_pred = res_phase$theta_pred[cell_use],
#'   r_cell = res_phase$r_cell[cell_use],
#'   strata_by = "condition",
#'   B = 100,
#'   n_worker = 1,
#'   setseed = 1,
#'   Verbose = FALSE
#' )
#'
#' res_module <- infer_rhythm_module(
#'   param_wide = res_param$param_wide,
#'   compare_by = "condition",
#'   group1 = "A",
#'   group2 = "B",
#'   out_dir = NULL,
#'   Verbose = FALSE
#' )
#'
#' p_polar <- plot_module_polar(
#'   gene_meta = res_module$module_result,
#'   module_id = 7,
#'   group1 = "A",
#'   group2 = "B"
#' )
#'
#' p_polar$plot
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_path geom_segment
#' @importFrom ggplot2 geom_point geom_text annotate coord_equal
#' @importFrom ggplot2 scale_fill_manual labs theme element_text element_blank
#' @importFrom ggplot2 margin
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats sd quantile
#' @export
plot_module_polar <- function(
    gene_meta,
    module_id,
    group1,
    group2,
    q_cut = 0.05,
    amp_cut = 0.10,
    n_gene = 30,
    n_label = 8,
    pt_size = 3.2,
    label_size = 3.2,
    title = NULL,
    palette = NULL,
    Verbose = FALSE
) {
  req_pkgs <- c(
    "dplyr",
    "tidyr",
    "tibble",
    "ggplot2",
    "ggrepel",
    "scales"
  )

  for (pkg in req_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for plot_module_polar().")
    }
  }

  if (!is.data.frame(gene_meta)) {
    stop("gene_meta must be a data.frame.")
  }

  if (!"gene" %in% colnames(gene_meta)) {
    stop("gene_meta must contain a 'gene' column.")
  }

  if (!"module_map" %in% colnames(gene_meta)) {
    stop("gene_meta must contain a 'module_map' column.")
  }

  if (!is.numeric(module_id) && !is.character(module_id)) {
    stop("module_id must be an integer or character value.")
  }

  if (length(module_id) != 1) {
    stop("module_id must be a single value.")
  }

  if (!is.character(group1) || length(group1) != 1) {
    stop("group1 must be a single character value.")
  }

  if (!is.character(group2) || length(group2) != 1) {
    stop("group2 must be a single character value.")
  }

  if (identical(group1, group2)) {
    stop("group1 and group2 must be different.")
  }

  if (!is.numeric(q_cut) || length(q_cut) != 1 || !is.finite(q_cut)) {
    stop("q_cut must be a single finite numeric value.")
  }

  if (!is.numeric(amp_cut) || length(amp_cut) != 1 || !is.finite(amp_cut)) {
    stop("amp_cut must be a single finite numeric value.")
  }

  if (!is.numeric(n_gene) || length(n_gene) != 1 || n_gene < 1) {
    stop("n_gene must be a positive integer.")
  }

  if (!is.numeric(n_label) || length(n_label) != 1 || n_label < 0) {
    stop("n_label must be a non-negative integer.")
  }

  if (!is.numeric(pt_size) || length(pt_size) != 1 || pt_size <= 0) {
    stop("pt_size must be a positive numeric value.")
  }

  if (!is.numeric(label_size) || length(label_size) != 1 || label_size <= 0) {
    stop("label_size must be a positive numeric value.")
  }

  if (!is.null(title) && (!is.character(title) || length(title) != 1)) {
    stop("title must be NULL or a single character value.")
  }

  n_gene <- as.integer(n_gene)
  n_label <- as.integer(n_label)
  module_id <- as.character(module_id)
  group_levels_use <- c(group1, group2)

  module_label_map <- c(
    "0" = "M0 Non-rhythmic",
    "1" = "M1 Stable Rhythm",
    "2" = "M2 Rhythm Gain",
    "3" = "M3 Rhythm Loss",
    "4" = "M4 Amp Enhanced",
    "5" = "M5 Amp Attenuated",
    "6" = "M6 Phase Advanced",
    "7" = "M7 Phase Delayed"
  )

  module_pal <- c(
    "M0 Non-rhythmic"   = "#BDBDBD",
    "M1 Stable Rhythm"  = "#AFC0D6",
    "M2 Rhythm Gain"    = "#EDC0A0",
    "M3 Rhythm Loss"    = "#7F9DBF",
    "M4 Amp Enhanced"   = "#DF8F79",
    "M5 Amp Attenuated" = "#E7EDF4",
    "M6 Phase Advanced" = "#FAEFD9",
    "M7 Phase Delayed"  = "#4973A8"
  )

  group_pal <- .srb_group_colors(group_levels_use)

  palette_default <- list(
    module = module_pal,
    group = group_pal
  )

  palette_use <- palette_default

  if (!is.null(palette)) {
    palette_use[names(palette)] <- palette
  }

  find_param_col <- function(prefix_candidates, group_name, all_cols) {
    candidate_vec <- unlist(lapply(prefix_candidates, function(px) {
      c(
        paste0(px, group_name),
        paste0(px, make.names(group_name)),
        paste0(px, gsub("\\.", "_", group_name)),
        paste0(px, gsub("_", ".", group_name))
      )
    }))

    candidate_vec <- unique(candidate_vec)
    hit <- candidate_vec[candidate_vec %in% all_cols]

    if (length(hit) == 0) {
      NA_character_
    } else {
      hit[1]
    }
  }

  infer_param_cols <- function(gene_meta, group_levels_use) {
    all_cols <- colnames(gene_meta)

    param_col_map <- lapply(group_levels_use, function(g) {
      list(
        phase = find_param_col(
          prefix_candidates = c("Phase_", "phi_"),
          group_name = g,
          all_cols = all_cols
        ),
        amplitude = find_param_col(
          prefix_candidates = c("Amplitude_", "A_"),
          group_name = g,
          all_cols = all_cols
        ),
        fdr = find_param_col(
          prefix_candidates = c("FDR_", "q_"),
          group_name = g,
          all_cols = all_cols
        )
      )
    })

    names(param_col_map) <- group_levels_use

    cols_need <- unlist(lapply(param_col_map[group_levels_use], function(x) {
      c(x$phase, x$amplitude, x$fdr)
    }))

    ok <- all(!is.na(cols_need)) && all(cols_need %in% colnames(gene_meta))

    if (!ok) {
      stop(
        paste0(
          "Cannot infer parameter columns from group1/group2.\n",
          "Expected columns similar to:\n",
          "  Phase_", group1, ", Amplitude_", group1, ", FDR_", group1, "\n",
          "  Phase_", group2, ", Amplitude_", group2, ", FDR_", group2, "\n",
          "or:\n",
          "  phi_", group1, ", A_", group1, ", q_", group1, "\n",
          "  phi_", group2, ", A_", group2, ", q_", group2, "\n",
          "Available columns in gene_meta:\n",
          paste(colnames(gene_meta), collapse = ", ")
        )
      )
    }

    param_col_map
  }

  param_col_map_use <- infer_param_cols(
    gene_meta = gene_meta,
    group_levels_use = group_levels_use
  )

  if (isTRUE(Verbose)) {
    message(
      "Using phase columns: ",
      paste(
        vapply(param_col_map_use, `[[`, character(1), "phase"),
        collapse = ", "
      )
    )

    message(
      "Using amplitude columns: ",
      paste(
        vapply(param_col_map_use, `[[`, character(1), "amplitude"),
        collapse = ", "
      )
    )

    message(
      "Using FDR columns: ",
      paste(
        vapply(param_col_map_use, `[[`, character(1), "fdr"),
        collapse = ", "
      )
    )
  }

  gene_meta_use <- gene_meta[
    as.character(gene_meta$module_map) == module_id,
    ,
    drop = FALSE
  ]

  if (nrow(gene_meta_use) == 0) {
    stop("No genes found for module_id = ", module_id, ".")
  }

  module_label_use <- unname(module_label_map[module_id])

  if (is.na(module_label_use) || length(module_label_use) == 0) {
    module_label_use <- paste0("Module ", module_id)
  }

  if (is.null(title)) {
    title_use <- paste0("M", module_id)
  } else {
    title_use <- title
  }

  wrap_phase <- function(x) {
    ((as.numeric(x) %% (2 * pi)) + 2 * pi) %% (2 * pi)
  }

  safe_scale <- function(x) {
    x <- as.numeric(x)
    sdx <- stats::sd(x, na.rm = TRUE)

    if (!is.finite(sdx) || sdx <= 0) {
      return(rep(0, length(x)))
    }

    out <- as.numeric(scale(x))
    out[!is.finite(out)] <- 0
    out
  }

  make_group_df <- function(group_name) {
    data.frame(
      gene = as.character(gene_meta_use$gene),
      group = rep(group_name, nrow(gene_meta_use)),
      Phase = as.numeric(
        gene_meta_use[[param_col_map_use[[group_name]]$phase]]
      ),
      Amplitude = as.numeric(
        gene_meta_use[[param_col_map_use[[group_name]]$amplitude]]
      ),
      FDR = as.numeric(
        gene_meta_use[[param_col_map_use[[group_name]]$fdr]]
      ),
      module_id = as.character(gene_meta_use$module_map),
      module_label = rep(module_label_use, nrow(gene_meta_use)),
      stringsAsFactors = FALSE
    )
  }

  plot_df0 <- do.call(
    rbind,
    lapply(group_levels_use, make_group_df)
  )

  plot_df0$phase_wrap <- wrap_phase(plot_df0$Phase)

  plot_df0 <- plot_df0[
    is.finite(plot_df0$Phase) |
      is.finite(plot_df0$Amplitude) |
      is.finite(plot_df0$FDR),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df0) == 0) {
    stop("No valid phase, amplitude, or FDR values are available for plotting.")
  }

  pair_base <- plot_df0 |>
    dplyr::select(
      .data$gene,
      .data$group,
      .data$phase_wrap,
      .data$Amplitude,
      .data$FDR,
      .data$module_id,
      .data$module_label
    ) |>
    tidyr::pivot_wider(
      names_from = .data$group,
      values_from = c(
        .data$phase_wrap,
        .data$Amplitude,
        .data$FDR
      )
    )

  phase_1 <- paste0("phase_wrap_", group1)
  phase_2 <- paste0("phase_wrap_", group2)
  amp_1 <- paste0("Amplitude_", group1)
  amp_2 <- paste0("Amplitude_", group2)
  fdr_1 <- paste0("FDR_", group1)
  fdr_2 <- paste0("FDR_", group2)

  pair_base <- pair_base |>
    dplyr::mutate(
      sig_1 = is.finite(.data[[fdr_1]]) & .data[[fdr_1]] < q_cut,
      sig_2 = is.finite(.data[[fdr_2]]) & .data[[fdr_2]] < q_cut,
      both_sig = .data$sig_1 & .data$sig_2,
      one_sig = xor(.data$sig_1, .data$sig_2),
      min_amp = pmin(.data[[amp_1]], .data[[amp_2]], na.rm = TRUE),
      max_amp = pmax(.data[[amp_1]], .data[[amp_2]], na.rm = TRUE),
      keep_gene = dplyr::case_when(
        .data$module_id %in% c("1", "4", "5", "6", "7") ~
          .data$both_sig &
          is.finite(.data$min_amp) &
          .data$min_amp > amp_cut,
        .data$module_id %in% c("2", "3") ~
          .data$one_sig &
          is.finite(.data$max_amp) &
          .data$max_amp > amp_cut,
        TRUE ~ FALSE
      )
    )

  gene_keep <- pair_base |>
    dplyr::filter(.data$keep_gene) |>
    dplyr::pull(.data$gene)

  if (length(gene_keep) == 0) {
    stop(
      paste0(
        "No genes remain after filtering for module_id = ",
        module_id,
        ". For M1/M4/M5/M6/M7, both groups must pass q_cut and amp_cut. ",
        "For M2/M3, exactly one group must pass q_cut and at least one group ",
        "must pass amp_cut. Consider increasing q_cut or decreasing amp_cut."
      )
    )
  }

  plot_df1 <- plot_df0 |>
    dplyr::filter(.data$gene %in% gene_keep)

  amp_cap <- stats::quantile(
    plot_df1$Amplitude,
    probs = 0.98,
    na.rm = TRUE,
    names = FALSE
  )

  if (!is.finite(amp_cap) || amp_cap <= 0) {
    stop("Cannot compute a valid amplitude scale from selected genes.")
  }

  plot_df2 <- plot_df1 |>
    dplyr::mutate(
      sig_flag = is.finite(.data$FDR) & .data$FDR < q_cut,
      amp_clip = pmin(.data$Amplitude, amp_cap),
      radius_raw = scales::rescale(
        .data$amp_clip,
        to = c(0.18, 0.95),
        from = c(0, amp_cap)
      ),
      radius = dplyr::case_when(
        .data$module_id %in% c("2", "3") & !.data$sig_flag ~ 0.12,
        TRUE ~ .data$radius_raw
      ),
      x = .data$radius * sin(.data$phase_wrap),
      y = .data$radius * cos(.data$phase_wrap)
    )

  pair_df <- plot_df2 |>
    dplyr::select(
      .data$gene,
      .data$group,
      .data$phase_wrap,
      .data$radius,
      .data$x,
      .data$y,
      .data$Amplitude,
      .data$FDR,
      .data$sig_flag,
      .data$module_id,
      .data$module_label
    ) |>
    tidyr::pivot_wider(
      names_from = .data$group,
      values_from = c(
        .data$phase_wrap,
        .data$radius,
        .data$x,
        .data$y,
        .data$Amplitude,
        .data$FDR,
        .data$sig_flag
      )
    )

  x_1 <- paste0("x_", group1)
  x_2 <- paste0("x_", group2)
  y_1 <- paste0("y_", group1)
  y_2 <- paste0("y_", group2)

  pair_df <- pair_df |>
    dplyr::mutate(
      delta_phase = atan2(
        sin(.data[[phase_2]] - .data[[phase_1]]),
        cos(.data[[phase_2]] - .data[[phase_1]])
      ),
      abs_delta_phase = abs(.data$delta_phase),
      delta_amp = .data[[amp_2]] - .data[[amp_1]],
      abs_delta_amp = abs(.data$delta_amp),
      mean_amp = (.data[[amp_1]] + .data[[amp_2]]) / 2,
      max_amp = pmax(.data[[amp_1]], .data[[amp_2]], na.rm = TRUE),
      sum_logfdr = -log10(.data[[fdr_1]] + 1e-300) +
        -log10(.data[[fdr_2]] + 1e-300)
    )

  pair_df <- pair_df |>
    dplyr::group_by(.data$module_id) |>
    dplyr::mutate(
      z_abs_delta_phase = safe_scale(.data$abs_delta_phase),
      z_abs_delta_amp = safe_scale(.data$abs_delta_amp),
      z_mean_amp = safe_scale(.data$mean_amp),
      z_max_amp = safe_scale(.data$max_amp),
      z_sum_logfdr = safe_scale(.data$sum_logfdr)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      rank_score = dplyr::case_when(
        .data$module_id == "1" ~
          .data$z_mean_amp + .data$z_sum_logfdr -
          .data$z_abs_delta_phase - .data$z_abs_delta_amp,
        .data$module_id == "2" ~
          .data$z_max_amp,
        .data$module_id == "3" ~
          .data$z_max_amp,
        .data$module_id == "4" ~
          .data$z_abs_delta_amp + 0.5 * .data$z_mean_amp +
          0.5 * .data$z_sum_logfdr - 0.5 * .data$z_abs_delta_phase,
        .data$module_id == "5" ~
          .data$z_abs_delta_amp + 0.5 * .data$z_mean_amp +
          0.5 * .data$z_sum_logfdr - 0.5 * .data$z_abs_delta_phase,
        .data$module_id == "6" ~
          .data$z_abs_delta_phase + 0.5 * .data$z_mean_amp +
          0.5 * .data$z_sum_logfdr,
        .data$module_id == "7" ~
          .data$z_abs_delta_phase + 0.5 * .data$z_mean_amp +
          0.5 * .data$z_sum_logfdr,
        TRUE ~
          .data$z_mean_amp
      )
    )

  top_gene_df <- pair_df |>
    dplyr::slice_max(
      order_by = .data$rank_score,
      n = n_gene,
      with_ties = FALSE
    ) |>
    dplyr::select(
      .data$gene,
      .data$module_id,
      .data$module_label,
      .data$rank_score,
      .data$abs_delta_phase,
      .data$abs_delta_amp,
      .data$mean_amp,
      .data$max_amp
    )

  plot_df_top <- plot_df2 |>
    dplyr::inner_join(
      top_gene_df |> dplyr::select(.data$gene, .data$module_id),
      by = c("gene", "module_id")
    )

  pair_df_top <- pair_df |>
    dplyr::inner_join(
      top_gene_df |> dplyr::select(.data$gene, .data$module_id),
      by = c("gene", "module_id")
    )

  if (nrow(pair_df_top) == 0) {
    stop("No genes are available for plotting after ranking.")
  }

  label_df <- plot_df_top[0, , drop = FALSE]

  if (n_label > 0) {
    label_gene <- pair_df_top |>
      dplyr::slice_max(
        order_by = .data$rank_score,
        n = n_label,
        with_ties = FALSE
      ) |>
      dplyr::pull(.data$gene)

    label_df <- plot_df_top |>
      dplyr::filter(.data$gene %in% label_gene) |>
      dplyr::group_by(.data$gene, .data$module_id) |>
      dplyr::slice_max(
        order_by = .data$Amplitude,
        n = 1,
        with_ties = FALSE
      ) |>
      dplyr::ungroup()
  }

  circle_df <- tibble::tibble(
    r = c(0.18, 0.36, 0.54, 0.72, 0.95)
  )

  circle_path_df <- do.call(
    rbind,
    lapply(circle_df$r, function(r_use) {
      t_use <- seq(0, 2 * pi, length.out = 400)

      data.frame(
        t = t_use,
        x = r_use * sin(t_use),
        y = r_use * cos(t_use),
        r = r_use
      )
    })
  )

  angle_df <- tibble::tibble(
    angle = c(0, pi / 2, pi, 3 * pi / 2),
    label = c("0", "\u03c0/2", "\u03c0", "3\u03c0/2")
  )

  angle_df$xend <- 0.98 * sin(angle_df$angle)
  angle_df$yend <- 0.98 * cos(angle_df$angle)
  angle_df$lab_x <- 1.16 * sin(angle_df$angle)
  angle_df$lab_y <- 1.16 * cos(angle_df$angle)

  amp_breaks_raw <- pretty(c(0, amp_cap), n = 4)
  amp_breaks_raw <- amp_breaks_raw[amp_breaks_raw >= 0]
  amp_breaks_raw <- unique(amp_breaks_raw)
  amp_breaks_raw <- amp_breaks_raw[amp_breaks_raw <= amp_cap]

  amp_tick_df <- tibble::tibble(
    amp_raw = amp_breaks_raw
  )

  amp_tick_df$r <- scales::rescale(
    pmin(amp_tick_df$amp_raw, amp_cap),
    to = c(0.18, 0.95),
    from = c(0, amp_cap)
  )

  amp_tick_df$x <- amp_tick_df$r
  amp_tick_df$y <- 0
  amp_tick_df$label <- format(round(amp_tick_df$amp_raw, 2), nsmall = 2)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey86",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "grey86",
      linewidth = 0.35
    ) +
    ggplot2::geom_path(
      data = circle_path_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$r
      ),
      color = "grey82",
      linewidth = 0.40
    ) +
    ggplot2::geom_segment(
      data = angle_df,
      ggplot2::aes(
        x = 0,
        y = 0,
        xend = .data$xend,
        yend = .data$yend
      ),
      color = "grey85",
      linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      data = pair_df_top,
      ggplot2::aes(
        x = .data[[x_1]],
        y = .data[[y_1]],
        xend = .data[[x_2]],
        yend = .data[[y_2]]
      ),
      color = scales::alpha("grey55", 0.60),
      linewidth = 0.55
    ) +
    ggplot2::geom_point(
      data = plot_df_top,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        fill = .data$group
      ),
      shape = 21,
      color = "black",
      stroke = 0.26,
      alpha = 0.95,
      size = pt_size
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        label = .data$gene
      ),
      size = label_size,
      box.padding = 0.28,
      point.padding = 0.18,
      segment.color = "black",
      segment.size = 0.30,
      max.overlaps = Inf,
      seed = 123
    ) +
    ggplot2::geom_text(
      data = angle_df,
      ggplot2::aes(
        x = .data$lab_x,
        y = .data$lab_y,
        label = .data$label
      ),
      size = 4.8,
      fontface = "bold"
    ) +
    ggplot2::geom_segment(
      data = amp_tick_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y - 0.015,
        xend = .data$x,
        yend = .data$y + 0.015
      ),
      linewidth = 0.30,
      color = "grey35"
    ) +
    ggplot2::geom_text(
      data = amp_tick_df,
      ggplot2::aes(
        x = .data$x + 0.04,
        y = .data$y,
        label = .data$label
      ),
      size = 3.1,
      hjust = 0,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = 0.78,
      y = 0.06,
      label = "Amplitude",
      size = 3.8,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = 0.54,
      y = 0.84,
      label = "Peak phase",
      angle = -50,
      size = 3.5,
      fontface = "bold"
    ) +
    ggplot2::coord_equal(
      xlim = c(-1.18, 1.18),
      ylim = c(-1.15, 1.18),
      clip = "off"
    ) +
    ggplot2::scale_fill_manual(
      values = palette_use$group,
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::labs(
      title = title_use
    ) +
    theme_srb_void(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 16,
        hjust = 0.5
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )

  out <- list(
    data = plot_df_top,
    pair = pair_df_top,
    label = label_df,
    top_gene = top_gene_df,
    param_col_map = param_col_map_use,
    module_label = module_label_use,
    palette = palette_use,
    plot = p
  )

  class(out) <- c("scRhythmModulePolar", class(out))

  invisible(out)
}
