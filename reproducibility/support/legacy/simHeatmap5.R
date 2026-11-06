suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})

CircadianPairedHeatmapByModule <- function(
    expr.matrix,
    cell.phase,
    gene.module,
    module_id,
    group.vec,
    group_levels = c("cond_1", "cond_2"),
    gene.phase.cond_1,
    gene.phase.cond_2 = NULL,
    sort_by = c("auto", "cond_1", "cond_2"),
    heatmap.pal = "RdYlBu",
    phase.style = c("bw6", "bw2"),
    title = NULL,
    use.raster = TRUE,
    raster.quality = 4,
    draw.separators = FALSE,
    separator.lwd = 1,
    zscore_rows = TRUE,
    clip_z = 2
) {
  phase.style <- match.arg(phase.style)
  sort_by <- match.arg(sort_by)

  stopifnot(!is.null(rownames(expr.matrix)))
  stopifnot(!is.null(colnames(expr.matrix)))
  stopifnot(length(cell.phase) == ncol(expr.matrix))
  stopifnot(length(group.vec) == ncol(expr.matrix))
  stopifnot(length(gene.module) == nrow(expr.matrix))
  stopifnot(!is.null(names(gene.module)))
  stopifnot(all(group_levels %in% unique(group.vec)))

  wrap_hour <- function(x) ((x %% 24) + 24) %% 24
  wrap_phase <- function(x) ((x %% (2 * pi)) + 2 * pi) %% (2 * pi)

  make_phase_annotation <- function(zt_hour, phase.style) {
    phase_out <- switch(
      phase.style,
      "bw6" = {
        phase_cat <- round(zt_hour / 4) * 4
        phase_cat <- phase_cat %% 24
        phase_cat <- factor(phase_cat, levels = c(0, 4, 8, 12, 16, 20))
        phase_cols <- c(
          "0"  = "#000000",
          "4"  = "#000000",
          "8"  = "#F5F5F5",
          "12" = "#F5F5F5",
          "16" = "#F5F5F5",
          "20" = "#000000"
        )
        list(cat = phase_cat, cols = phase_cols)
      },
      "bw2" = {
        phase_cat <- c("Night", "Day")[1 + (zt_hour >= 6 & zt_hour < 18)]
        phase_cat <- factor(phase_cat, levels = c("Night", "Day"))
        phase_cols <- c(
          "Night" = "#000000",
          "Day"   = "#F5F5F5"
        )
        list(cat = phase_cat, cols = phase_cols)
      }
    )
    phase_out
  }

  gene.module <- gene.module[rownames(expr.matrix)]
  features_use <- names(gene.module)[gene.module == module_id]
  features_use <- intersect(features_use, rownames(expr.matrix))
  stopifnot(length(features_use) > 0)

  auto_sort_map <- c(
    "0" = "cond_1",
    "1" = "cond_1",
    "2" = "cond_2",
    "3" = "cond_1",
    "4" = "cond_2",
    "5" = "cond_1",
    "6" = "cond_1",
    "7" = "cond_1"
  )

  sort_by_use <- c(
    unname(auto_sort_map[as.character(module_id)]),
    sort_by
  )[1 + (sort_by != "auto")]

  gene.phase.list <- list(
    cond_1 = gene.phase.cond_1,
    cond_2 = gene.phase.cond_2
  )

  stopifnot(!is.null(gene.phase.list[[sort_by_use]]))

  gene.phase.use <- gene.phase.list[[sort_by_use]]
  gene.phase.use <- gene.phase.use[features_use]
  gene.phase.use <- as.numeric(gene.phase.use)
  names(gene.phase.use) <- features_use

  ok_idx <- is.finite(gene.phase.use)

  gene.phase.ok <- gene.phase.use[ok_idx]
  is_gene_phase_radian <- max(gene.phase.ok, na.rm = TRUE) <= (2 * pi + 1e-8)

  gene.phase.hour <- c(
    wrap_hour(gene.phase.ok),
    wrap_hour(gene.phase.ok / (2 * pi) * 24)
  )[seq_along(gene.phase.ok) + length(gene.phase.ok) * is_gene_phase_radian]

  names(gene.phase.hour) <- names(gene.phase.ok)

  row_order <- c(
    names(sort(gene.phase.hour)),
    features_use[!ok_idx]
  )

  row_order <- unique(row_order)
  stopifnot(length(row_order) > 0)

  heat_cols <- switch(
    heatmap.pal,
    "RdYlBu"   = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    "Spectral" = rev(RColorBrewer::brewer.pal(11, "Spectral")),
    "RdBu"     = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  )

  col_fun <- circlize::colorRamp2(
    seq(-clip_z, clip_z, length.out = 11),
    heat_cols
  )

  condition_cols <- setNames(
    c("#4DBBD5", "#E64B35"),
    pair_levels
  )

  make_one_group_ht <- function(group_name, row_order, hm_name, phase_anno_name) {
    keep_cells <- which(group.vec == group_name)
    stopifnot(length(keep_cells) > 0)

    mat <- expr.matrix[row_order, keep_cells, drop = FALSE]
    phase_sub <- cell.phase[keep_cells]

    is_radian <- max(phase_sub, na.rm = TRUE) <= (2 * pi + 1e-8)
    zt_hour <- c(
      wrap_hour(as.numeric(phase_sub)),
      wrap_hour(as.numeric(phase_sub) / (2 * pi) * 24)
    )[seq_along(phase_sub) + length(phase_sub) * is_radian]

    col_order <- order(zt_hour)
    mat <- mat[, col_order, drop = FALSE]
    zt_hour <- zt_hour[col_order]

    mat <- as.matrix(mat)
    mat[!is.finite(mat)] <- 0

    z_mat <- t(scale(t(mat)))
    z_mat[!is.finite(z_mat)] <- 0
    mat <- list(raw = mat, z = z_mat)[[c("raw", "z")[1 + zscore_rows]]]

    mat[mat > clip_z] <- clip_z
    mat[mat < -clip_z] <- -clip_z

    phase_info <- make_phase_annotation(zt_hour, phase.style)
    condition_anno <- factor(rep(group_name, length(phase_info$cat)), levels = group_levels)

    anno_name_map <- list(
      left  = c("Cell phase_left", "Condition_left"),
      right = c("Cell phase", "Condition")
    )
    show_name_map <- list(
      left  = c(FALSE, FALSE),
      right = c(TRUE, TRUE)
    )
    show_legend_map <- list(
      left  = c(FALSE, FALSE),
      right = c(TRUE, TRUE)
    )

    side_key <- c("right", "left")[1 + (phase_anno_name == "Cell phase_left")]
    anno_names <- anno_name_map[[side_key]]

    anno_df <- data.frame(
      phase = phase_info$cat,
      condition = condition_anno,
      check.names = FALSE
    )
    colnames(anno_df) <- anno_names

    ha_top <- HeatmapAnnotation(
      df = anno_df,
      col = structure(
        list(phase_info$cols, condition_cols),
        names = anno_names
      ),
      show_annotation_name = show_name_map[[side_key]],
      show_legend = show_legend_map[[side_key]],
      annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
      simple_anno_size = unit(0.45, "cm"),
      gap = unit(0, "mm"),
      border = FALSE,
      gp = gpar(col = NA, lwd = 0)
    )

    ht <- Heatmap(
      mat,
      name = hm_name,
      col = col_fun,
      top_annotation = ha_top,
      show_row_names = FALSE,
      show_column_names = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      use_raster = use.raster,
      raster_quality = raster.quality,
      border = FALSE,
      rect_gp = gpar(col = NA, lwd = 0),
      heatmap_legend_param = list(
        title = "z-score",
        title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 10)
      ),
      column_title = paste0("Module ", module_id, " | ", group_name),
      column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
      show_heatmap_legend = (group_name == group_levels[2])
    )

    list(
      ht = ht,
      mat = mat,
      phase_annotation = phase_info$cat,
      phase_annotation_name = anno_names[1]
    )
  }

  left_obj <- make_one_group_ht(
    group_name = group_levels[1],
    row_order = row_order,
    hm_name = "z-score-left",
    phase_anno_name = "Cell phase_left"
  )

  right_obj <- make_one_group_ht(
    group_name = group_levels[2],
    row_order = row_order,
    hm_name = "z-score-right",
    phase_anno_name = "Cell phase"
  )

  ht_pair <- left_obj$ht + right_obj$ht

  title_use <- c(
    paste0("Module ", module_id, " | paired heatmap"),
    title
  )[1 + (!is.null(title))]

  out <- list(
    heatmap = ht_pair,
    matrix_left = left_obj$mat,
    matrix_right = right_obj$mat,
    phase_annotation_left = left_obj$phase_annotation,
    phase_annotation_right = right_obj$phase_annotation,
    phase_annotation_name_left = left_obj$phase_annotation_name,
    phase_annotation_name_right = right_obj$phase_annotation_name,
    row_order = row_order,
    module_id = module_id,
    draw_separators = draw.separators,
    separator_lwd = separator.lwd,
    title = title_use,
    group_levels = group_levels,
    sort_by_use = sort_by_use
  )

  class(out) <- "CircadianPairedHeatmap"
  return(out)
}

DrawCircadianPairedHeatmap <- function(x, merge_legend = TRUE) {
  stopifnot(inherits(x, "CircadianPairedHeatmap"))

  ht_drawn <- draw(
    x$heatmap,
    merge_legends = merge_legend,
    column_title = x$title,
    column_title_gp = grid::gpar(fontsize = 16, fontface = "bold")
  )

  phase_num_left <- as.character(x$phase_annotation_left)
  boundary_idx_left <- which(phase_num_left[-1] != phase_num_left[-length(phase_num_left)])
  x_pos_left <- boundary_idx_left / pmax(length(phase_num_left), 1)

  phase_num_right <- as.character(x$phase_annotation_right)
  boundary_idx_right <- which(phase_num_right[-1] != phase_num_right[-length(phase_num_right)])
  x_pos_right <- boundary_idx_right / pmax(length(phase_num_right), 1)

  draw_sep <- function(anno_name, xpos) {
    decorate_annotation(anno_name, {
      for (xline in xpos) {
        grid::grid.lines(
          x = unit(c(xline, xline), "npc"),
          y = unit(c(0, 1), "npc"),
          gp = gpar(col = "white", lwd = x$separator_lwd)
        )
      }
    })
  }

  sep_list <- list(
    left = list(name = x$phase_annotation_name_left, pos = x_pos_left),
    right = list(name = x$phase_annotation_name_right, pos = x_pos_right)
  )

  invisible(c(
    NULL,
    draw_sep(sep_list$left$name, sep_list$left$pos),
    draw_sep(sep_list$right$name, sep_list$right$pos)
  )[1 + 2 * isTRUE(x$draw_separators)])

  invisible(ht_drawn)
}

SaveCircadianPairedHeatmap <- function(
    x,
    pdf.file = NA_character_,
    png.file = NA_character_,
    width = 12,
    height = 7,
    res = 220
) {
  stopifnot(inherits(x, "CircadianPairedHeatmap"))

  save_pdf <- function() {
    pdf(pdf.file, width = width, height = height)
    DrawCircadianPairedHeatmap(x)
    dev.off()
    invisible(NULL)
  }

  save_png <- function() {
    png(png.file, width = width * res, height = height * res, res = res)
    DrawCircadianPairedHeatmap(x)
    dev.off()
    invisible(NULL)
  }

  c(NULL, save_pdf())[1 + (!is.na(pdf.file))]
  c(NULL, save_png())[1 + (!is.na(png.file))]

  invisible(NULL)
}
