#' Plot cell-phase rho by cell type
#'
#' This function plots circular concentration rho values from a
#' \code{scRhythmPhaseRho} object returned by \code{run_phase_rho()}.
#'
#' If sample-level rho results are available, the function draws boxplots across
#' samples within each cell type. If no sample-level results are available, the
#' function draws a group-level dot plot with one rho value per cell type and
#' group. Group-level rho plots are descriptive only and do not support
#' statistical testing.
#'
#' @param rho_res A \code{scRhythmPhaseRho} object returned by
#' \code{run_phase_rho()}.
#' @param group_colors A named character vector specifying colors for groups. If
#' \code{NULL}, default scRhythmBayes group colors are used. Defaults to
#' \code{NULL}.
#' @param celltype_order A character vector specifying the order of cell types.
#' If \code{NULL}, cell types are ordered by mean rho. Defaults to
#' \code{NULL}.
#' @param show_sig A logical value indicating whether significance labels should
#' be displayed when sample-level significance results are available. Defaults
#' to \code{TRUE}.
#' @param facet_scales A character vector specifying facet scales for
#' sample-level faceted boxplots. Defaults to \code{"free_y"}.
#' @param point_size A numeric value specifying point size. Defaults to
#' \code{2.2}.
#' @param point_alpha A numeric value specifying point transparency. Defaults to
#' \code{0.85}.
#' @param title A character vector specifying the plot title. If \code{NULL}, no
#' title is shown. Defaults to \code{NULL}.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{\code{plot}}{The ggplot object.}
#'   \item{\code{data}}{The plotting data.}
#'   \item{\code{sig}}{The significance-label data, if available.}
#'   \item{\code{mode}}{The plotting mode, either \code{"sample_boxplot"} or \code{"group_dotplot"}.}
#'   \item{\code{params}}{Plotting parameters used by the function.}
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
#' rho_res <- run_phase_rho(
#'   res_phase = res_phase,
#'   metadata = srb_example_data$meta,
#'   celltype_by = "celltype",
#'   group_by = "condition",
#'   theta_by = "theta_pred",
#'   group_levels = c("A", "B"),
#'   min_cell = 10,
#'   Verbose = FALSE
#' )
#'
#' p_rho <- plot_phase_rho(rho_res)
#' p_rho$plot
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_text geom_point
#' @importFrom ggplot2 scale_fill_manual scale_color_manual
#' @importFrom ggplot2 scale_y_continuous expansion
#' @importFrom ggplot2 coord_cartesian coord_flip facet_wrap labs theme
#' @importFrom ggplot2 element_blank element_line element_text margin
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @export
plot_phase_rho <- function(
    rho_res,
    group_colors = NULL,
    celltype_order = NULL,
    show_sig = TRUE,
    facet_scales = "free_y",
    point_size = 2.2,
    point_alpha = 0.85,
    title = NULL
) {
  req_pkgs <- c(
    "ggplot2",
    "dplyr"
  )

  for (pkg in req_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for plot_phase_rho().")
    }
  }

  if (!inherits(rho_res, "scRhythmPhaseRho")) {
    stop("rho_res must be a scRhythmPhaseRho object returned by run_phase_rho().")
  }

  if (is.null(rho_res$rho_result)) {
    stop("rho_res must contain rho_result.")
  }

  if (!is.logical(show_sig) || length(show_sig) != 1) {
    stop("show_sig must be TRUE or FALSE.")
  }

  if (!is.character(facet_scales) || length(facet_scales) != 1) {
    stop("facet_scales must be a single character value.")
  }

  if (!is.numeric(point_size) || length(point_size) != 1 ||
      !is.finite(point_size) || point_size <= 0) {
    stop("point_size must be a positive numeric value.")
  }

  if (!is.numeric(point_alpha) || length(point_alpha) != 1 ||
      !is.finite(point_alpha) || point_alpha < 0 || point_alpha > 1) {
    stop("point_alpha must be a numeric value between 0 and 1.")
  }

  if (!is.null(title) && (!is.character(title) || length(title) != 1)) {
    stop("title must be NULL or a single character value.")
  }

  group_levels <- rho_res$params$group_levels

  if (is.null(group_levels) || length(group_levels) == 0) {
    group_levels <- unique(as.character(rho_res$rho_result$Group))
  }

  group_levels <- as.character(group_levels)

  if (is.null(group_colors)) {
    group_colors <- .srb_group_colors(group_levels)
  } else {
    if (is.null(names(group_colors))) {
      if (length(group_colors) < length(group_levels)) {
        stop("group_colors must contain enough colors for all groups.")
      }

      group_colors <- stats::setNames(
        group_colors[seq_along(group_levels)],
        group_levels
      )
    } else {
      miss_cols <- setdiff(group_levels, names(group_colors))

      if (length(miss_cols) > 0) {
        stop(
          "group_colors is missing colors for groups: ",
          paste(miss_cols, collapse = ", ")
        )
      }

      group_colors <- group_colors[group_levels]
    }
  }

  has_sample <- !is.null(rho_res$rho_sample) &&
    is.data.frame(rho_res$rho_sample) &&
    nrow(rho_res$rho_sample) > 0

  if (isTRUE(has_sample)) {
    plot_mode <- "sample_boxplot"
    plot_df <- as.data.frame(rho_res$rho_sample, stringsAsFactors = FALSE)

    need_cols <- c("cell_type", "Group", "rho")
    miss_cols <- setdiff(need_cols, colnames(plot_df))

    if (length(miss_cols) > 0) {
      stop(
        "rho_res$rho_sample is missing required columns: ",
        paste(miss_cols, collapse = ", ")
      )
    }

    plot_df <- plot_df[
      is.finite(plot_df$rho) &
        !is.na(plot_df$cell_type) &
        !is.na(plot_df$Group),
      ,
      drop = FALSE
    ]

    if (nrow(plot_df) == 0) {
      stop("No finite sample-level rho values are available for plotting.")
    }

    if (is.null(celltype_order)) {
      if (!is.null(rho_res$rho_sample_summary)) {
        celltype_order <- as.data.frame(
          rho_res$rho_sample_summary,
          stringsAsFactors = FALSE
        ) |>
          dplyr::group_by(.data$cell_type) |>
          dplyr::summarise(
            rho_mean_use = mean(.data$rho_mean, na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::arrange(dplyr::desc(.data$rho_mean_use)) |>
          dplyr::pull(.data$cell_type) |>
          as.character()
      } else {
        celltype_order <- plot_df |>
          dplyr::group_by(.data$cell_type) |>
          dplyr::summarise(
            rho_mean_use = mean(.data$rho, na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::arrange(dplyr::desc(.data$rho_mean_use)) |>
          dplyr::pull(.data$cell_type) |>
          as.character()
      }
    } else {
      celltype_order <- as.character(celltype_order)
    }

    plot_df$cell_type <- factor(plot_df$cell_type, levels = celltype_order)
    plot_df$Group <- factor(plot_df$Group, levels = group_levels)

    sig_df <- rho_res$sig_result

    if (!is.null(sig_df)) {
      sig_df <- as.data.frame(sig_df, stringsAsFactors = FALSE)

      if (!"cell_type" %in% colnames(sig_df)) {
        sig_df <- NULL
      } else {
        sig_df$cell_type <- factor(sig_df$cell_type, levels = celltype_order)

        if (!"x_pos" %in% colnames(sig_df)) {
          sig_df$x_pos <- mean(seq_along(group_levels))
        }

        if (!"p_label" %in% colnames(sig_df)) {
          sig_df$p_label <- NA_character_
        }

        y_df <- plot_df |>
          dplyr::group_by(.data$cell_type) |>
          dplyr::summarise(
            y_max = max(.data$rho, na.rm = TRUE),
            y_min = min(.data$rho, na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::mutate(
            y_range = pmax(.data$y_max - .data$y_min, 1e-3),
            y_pos = .data$y_max + 0.10 * .data$y_range
          )

        sig_df <- sig_df |>
          dplyr::left_join(y_df, by = "cell_type")
      }
    }

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = Group, y = rho, fill = Group)
    ) +
      ggplot2::geom_boxplot(
        width = 0.65,
        outlier.shape = NA,
        linewidth = 0.7,
        color = "black"
      ) +
      ggplot2::geom_jitter(
        width = 0.12,
        size = point_size,
        alpha = point_alpha,
        shape = 21,
        color = "black",
        stroke = 0.25
      ) +
      ggplot2::scale_fill_manual(
        values = group_colors,
        breaks = group_levels,
        guide = "none",
        drop = FALSE
      ) +
      ggplot2::facet_wrap(
        ~cell_type,
        scales = facet_scales
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0.06, 0.18))
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = title,
        x = NULL,
        y = "Rho"
      ) +
      theme_srb(base_size = 14) +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(
          color = "black",
          linewidth = 0.5
        ),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 14, face = "bold"),
        axis.title = ggplot2::element_text(size = 16, face = "bold"),
        axis.text = ggplot2::element_text(size = 12, face = "bold"),
        plot.title = ggplot2::element_text(
          size = 17,
          face = "bold",
          hjust = 0.5
        ),
        plot.margin = ggplot2::margin(8, 8, 8, 8)
      )

    if (isTRUE(show_sig) && !is.null(sig_df)) {
      p <- p +
        ggplot2::geom_text(
          data = sig_df,
          ggplot2::aes(
            x = x_pos,
            y = y_pos,
            label = p_label
          ),
          inherit.aes = FALSE,
          vjust = 0,
          size = 5,
          fontface = "bold"
        )
    }
  } else {
    plot_mode <- "group_dotplot"
    sig_df <- NULL

    plot_df <- as.data.frame(rho_res$rho_result, stringsAsFactors = FALSE)

    need_cols <- c("cell_type", "Group", "rho")
    miss_cols <- setdiff(need_cols, colnames(plot_df))

    if (length(miss_cols) > 0) {
      stop(
        "rho_res$rho_result is missing required columns: ",
        paste(miss_cols, collapse = ", ")
      )
    }

    plot_df <- plot_df[
      is.finite(plot_df$rho) &
        !is.na(plot_df$cell_type) &
        !is.na(plot_df$Group),
      ,
      drop = FALSE
    ]

    if (nrow(plot_df) == 0) {
      stop("No finite group-level rho values are available for plotting.")
    }

    if (is.null(celltype_order)) {
      if ("rho_mean" %in% colnames(plot_df)) {
        celltype_order <- plot_df |>
          dplyr::distinct(.data$cell_type, .data$rho_mean) |>
          dplyr::arrange(dplyr::desc(.data$rho_mean)) |>
          dplyr::pull(.data$cell_type) |>
          as.character()
      } else {
        celltype_order <- plot_df |>
          dplyr::group_by(.data$cell_type) |>
          dplyr::summarise(
            rho_mean_use = mean(.data$rho, na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::arrange(dplyr::desc(.data$rho_mean_use)) |>
          dplyr::pull(.data$cell_type) |>
          as.character()
      }
    } else {
      celltype_order <- as.character(celltype_order)
    }

    plot_df$cell_type <- factor(plot_df$cell_type, levels = rev(celltype_order))
    plot_df$Group <- factor(plot_df$Group, levels = group_levels)

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = cell_type, y = rho, color = Group)
    ) +
      ggplot2::geom_point(
        size = point_size + 1.2,
        alpha = point_alpha
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = sprintf("%.3f", rho),
          color = Group
        ),
        hjust = -0.20,
        size = 4.2,
        fontface = "bold",
        show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = group_colors,
        breaks = group_levels,
        name = NULL,
        drop = FALSE
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0.02, 0.18))
      ) +
      ggplot2::labs(
        title = title,
        x = NULL,
        y = "Rho"
      ) +
      theme_srb(base_size = 14) +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(
          color = "black",
          linewidth = 0.5
        ),
        axis.title = ggplot2::element_text(size = 16, face = "bold"),
        axis.text = ggplot2::element_text(size = 12, face = "bold"),
        legend.position = "top",
        legend.text = ggplot2::element_text(size = 13, face = "bold"),
        plot.title = ggplot2::element_text(
          size = 17,
          face = "bold",
          hjust = 0.5
        ),
        plot.margin = ggplot2::margin(8, 12, 8, 8)
      )
  }

  invisible(list(
    plot = p,
    data = plot_df,
    sig = sig_df,
    mode = plot_mode,
    params = list(
      group_colors = group_colors,
      group_levels = group_levels,
      celltype_order = celltype_order,
      show_sig = show_sig,
      facet_scales = facet_scales,
      point_size = point_size,
      point_alpha = point_alpha
    )
  ))
}
