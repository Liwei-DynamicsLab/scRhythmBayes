#' Plot model-fitted Target-versus-Rest rhythmic curves
#'
#' Reproduces the model-fitted Target-versus-Rest curve design used for the
#' original OVR figure. Target and pooled-Rest trajectories are calculated from
#' the fitted mesor, amplitude, and phase returned by
#' \code{infer_ovr_rhythm()}. Each trajectory is divided by its mean so that the
#' horizontal reference value is one. No additional smoother is fitted.
#'
#' @param ovr A result returned by \code{infer_ovr_rhythm()}, or its
#' \code{ovr_result} table.
#' @param genes Optional character vector of genes to plot. If \code{NULL}, all
#' OVR markers for \code{target_group} are used.
#' @param target_group Single target group to plot against pooled Rest.
#' @param group_by Group column in \code{ovr_result}. It is read from an OVR
#' result list or inferred when the table contains exactly one non-result column.
#' @param n_grid Number of phase-grid points over 24 hours.
#' @param target_color Target-curve color. If \code{NULL}, the original OVR
#' cell-type palette is matched to the target group order.
#' @param rest_color Color for the pooled-Rest curve.
#' @param ncol Optional number of facet columns for multiple genes.
#'
#' @return A \code{ggplot} object. The model-derived values used for plotting
#' are stored in the \code{curve_data} attribute.
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
#'   target_groups = "1",
#'   sample_by = "sample",
#'   B = 3,
#'   n_worker = 1,
#'   Verbose = FALSE
#' )
#'
#' plot_ovr_curve(
#'   res_ovr,
#'   genes = res_ovr$ovr_markers$gene[1],
#'   target_group = "1"
#' )
#' }
#'
#' @export
plot_ovr_curve <- function(
    ovr,
    genes = NULL,
    target_group,
    group_by = NULL,
    n_grid = 400,
    target_color = NULL,
    rest_color = "#707070",
    ncol = NULL
) {
  if (!is.character(target_group) || length(target_group) != 1L ||
      is.na(target_group) || !nzchar(target_group)) {
    stop("target_group must be one non-empty character value.")
  }

  group_levels <- NULL
  if (is.list(ovr) && !is.data.frame(ovr)) {
    if (is.null(ovr$ovr_result)) stop("ovr does not contain an ovr_result table.")
    if (is.null(group_by) && is.character(ovr$group_by)) group_by <- ovr$group_by
    if (is.character(ovr$settings$groups)) group_levels <- ovr$settings$groups
    result <- ovr$ovr_result
  } else {
    result <- ovr
  }

  if (!is.data.frame(result)) stop("ovr must be an OVR result list or data.frame.")

  result_cols <- c(
    "gene", "Mesor_target", "Mesor_rest", "Delta_Mesor",
    "Mean_target", "Mean_rest", "Delta_Mean",
    "Amplitude_target", "Amplitude_rest", "Delta_Amp",
    "Phase_target", "Phase_rest", "Delta_Phase",
    "FDR_target", "FDR_rest", "ovr_marker"
  )

  if (is.null(group_by)) {
    candidates <- setdiff(names(result), result_cols)
    if (length(candidates) != 1L) {
      stop("group_by could not be inferred; supply its column name.")
    }
    group_by <- candidates
  }

  if (!is.character(group_by) || length(group_by) != 1L ||
      !group_by %in% names(result)) {
    stop("group_by must name one column in ovr_result.")
  }

  required <- c(
    "gene", group_by, "Mesor_target", "Mesor_rest",
    "Amplitude_target", "Amplitude_rest", "Phase_target", "Phase_rest",
    "ovr_marker"
  )
  missing <- setdiff(required, names(result))
  if (length(missing)) {
    stop("ovr_result is missing required columns: ", paste(missing, collapse = ", "))
  }

  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 25) {
    stop("n_grid must be a finite number of at least 25.")
  }
  n_grid <- as.integer(n_grid)

  use <- result[as.character(result[[group_by]]) == target_group, , drop = FALSE]
  if (!nrow(use)) stop("target_group is absent from ovr_result.")

  if (is.null(genes)) {
    genes <- as.character(use$gene[!is.na(use$ovr_marker) & use$ovr_marker])
    if (!length(genes)) stop("No OVR markers are available for target_group; supply genes explicitly.")
  } else {
    if (!is.character(genes) || anyNA(genes) || any(!nzchar(genes))) {
      stop("genes must be a non-empty character vector.")
    }
    genes <- unique(genes)
  }

  missing_genes <- setdiff(genes, as.character(use$gene))
  if (length(missing_genes)) {
    stop("Genes absent for target_group: ", paste(missing_genes, collapse = ", "))
  }

  use <- use[match(genes, use$gene), , drop = FALSE]
  parameter_cols <- c(
    "Mesor_target", "Mesor_rest", "Amplitude_target", "Amplitude_rest",
    "Phase_target", "Phase_rest"
  )
  if (any(!is.finite(as.matrix(use[parameter_cols])))) {
    stop("Selected genes contain non-finite fitted curve parameters.")
  }

  if (is.null(target_color)) {
    source_palette <- c(
      "#E84D33", "#55B3CF", "#455A8B", "#9ED3C9", "#13A88E",
      "#F29C7E", "#6A78A8"
    )
    if (is.null(group_levels)) group_levels <- sort(unique(as.character(result[[group_by]])))
    color_id <- match(target_group, group_levels)
    if (is.na(color_id)) color_id <- 1L
    target_color <- source_palette[(color_id - 1L) %% length(source_palette) + 1L]
  }

  if (!is.character(target_color) || length(target_color) != 1L || is.na(target_color) ||
      !is.character(rest_color) || length(rest_color) != 1L || is.na(rest_color)) {
    stop("target_color and rest_color must each be one valid color string.")
  }

  theta_hour <- seq(0, 24, length.out = n_grid)
  theta <- theta_hour / 24 * 2 * pi

  curve_data <- do.call(rbind, lapply(seq_len(nrow(use)), function(i) {
    target_absolute <- exp(
      use$Mesor_target[i] +
        use$Amplitude_target[i] * cos(theta - use$Phase_target[i])
    )
    rest_absolute <- exp(
      use$Mesor_rest[i] +
        use$Amplitude_rest[i] * cos(theta - use$Phase_rest[i])
    )

    data.frame(
      gene = rep(as.character(use$gene[i]), 2L * n_grid),
      group = rep(c("Rest", "Target"), each = n_grid),
      theta = rep(theta, 2L),
      phase_hour = rep(theta_hour, 2L),
      absolute_mean = c(rest_absolute, target_absolute),
      fitted_mean = c(
        rest_absolute / mean(rest_absolute),
        target_absolute / mean(target_absolute)
      ),
      stringsAsFactors = FALSE
    )
  }))

  curve_data$gene <- factor(curve_data$gene, levels = genes)
  curve_data$group <- factor(curve_data$group, levels = c("Rest", "Target"))

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = 1,
      linewidth = 0.35,
      colour = "#D9D9D9"
    ) +
    ggplot2::geom_line(
      data = curve_data[curve_data$group == "Rest", , drop = FALSE],
      ggplot2::aes(
        x = .data$phase_hour,
        y = .data$fitted_mean,
        group = .data$gene
      ),
      colour = rest_color,
      linewidth = 1.25,
      linetype = "22",
      lineend = "round"
    ) +
    ggplot2::geom_line(
      data = curve_data[curve_data$group == "Target", , drop = FALSE],
      ggplot2::aes(
        x = .data$phase_hour,
        y = .data$fitted_mean,
        group = .data$gene
      ),
      colour = target_color,
      linewidth = 2.1,
      lineend = "round"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 24),
      breaks = c(0, 6, 12, 18, 24),
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.12))
    ) +
    ggplot2::labs(
      title = if (length(genes) == 1L) genes else NULL,
      subtitle = target_group,
      x = "Predicted phase (h)",
      y = "Model-fitted relative expression"
    ) +
    ggplot2::theme_classic(base_size = 11, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0.5, size = 20,
        margin = ggplot2::margin(b = 1)
      ),
      plot.subtitle = ggplot2::element_text(
        face = "bold", colour = target_color, hjust = 0.5, size = 10.5,
        margin = ggplot2::margin(b = 4)
      ),
      axis.title.x = ggplot2::element_text(
        size = 10.5, face = "bold", margin = ggplot2::margin(t = 5)
      ),
      axis.title.y = ggplot2::element_text(
        size = 10.5, face = "bold", margin = ggplot2::margin(r = 5)
      ),
      axis.text = ggplot2::element_text(size = 9, colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(7, 7, 7, 7)
    )

  if (length(genes) > 1L) {
    p <- p +
      ggplot2::facet_wrap(ggplot2::vars(.data$gene), scales = "free_y", ncol = ncol) +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold", size = 12)
      )
  }

  attr(p, "curve_data") <- tibble::as_tibble(curve_data)
  p
}
