#' Plot predicted cell phase against true time
#'
#' This function plots predicted cell phase against true cell time for simulation
#' or example data with known ground-truth cell time. Cells are matched between
#' \code{res_phase} and \code{metadata} by the \code{cell} column.
#'
#' @param res_phase A data.frame returned by \code{infer_cell_phase()}.
#' @param metadata A data.frame containing cell-level metadata.
#' @param truth_by A character vector specifying the column name in
#' \code{metadata} that contains true cell time in hours. Defaults to
#' \code{"T_true"}.
#' @param pred_by A character vector specifying the column name in
#' \code{res_phase} that contains predicted cell time in hours. Defaults to
#' \code{"ZT_pred"}.
#' @param point_color A character vector specifying the point color. Defaults to
#' \code{"#A8C7E6"}.
#' @param point_size A numeric value specifying the point size. Defaults to
#' \code{0.45}.
#' @param point_alpha A numeric value specifying the point transparency.
#' Defaults to \code{0.60}.
#' @param title A character vector specifying the plot title. Defaults to
#' \code{"Cell phase inference"}.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{\code{plot}}{A ggplot object.}
#'   \item{\code{data}}{The processed plotting data.}
#'   \item{\code{metrics}}{A data.frame containing cell number, mean absolute
#'   circular error, and median absolute circular error.}
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
#' p_scatter <- plot_phase_scatter(
#'   res_phase = res_phase,
#'   metadata = srb_example_data$meta,
#'   truth_by = "T_true"
#' )
#'
#' p_scatter$plot
#' p_scatter$metrics
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_abline geom_line geom_point
#' @importFrom ggplot2 annotate scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 coord_cartesian labs theme element_text margin
#' @importFrom stats median
#' @export
plot_phase_scatter <- function(
    res_phase,
    metadata,
    truth_by = "T_true",
    pred_by = "ZT_pred",
    point_color = "#A8C7E6",
    point_size = 0.45,
    point_alpha = 0.60,
    title = "Cell phase inference"
) {
  if (!is.data.frame(res_phase)) {
    stop("res_phase must be a data.frame.")
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame.")
  }

  if (!"cell" %in% colnames(res_phase)) {
    stop("res_phase must contain a 'cell' column.")
  }

  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!"cell" %in% colnames(meta)) {
    if (!is.null(rownames(meta))) {
      meta$cell <- rownames(meta)
    } else {
      stop("metadata must contain a 'cell' column or rownames.")
    }
  }

  if (!truth_by %in% colnames(meta)) {
    stop("truth_by is not found in metadata.")
  }

  if (!pred_by %in% colnames(res_phase)) {
    stop("pred_by is not found in res_phase.")
  }

  if (!is.numeric(point_size) || length(point_size) != 1 ||
      !is.finite(point_size) || point_size <= 0) {
    stop("point_size must be a positive numeric value.")
  }

  if (!is.numeric(point_alpha) || length(point_alpha) != 1 ||
      !is.finite(point_alpha) || point_alpha < 0 || point_alpha > 1) {
    stop("point_alpha must be a numeric value between 0 and 1.")
  }

  idx <- match(res_phase$cell, meta$cell)

  if (any(is.na(idx))) {
    stop("Some cells in res_phase are not found in metadata.")
  }

  plot_df <- data.frame(
    cell = as.character(res_phase$cell),
    true_time = as.numeric(meta[[truth_by]][idx]),
    pred_time = as.numeric(res_phase[[pred_by]]),
    stringsAsFactors = FALSE
  )

  plot_df$true_time <- plot_df$true_time %% 24
  plot_df$pred_time <- plot_df$pred_time %% 24

  plot_df <- plot_df[
    is.finite(plot_df$true_time) & is.finite(plot_df$pred_time),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) == 0) {
    stop("No finite values available for plotting.")
  }

  err <- abs(plot_df$pred_time - plot_df$true_time)
  err <- pmin(err, 24 - err)

  metrics <- data.frame(
    n_cell = nrow(plot_df),
    mean_error = mean(err),
    median_error = stats::median(err),
    stringsAsFactors = FALSE
  )

  line_df <- data.frame(
    x = seq(-1, 26, by = 0.1)
  )

  metric_lab <- paste0(
    "N = ", metrics$n_cell,
    "\nMedAE = ", sprintf("%.2f h", metrics$median_error),
    "\nMeanAE = ", sprintf("%.2f h", metrics$mean_error)
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = line_df,
      ggplot2::aes(x = x, ymin = x - 2, ymax = x + 2),
      fill = "#E9E7E7",
      alpha = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_ribbon(
      data = line_df,
      ggplot2::aes(x = x, ymin = x + 2, ymax = x + 4),
      fill = "#FCEBEC",
      alpha = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_ribbon(
      data = line_df,
      ggplot2::aes(x = x, ymin = x - 4, ymax = x - 2),
      fill = "#FCEBEC",
      alpha = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_ribbon(
      data = line_df,
      ggplot2::aes(x = x, ymin = x + 22, ymax = x + 26),
      fill = "#E9E7E7",
      alpha = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_ribbon(
      data = line_df,
      ggplot2::aes(x = x, ymin = x - 26, ymax = x - 22),
      fill = "#E9E7E7",
      alpha = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_ribbon(
      data = line_df,
      ggplot2::aes(x = x, ymin = x + 20, ymax = x + 22),
      fill = "#FCEBEC",
      alpha = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_ribbon(
      data = line_df,
      ggplot2::aes(x = x, ymin = x - 22, ymax = x - 20),
      fill = "#FCEBEC",
      alpha = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      linewidth = 0.5,
      color = "#B1C2E5"
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = x + 2),
      linetype = 2,
      linewidth = 0.5,
      color = "#DEBCD9"
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = x - 2),
      linetype = 2,
      linewidth = 0.5,
      color = "#DEBCD9"
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = x + 4),
      linetype = 3,
      linewidth = 0.5,
      color = "#B1C2E5"
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = x - 4),
      linetype = 3,
      linewidth = 0.5,
      color = "#B1C2E5"
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = x + 24),
      linetype = 2,
      linewidth = 0.5,
      color = "#DEBCD9"
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = x - 24),
      linetype = 2,
      linewidth = 0.5,
      color = "#DEBCD9"
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = x + 22),
      linetype = 2,
      linewidth = 0.5,
      color = "#DEBCD9"
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = x - 22),
      linetype = 2,
      linewidth = 0.5,
      color = "#DEBCD9"
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = x + 20),
      linetype = 3,
      linewidth = 0.5,
      color = "#B1C2E5"
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = x - 20),
      linetype = 3,
      linewidth = 0.5,
      color = "#B1C2E5"
    ) +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(x = true_time, y = pred_time),
      color = point_color,
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::annotate(
      "text",
      x = 23.6,
      y = 1.2,
      label = metric_lab,
      hjust = 1,
      vjust = 0,
      size = 3.0,
      color = "black"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 24, by = 2),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 24, by = 2),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-1, 24.5),
      ylim = c(0, 24),
      clip = "on"
    ) +
    ggplot2::labs(
      title = title,
      x = "True Time (h)",
      y = "Predicted Time (h)"
    ) +
    theme_srb(base_size = 12, aspect.ratio = 1) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )

  invisible(list(
    plot = p,
    data = plot_df,
    metrics = metrics
  ))
}


#' Plot circular cell-phase density by sampling time
#'
#' This function plots circular density summaries of inferred cell phases across
#' sampling Zeitgeber time groups. Cells are matched between \code{res_phase} and
#' \code{metadata} by the \code{cell} column.
#'
#' @param res_phase A data.frame returned by \code{infer_cell_phase()}.
#' @param metadata A data.frame containing cell-level metadata.
#' @param zt_by A character vector specifying the column name in \code{metadata}
#' that contains sampling Zeitgeber time in hours. Defaults to \code{"ZT"}.
#' @param group_by A character vector specifying the column name in
#' \code{metadata} used for faceting. If \code{NULL}, all cells are plotted in
#' one panel. Defaults to \code{"celltype"}.
#' @param theta_by A character vector specifying the column name in
#' \code{res_phase} that contains inferred cell phase in radians. Defaults to
#' \code{"theta_pred"}.
#' @param group_levels A character vector specifying the order of facet groups.
#' If \code{NULL}, groups are sorted automatically. Defaults to \code{NULL}.
#' @param zt_levels A numeric vector specifying the order of ZT groups. If
#' \code{NULL}, ZT levels are sorted automatically. Defaults to \code{NULL}.
#' @param zt_colors A named character vector specifying colors for ZT groups.
#' Names can be formatted ZT labels, such as \code{"ZT00"} and \code{"ZT04"},
#' or numeric hour labels, such as \code{"0"} and \code{"4"}. If \code{NULL},
#' the default scRhythmBayes phase palette is used. Defaults to \code{NULL}.
#' @param title A character vector specifying the plot title. If \code{NULL}, a
#' default title is used. Defaults to \code{NULL}.
#' @param Verbose A logical value indicating whether progress messages should be
#' printed. Defaults to \code{FALSE}.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{\code{plot}}{A ggplot object.}
#'   \item{\code{data}}{The processed cell-level plotting data.}
#'   \item{\code{density}}{The circular density data.}
#'   \item{\code{ring}}{The ZT ring annotation data.}
#'   \item{\code{cell_count}}{Cell-count summary by group.}
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
#' p_flower <- plot_phase_flower(
#'   res_phase = res_phase,
#'   metadata = srb_example_data$meta,
#'   zt_by = "ZT",
#'   group_by = "celltype"
#' )
#'
#' p_flower$plot
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_path geom_segment annotate
#' @importFrom ggplot2 geom_text coord_polar scale_fill_manual
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous facet_wrap
#' @importFrom ggplot2 labs theme element_text margin
#' @importFrom stats density ave
#' @export
plot_phase_flower <- function(
    res_phase,
    metadata,
    zt_by = "ZT",
    group_by = "celltype",
    theta_by = "theta_pred",
    group_levels = NULL,
    zt_levels = NULL,
    zt_colors = NULL,
    title = NULL,
    Verbose = FALSE
) {
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

  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!"cell" %in% colnames(meta)) {
    if (!is.null(rownames(meta))) {
      meta$cell <- rownames(meta)
    } else {
      stop("metadata must contain a 'cell' column or rownames.")
    }
  }

  if (!zt_by %in% colnames(meta)) {
    stop("zt_by is not found in metadata.")
  }

  if (!is.null(group_by) && !group_by %in% colnames(meta)) {
    stop("group_by is not found in metadata.")
  }

  bw <- 0.35
  n_grid <- 720

  ring_inner <- 0.98
  ring_outer <- 1.18
  density_inner <- 1.10
  density_scale <- 0.95

  wrap_phase <- function(x) {
    ((as.numeric(x) %% (2 * pi)) + 2 * pi) %% (2 * pi)
  }

  hour_to_rad <- function(h) {
    (as.numeric(h) %% 24) / 24 * 2 * pi
  }

  idx <- match(res_phase$cell, meta$cell)

  if (any(is.na(idx))) {
    stop("Some cells in res_phase are not found in metadata.")
  }

  group_vec <- if (is.null(group_by)) {
    rep("All", nrow(res_phase))
  } else {
    as.character(meta[[group_by]][idx])
  }

  plot_df <- data.frame(
    cell = as.character(res_phase$cell),
    theta = wrap_phase(res_phase[[theta_by]]),
    ZT = as.numeric(meta[[zt_by]][idx]),
    Group = group_vec,
    stringsAsFactors = FALSE
  )

  plot_df <- plot_df[
    is.finite(plot_df$theta) &
      is.finite(plot_df$ZT) &
      !is.na(plot_df$Group),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) == 0) {
    stop("No finite values available for plotting.")
  }

  if (is.null(zt_levels)) {
    zt_levels <- sort(unique(plot_df$ZT))
  } else {
    zt_levels <- as.numeric(zt_levels)
  }

  zt_levels <- zt_levels[is.finite(zt_levels)]

  if (length(zt_levels) == 0) {
    stop("zt_levels must contain at least one finite value.")
  }

  if (is.null(group_levels)) {
    group_levels <- sort(unique(plot_df$Group))
  } else {
    group_levels <- as.character(group_levels)
  }

  if (length(group_levels) == 0) {
    stop("group_levels must contain at least one group.")
  }

  zt_keys <- .srb_zt_key(zt_levels)
  zt_labels <- .srb_zt_label(zt_levels)
  zt_colors <- .srb_match_zt_colors(zt_levels, zt_colors)

  plot_df$ZT_group <- factor(
    .srb_zt_key(plot_df$ZT),
    levels = zt_keys
  )

  plot_df$ZT_label <- factor(
    .srb_zt_label(plot_df$ZT),
    levels = zt_labels
  )

  plot_df$Group <- factor(
    plot_df$Group,
    levels = group_levels
  )

  plot_df <- plot_df[
    !is.na(plot_df$ZT_group) & !is.na(plot_df$Group),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) == 0) {
    stop("No cells remain after applying zt_levels and group_levels.")
  }

  if (isTRUE(Verbose)) {
    message("Use ", nrow(plot_df), " cells for circular phase-density plotting.")
    message("Use ", length(group_levels), " group(s): ", paste(group_levels, collapse = ", "))
    message("Use ", length(zt_keys), " ZT group(s): ", paste(zt_labels, collapse = ", "))
  }

  circ_density_df <- function(theta_vec, zt_key, group_name) {
    theta_vec <- as.numeric(theta_vec)
    theta_vec <- theta_vec[is.finite(theta_vec)]
    theta_vec <- wrap_phase(theta_vec)

    if (length(theta_vec) < 2) {
      return(data.frame(
        theta = numeric(0),
        density = numeric(0),
        ZT_group = character(0),
        Group = character(0),
        stringsAsFactors = FALSE
      ))
    }

    theta_ext <- c(theta_vec - 2 * pi, theta_vec, theta_vec + 2 * pi)

    dens <- stats::density(
      theta_ext,
      from = 0,
      to = 2 * pi,
      n = n_grid,
      bw = bw
    )

    data.frame(
      theta = dens$x,
      density = dens$y,
      ZT_group = zt_key,
      Group = group_name,
      stringsAsFactors = FALSE
    )
  }

  density_list <- list()
  k <- 1L

  for (g in group_levels) {
    for (zt_key in zt_keys) {
      theta_sub <- plot_df$theta[
        as.character(plot_df$Group) == g &
          as.character(plot_df$ZT_group) == zt_key
      ]

      density_list[[k]] <- circ_density_df(
        theta_vec = theta_sub,
        zt_key = zt_key,
        group_name = g
      )

      k <- k + 1L
    }
  }

  density_df <- do.call(rbind, density_list)

  if (is.null(density_df) || nrow(density_df) == 0) {
    stop("Not enough cells per group/ZT to estimate circular density.")
  }

  density_df$density_scaled <- stats::ave(
    density_df$density,
    density_df$Group,
    density_df$ZT_group,
    FUN = function(x) {
      mx <- max(x, na.rm = TRUE)

      if (!is.finite(mx) || mx <= 0) {
        rep(0, length(x))
      } else {
        x / mx
      }
    }
  )

  density_df$r_inner <- density_inner
  density_df$r_outer <- density_inner + density_df$density_scaled * density_scale
  density_df$Group <- factor(density_df$Group, levels = group_levels)
  density_df$ZT_group <- factor(density_df$ZT_group, levels = zt_keys)

  make_ring_df <- function(zt_values, group_levels) {
    zt_values <- sort(unique(as.numeric(zt_values)))
    zt_values <- zt_values[is.finite(zt_values)]

    if (length(zt_values) == 1) {
      lefts <- zt_values - 12
      rights <- zt_values + 12
    } else {
      zt_ext <- c(
        zt_values[length(zt_values)] - 24,
        zt_values,
        zt_values[1] + 24
      )

      mids <- (zt_ext[-length(zt_ext)] + zt_ext[-1]) / 2
      lefts <- mids[-length(mids)]
      rights <- mids[-1]
    }

    ring_list <- vector("list", length(zt_values))

    for (i in seq_along(zt_values)) {
      hour_raw <- seq(lefts[i], rights[i], length.out = 240)
      hour <- hour_raw %% 24

      ring_list[[i]] <- data.frame(
        hour = hour,
        theta = hour_to_rad(hour),
        ZT_group = .srb_zt_key(zt_values[i]),
        r0 = ring_inner,
        r1 = ring_outer,
        stringsAsFactors = FALSE
      )
    }

    ring_one <- do.call(rbind, ring_list)

    out <- do.call(rbind, lapply(group_levels, function(g) {
      tmp <- ring_one
      tmp$Group <- g
      tmp
    }))

    out$Group <- factor(out$Group, levels = group_levels)
    out$ZT_group <- factor(out$ZT_group, levels = zt_keys)

    out
  }

  ring_df <- make_ring_df(
    zt_values = zt_levels,
    group_levels = group_levels
  )

  circle_df <- data.frame(
    theta = seq(0, 2 * pi, length.out = 500),
    r = 2.05
  )

  major_tick_hours <- seq(0, 20, by = 4)
  minor_tick_hours <- seq(2, 22, by = 4)

  major_tick_df <- data.frame(
    hour = major_tick_hours,
    theta = hour_to_rad(major_tick_hours)
  )

  minor_tick_df <- data.frame(
    hour = minor_tick_hours,
    theta = hour_to_rad(minor_tick_hours)
  )

  label_df_major <- data.frame(
    hour = major_tick_hours,
    theta = hour_to_rad(major_tick_hours),
    r = 2.18,
    lab = paste0(major_tick_hours, "h")
  )

  label_df_minor <- data.frame(
    hour = minor_tick_hours,
    theta = hour_to_rad(minor_tick_hours),
    r = 2.11,
    lab = paste0(minor_tick_hours, "h")
  )

  n_df <- as.data.frame(table(plot_df$Group), stringsAsFactors = FALSE)
  colnames(n_df) <- c("Group", "n_cells")
  n_df$Group <- factor(n_df$Group, levels = group_levels)
  n_df$lab <- paste0("n = ", n_df$n_cells)

  if (is.null(title)) {
    title <- "ZT Round Distribution of Predicted Phase"
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = density_df,
      ggplot2::aes(
        x = theta,
        ymin = r_inner,
        ymax = r_outer,
        group = interaction(Group, ZT_group),
        fill = ZT_group
      ),
      alpha = 0.35,
      color = NA
    ) +
    ggplot2::geom_ribbon(
      data = ring_df,
      ggplot2::aes(
        x = theta,
        ymin = r0,
        ymax = r1,
        fill = ZT_group
      ),
      alpha = 0.95,
      color = NA
    ) +
    ggplot2::geom_path(
      data = circle_df,
      ggplot2::aes(x = theta, y = r),
      color = "grey45",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = major_tick_df,
      ggplot2::aes(
        x = theta,
        xend = theta,
        y = 0,
        yend = 2.05
      ),
      color = "grey75",
      linewidth = 0.45
    ) +
    ggplot2::geom_segment(
      data = minor_tick_df,
      ggplot2::aes(
        x = theta,
        xend = theta,
        y = 1.93,
        yend = 2.05
      ),
      color = "grey45",
      linewidth = 0.45
    ) +
    ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = 0,
      ymax = 0.92,
      fill = "white",
      color = NA
    ) +
    ggplot2::geom_text(
      data = label_df_major,
      ggplot2::aes(x = theta, y = r, label = lab),
      size = 4.8,
      fontface = "bold",
      color = "grey20"
    ) +
    ggplot2::geom_text(
      data = label_df_minor,
      ggplot2::aes(x = theta, y = r, label = lab),
      size = 3.6,
      fontface = "plain",
      color = "grey35"
    ) +
    ggplot2::coord_polar(start = 0, direction = 1, clip = "off") +
    ggplot2::scale_fill_manual(
      values = zt_colors,
      breaks = zt_keys,
      labels = zt_labels,
      name = "ZT group"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 2 * pi)) +
    ggplot2::scale_y_continuous(limits = c(0, 2.20)) +
    ggplot2::facet_wrap(~Group, nrow = 1) +
    theme_srb_void(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      strip.text = ggplot2::element_text(size = 16, face = "bold"),
      plot.title = ggplot2::element_text(
        size = 20,
        face = "bold",
        hjust = 0.5
      ),
      plot.margin = ggplot2::margin(20, 30, 20, 20)
    ) +
    ggplot2::labs(
      title = title
    )

  invisible(list(
    plot = p,
    data = plot_df,
    density = density_df,
    ring = ring_df,
    cell_count = n_df,
    params = list(
      zt_by = zt_by,
      group_by = group_by,
      theta_by = theta_by,
      group_levels = group_levels,
      zt_levels = zt_levels,
      zt_keys = zt_keys,
      zt_labels = zt_labels,
      zt_colors = zt_colors
    )
  ))
}


#' Plot cell phases on a cosine-sine circle
#'
#' This function plots inferred cell phases on the unit circle using cosine and
#' sine coordinates. Cells are matched between \code{res_phase} and
#' \code{metadata} by the \code{cell} column.
#'
#' @param res_phase A data.frame returned by \code{infer_cell_phase()}.
#' @param metadata A data.frame containing cell-level metadata.
#' @param zt_by A character vector specifying the column name in \code{metadata}
#' that contains sampling Zeitgeber time in hours. Defaults to \code{"ZT"}.
#' @param group_by A character vector specifying the column name in
#' \code{metadata} used for faceting. If \code{NULL}, all cells are plotted in
#' one panel. Defaults to \code{"celltype"}.
#' @param theta_by A character vector specifying the column name in
#' \code{res_phase} that contains inferred cell phase in radians. Defaults to
#' \code{"theta_pred"}.
#' @param group_levels A character vector specifying the order of facet groups.
#' If \code{NULL}, groups are sorted automatically. Defaults to \code{NULL}.
#' @param zt_levels A numeric vector specifying the order of ZT groups. If
#' \code{NULL}, ZT levels are sorted automatically. Defaults to \code{NULL}.
#' @param zt_colors A named character vector specifying colors for ZT groups.
#' Names can be formatted ZT labels, such as \code{"ZT00"} and \code{"ZT04"},
#' or numeric hour labels, such as \code{"0"} and \code{"4"}. If \code{NULL},
#' the default scRhythmBayes phase palette is used. Defaults to \code{NULL}.
#' @param point_size A numeric value specifying the point size. Defaults to
#' \code{3.5}.
#' @param point_alpha A numeric value specifying the point transparency.
#' Defaults to \code{0.40}.
#' @param title A character vector specifying the plot title. If \code{NULL}, no
#' title is shown. Defaults to \code{NULL}.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{\code{plot}}{A ggplot object.}
#'   \item{\code{data}}{The processed plotting data.}
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
#' p_circle <- plot_phase_circle(
#'   res_phase = res_phase,
#'   metadata = srb_example_data$meta,
#'   zt_by = "ZT",
#'   group_by = "celltype"
#' )
#'
#' p_circle$plot
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_path geom_hline geom_vline geom_point
#' @importFrom ggplot2 geom_text annotate coord_equal scale_color_manual
#' @importFrom ggplot2 facet_wrap labs theme element_text margin arrow
#' @importFrom grid unit
#' @export
plot_phase_circle <- function(
    res_phase,
    metadata,
    zt_by = "ZT",
    group_by = "celltype",
    theta_by = "theta_pred",
    group_levels = NULL,
    zt_levels = NULL,
    zt_colors = NULL,
    point_size = 3.5,
    point_alpha = 0.40,
    title = NULL
) {
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

  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!"cell" %in% colnames(meta)) {
    if (!is.null(rownames(meta))) {
      meta$cell <- rownames(meta)
    } else {
      stop("metadata must contain a 'cell' column or rownames.")
    }
  }

  if (!zt_by %in% colnames(meta)) {
    stop("zt_by is not found in metadata.")
  }

  if (!is.null(group_by) && !group_by %in% colnames(meta)) {
    stop("group_by is not found in metadata.")
  }

  if (!is.numeric(point_size) || length(point_size) != 1 ||
      !is.finite(point_size) || point_size <= 0) {
    stop("point_size must be a positive numeric value.")
  }

  if (!is.numeric(point_alpha) || length(point_alpha) != 1 ||
      !is.finite(point_alpha) || point_alpha < 0 || point_alpha > 1) {
    stop("point_alpha must be a numeric value between 0 and 1.")
  }

  wrap_phase <- function(x) {
    ((as.numeric(x) %% (2 * pi)) + 2 * pi) %% (2 * pi)
  }

  idx <- match(res_phase$cell, meta$cell)

  if (any(is.na(idx))) {
    stop("Some cells in res_phase are not found in metadata.")
  }

  group_vec <- if (is.null(group_by)) {
    rep("All", nrow(res_phase))
  } else {
    as.character(meta[[group_by]][idx])
  }

  plot_df <- data.frame(
    cell = as.character(res_phase$cell),
    theta = wrap_phase(res_phase[[theta_by]]),
    ZT = as.numeric(meta[[zt_by]][idx]),
    Group = group_vec,
    stringsAsFactors = FALSE
  )

  plot_df$x <- cos(plot_df$theta)
  plot_df$y <- sin(plot_df$theta)

  plot_df <- plot_df[
    is.finite(plot_df$theta) &
      is.finite(plot_df$ZT) &
      is.finite(plot_df$x) &
      is.finite(plot_df$y) &
      !is.na(plot_df$Group),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) == 0) {
    stop("No finite values available for plotting.")
  }

  if (is.null(zt_levels)) {
    zt_levels <- sort(unique(plot_df$ZT))
  } else {
    zt_levels <- as.numeric(zt_levels)
  }

  zt_levels <- zt_levels[is.finite(zt_levels)]

  if (length(zt_levels) == 0) {
    stop("zt_levels must contain at least one finite value.")
  }

  if (is.null(group_levels)) {
    group_levels <- sort(unique(plot_df$Group))
  } else {
    group_levels <- as.character(group_levels)
  }

  if (length(group_levels) == 0) {
    stop("group_levels must contain at least one group.")
  }

  zt_keys <- .srb_zt_key(zt_levels)
  zt_labels <- .srb_zt_label(zt_levels)
  zt_colors <- .srb_match_zt_colors(zt_levels, zt_colors)

  plot_df$ZT_group <- factor(
    .srb_zt_key(plot_df$ZT),
    levels = zt_keys
  )

  plot_df$ZT_label <- factor(
    .srb_zt_label(plot_df$ZT),
    levels = zt_labels
  )

  plot_df$Group <- factor(
    plot_df$Group,
    levels = group_levels
  )

  plot_df <- plot_df[
    !is.na(plot_df$ZT_group) & !is.na(plot_df$Group),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) == 0) {
    stop("No cells remain after applying zt_levels and group_levels.")
  }

  circle_df <- data.frame(
    t = seq(0, 2 * pi, length.out = 500)
  )
  circle_df$x <- cos(circle_df$t)
  circle_df$y <- sin(circle_df$t)

  axis_df <- data.frame(
    x = c(1.14, 0, -1.14, 0),
    y = c(0, 1.14, 0, -1.14),
    lab = c("0", "pi/2", "pi", "3pi/2"),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_path(
      data = circle_df,
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      color = "grey45",
      linewidth = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey78",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "grey78",
      linewidth = 0.35
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = ZT_group),
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::scale_color_manual(
      values = zt_colors,
      breaks = zt_keys,
      labels = zt_labels,
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::geom_text(
      data = axis_df,
      ggplot2::aes(x = x, y = y, label = lab),
      inherit.aes = FALSE,
      size = 5.0,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "segment",
      x = -1.00,
      y = -1.12,
      xend = -0.60,
      yend = -1.12,
      linewidth = 0.55,
      arrow = ggplot2::arrow(length = grid::unit(0.12, "inches"))
    ) +
    ggplot2::annotate(
      "segment",
      x = -1.00,
      y = -1.12,
      xend = -1.00,
      yend = -0.72,
      linewidth = 0.55,
      arrow = ggplot2::arrow(length = grid::unit(0.12, "inches"))
    ) +
    ggplot2::annotate(
      "text",
      x = -0.56,
      y = -1.24,
      label = "Cos(phi)",
      size = 5.4,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = -1.16,
      y = -0.66,
      label = "Sin(phi)",
      angle = 90,
      size = 5.4,
      fontface = "bold"
    ) +
    ggplot2::coord_equal(
      xlim = c(-1.28, 1.18),
      ylim = c(-1.28, 1.18),
      clip = "off"
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    theme_srb_void(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 15
      ),
      plot.margin = ggplot2::margin(18, 30, 18, 22)
    )

  if (!is.null(group_by)) {
    p <- p +
      ggplot2::facet_wrap(~Group) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(
          size = 13,
          face = "bold"
        )
      )
  }

  invisible(list(
    plot = p,
    data = plot_df,
    params = list(
      zt_by = zt_by,
      group_by = group_by,
      theta_by = theta_by,
      group_levels = group_levels,
      zt_levels = zt_levels,
      zt_keys = zt_keys,
      zt_labels = zt_labels,
      zt_colors = zt_colors
    )
  ))
}

#' Plot cell phase rose histogram
#'
#' This function plots polar histograms of inferred cell phase. Cells are
#' matched between \code{res_phase} and \code{metadata} by the \code{cell}
#' column.
#'
#' @param res_phase A data.frame returned by \code{infer_cell_phase()}.
#' @param metadata A data.frame containing cell-level metadata.
#' @param group_by A character vector specifying the column name in
#' \code{metadata} used for grouping. If \code{NULL}, all cells are plotted
#' together. Defaults to \code{NULL}.
#' @param split_by A character vector specifying the column name in
#' \code{metadata} used for faceting. If \code{NULL}, no facet is used.
#' Defaults to \code{NULL}.
#' @param theta_by A character vector specifying the column name in
#' \code{res_phase} that contains inferred cell phase in radians. Defaults to
#' \code{"theta_pred"}.
#' @param n_bin A numeric value specifying the number of circular bins. Defaults
#' to \code{24}.
#' @param group_colors A named character vector specifying colors for groups. If
#' \code{NULL}, the default scRhythmBayes group palette is used. Defaults to
#' \code{NULL}.
#' @param show_count A logical value indicating whether bin counts should be
#' displayed. Defaults to \code{TRUE}.
#' @param same_radius A logical value indicating whether all panels should use
#' the same radial scale. Defaults to \code{TRUE}.
#' @param title A character vector specifying the plot title. If \code{NULL}, no
#' title is shown. Defaults to \code{NULL}.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{\code{plot}}{A ggplot object.}
#'   \item{\code{data}}{The processed cell-level plotting data.}
#'   \item{\code{histogram}}{The phase-bin histogram data.}
#'   \item{\code{cell_count}}{Cell-count summary by group and facet.}
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
#' p_rose <- plot_phase_rose(
#'   res_phase = res_phase,
#'   metadata = srb_example_data$meta,
#'   group_by = "condition",
#'   split_by = "celltype"
#' )
#'
#' p_rose$plot
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text coord_polar
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_fill_manual
#' @importFrom ggplot2 facet_grid facet_wrap vars labs theme element_blank
#' @importFrom ggplot2 element_text element_line margin
#' @importFrom graphics hist
#' @importFrom stats setNames ave
#' @export
plot_phase_rose <- function(
    res_phase,
    metadata,
    group_by = NULL,
    split_by = NULL,
    theta_by = "theta_pred",
    n_bin = 24,
    group_colors = NULL,
    show_count = TRUE,
    same_radius = TRUE,
    title = NULL
) {
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

  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!"cell" %in% colnames(meta)) {
    if (!is.null(rownames(meta))) {
      meta$cell <- rownames(meta)
    } else {
      stop("metadata must contain a 'cell' column or rownames.")
    }
  }

  if (!is.null(group_by) && !group_by %in% colnames(meta)) {
    stop("group_by is not found in metadata.")
  }

  if (!is.null(split_by) && !split_by %in% colnames(meta)) {
    stop("split_by is not found in metadata.")
  }

  if (!is.numeric(n_bin) || length(n_bin) != 1 ||
      !is.finite(n_bin) || n_bin < 4) {
    stop("n_bin must be a numeric value greater than or equal to 4.")
  }

  n_bin <- as.integer(n_bin)

  if (!is.logical(show_count) || length(show_count) != 1) {
    stop("show_count must be TRUE or FALSE.")
  }

  if (!is.logical(same_radius) || length(same_radius) != 1) {
    stop("same_radius must be TRUE or FALSE.")
  }

  wrap_phase <- function(x) {
    ((as.numeric(x) %% (2 * pi)) + 2 * pi) %% (2 * pi)
  }

  idx <- match(res_phase$cell, meta$cell)

  if (any(is.na(idx))) {
    stop("Some cells in res_phase are not found in metadata.")
  }

  group_vec <- if (is.null(group_by)) {
    rep("All", nrow(res_phase))
  } else {
    as.character(meta[[group_by]][idx])
  }

  split_vec <- if (is.null(split_by)) {
    rep("All", nrow(res_phase))
  } else {
    as.character(meta[[split_by]][idx])
  }

  plot_df <- data.frame(
    cell = as.character(res_phase$cell),
    theta = wrap_phase(res_phase[[theta_by]]),
    Group = group_vec,
    Split = split_vec,
    stringsAsFactors = FALSE
  )

  plot_df <- plot_df[
    is.finite(plot_df$theta) &
      !is.na(plot_df$Group) &
      !is.na(plot_df$Split),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) == 0) {
    stop("No finite values available for plotting.")
  }

  group_levels <- sort(unique(plot_df$Group))
  split_levels <- sort(unique(plot_df$Split))

  plot_df$Group <- factor(plot_df$Group, levels = group_levels)
  plot_df$Split <- factor(plot_df$Split, levels = split_levels)

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

  breaks_theta <- seq(0, 2 * pi, length.out = n_bin + 1)
  mids_theta <- (breaks_theta[-1] + breaks_theta[-length(breaks_theta)]) / 2
  width_theta <- 2 * pi / n_bin

  hist_one <- function(df_one) {
    h <- graphics::hist(
      df_one$theta,
      breaks = breaks_theta,
      plot = FALSE,
      include.lowest = TRUE,
      right = FALSE
    )

    data.frame(
      theta_mid = mids_theta,
      count = as.numeric(h$counts),
      prop = as.numeric(h$counts) / sum(h$counts),
      stringsAsFactors = FALSE
    )
  }

  hist_list <- list()
  k <- 1L

  for (sp in split_levels) {
    for (gp in group_levels) {
      df_one <- plot_df[
        as.character(plot_df$Split) == sp &
          as.character(plot_df$Group) == gp,
        ,
        drop = FALSE
      ]

      if (nrow(df_one) == 0) {
        next
      }

      tmp <- hist_one(df_one)
      tmp$Split <- sp
      tmp$Group <- gp
      tmp$n_cell <- nrow(df_one)

      hist_list[[k]] <- tmp
      k <- k + 1L
    }
  }

  hist_df <- do.call(rbind, hist_list)

  if (is.null(hist_df) || nrow(hist_df) == 0) {
    stop("No histogram data available for plotting.")
  }

  hist_df$Split <- factor(hist_df$Split, levels = split_levels)
  hist_df$Group <- factor(hist_df$Group, levels = group_levels)
  hist_df$label <- ifelse(hist_df$count > 0, as.character(hist_df$count), "")

  if (isTRUE(same_radius)) {
    ymax_use <- max(hist_df$prop, na.rm = TRUE)
    hist_df$ymax <- ymax_use
  } else {
    hist_df$ymax <- stats::ave(
      hist_df$prop,
      hist_df$Split,
      hist_df$Group,
      FUN = function(x) {
        max(x, na.rm = TRUE)
      }
    )
  }

  hist_df$ymax <- pmax(hist_df$ymax, 1e-6)

  n_df <- unique(hist_df[, c("Split", "Group", "n_cell", "ymax")])
  n_df$x <- 0.35
  n_df$y <- n_df$ymax * 0.88
  n_df$label <- paste0("N = ", n_df$n_cell)

  y_max_global <- max(hist_df$ymax, na.rm = TRUE)

  p <- ggplot2::ggplot(
    hist_df,
    ggplot2::aes(x = theta_mid, y = prop, fill = Group)
  ) +
    ggplot2::geom_col(
      width = width_theta * 0.94,
      alpha = 0.90,
      colour = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      data = n_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 4.8,
      hjust = 0,
      vjust = 1,
      fontface = "bold"
    ) +
    ggplot2::coord_polar(
      start = 0,
      direction = 1,
      clip = "off"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 2 * pi),
      breaks = seq(0, 2 * pi - 2 * pi / 8, length.out = 8),
      labels = c(
        "0:00", "3:00", "6:00", "9:00",
        "12:00", "15:00", "18:00", "21:00"
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-y_max_global * 0.22, y_max_global * 1.18),
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_manual(
      values = group_colors,
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    theme_srb_void(base_size = 15) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 13, face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        linewidth = 0.35,
        colour = "grey80"
      ),
      panel.grid.major.y = ggplot2::element_line(
        linewidth = 0.35,
        colour = "grey88"
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 17
      ),
      strip.text = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      legend.position = "right",
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )

  if (isTRUE(show_count)) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(
          y = prop + ymax * 0.035,
          label = label
        ),
        size = 3.0,
        fontface = "bold",
        show.legend = FALSE
      )
  }

  if (!is.null(group_by) && !is.null(split_by)) {
    p <- p +
      ggplot2::facet_grid(
        rows = ggplot2::vars(Split),
        cols = ggplot2::vars(Group)
      )
  } else if (!is.null(group_by)) {
    p <- p +
      ggplot2::facet_wrap(~Group)
  } else if (!is.null(split_by)) {
    p <- p +
      ggplot2::facet_wrap(~Split)
  }

  invisible(list(
    plot = p,
    data = plot_df,
    histogram = hist_df,
    cell_count = n_df,
    params = list(
      group_by = group_by,
      split_by = split_by,
      theta_by = theta_by,
      n_bin = n_bin,
      group_levels = group_levels,
      split_levels = split_levels,
      group_colors = group_colors,
      show_count = show_count,
      same_radius = same_radius
    )
  ))
}
