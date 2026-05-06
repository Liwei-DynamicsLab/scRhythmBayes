#' Run cell-phase rho summary
#'
#' This function summarizes the circular concentration of inferred cell phases
#' within each cell type and group. If a sample-level column is provided, rho is
#' also summarized at the cell type, group, and sample level, and pairwise group
#' comparisons are performed within each cell type.
#'
#' Cells are matched between \code{res_phase} and \code{metadata} by the
#' \code{cell} column.
#'
#' @param res_phase A data.frame returned by \code{infer_cell_phase()}. It must
#' contain \code{"cell"} and the column specified by \code{theta_by}.
#' @param metadata A data.frame containing cell-level metadata.
#' @param celltype_by A character vector specifying the column name in
#' \code{metadata} used as the cell-type annotation. Defaults to
#' \code{"celltype"}.
#' @param group_by A character vector specifying the column name in
#' \code{metadata} used as the group or condition annotation. Defaults to
#' \code{"condition"}.
#' @param sample_by A character vector specifying the column name in
#' \code{metadata} used as the sample annotation. If \code{NULL}, only
#' group-level rho is computed. Defaults to \code{NULL}.
#' @param theta_by A character vector specifying the column name in
#' \code{res_phase} containing inferred cell phase in radians. Defaults to
#' \code{"theta_pred"}.
#' @param celltype_levels A character vector specifying the order of cell types.
#' If \code{NULL}, levels are sorted automatically. Defaults to \code{NULL}.
#' @param group_levels A character vector specifying the order of groups. If
#' \code{NULL}, levels are sorted automatically. Defaults to \code{NULL}.
#' @param compare_groups A character vector of length two specifying two groups
#' to compare in sample-level rho tests. If \code{NULL}, the first two group
#' levels are used. Defaults to \code{NULL}.
#' @param min_cell A numeric value specifying the minimum number of cells
#' required for each rho estimate. Defaults to \code{10}.
#' @param p_adjust_method A character vector specifying the multiple-testing
#' correction method passed to \code{\link[stats]{p.adjust}}. Defaults to
#' \code{"BH"}.
#' @param Verbose A logical value indicating whether progress messages should be
#' printed. Defaults to \code{FALSE}.
#'
#' @return An object of class \code{scRhythmPhaseRho}, containing:
#' \describe{
#'   \item{\code{data}}{Cleaned cell-level phase and metadata table.}
#'   \item{\code{rho_result}}{Group-level rho estimates for each cell type and group.}
#'   \item{\code{rho_sample}}{Sample-level rho estimates, returned only when \code{sample_by} is provided.}
#'   \item{\code{rho_sample_summary}}{Sample-level rho summary by cell type and group, returned only when \code{sample_by} is provided.}
#'   \item{\code{sig_result}}{Sample-level group-comparison results, returned only when \code{sample_by} is provided and comparison groups are available.}
#'   \item{\code{params}}{Parameters used for rho summarization.}
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
#'   group_levels = c("A", "B"),
#'   min_cell = 10,
#'   Verbose = TRUE
#' )
#'
#' rho_res
#' head(rho_res$rho_result)
#' }
#'
#' @importFrom rlang .data
#' @importFrom stats median p.adjust quantile sd wilcox.test
#' @export
run_phase_rho <- function(
    res_phase,
    metadata,
    celltype_by = "celltype",
    group_by = "condition",
    sample_by = NULL,
    theta_by = "theta_pred",
    celltype_levels = NULL,
    group_levels = NULL,
    compare_groups = NULL,
    min_cell = 10,
    p_adjust_method = "BH",
    Verbose = FALSE
) {
  req_pkgs <- c(
    "dplyr",
    "tibble",
    "circular"
  )

  for (pkg in req_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for run_phase_rho().")
    }
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

  if (!is.character(celltype_by) || length(celltype_by) != 1) {
    stop("celltype_by must be a single character value.")
  }

  if (!is.character(group_by) || length(group_by) != 1) {
    stop("group_by must be a single character value.")
  }

  if (!is.null(sample_by) &&
      (!is.character(sample_by) || length(sample_by) != 1)) {
    stop("sample_by must be NULL or a single character value.")
  }

  if (!is.null(compare_groups) &&
      (!is.character(compare_groups) || length(compare_groups) != 2)) {
    stop("compare_groups must be NULL or a character vector of length two.")
  }

  if (!is.numeric(min_cell) || length(min_cell) != 1 ||
      !is.finite(min_cell) || min_cell < 1) {
    stop("min_cell must be a positive numeric value.")
  }

  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1) {
    stop("p_adjust_method must be a single character value.")
  }

  min_cell <- as.integer(min_cell)

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

  if (!is.null(sample_by) && !sample_by %in% colnames(meta)) {
    stop("sample_by is not found in metadata.")
  }

  wrap_phase <- function(x) {
    ((as.numeric(x) %% (2 * pi)) + 2 * pi) %% (2 * pi)
  }

  rad_to_hour <- function(theta) {
    ((wrap_phase(theta) / (2 * pi)) * 24) %% 24
  }

  calc_rho <- function(theta_vec) {
    theta_vec <- as.numeric(theta_vec)
    theta_vec <- theta_vec[is.finite(theta_vec)]

    if (length(theta_vec) < min_cell) {
      return(list(
        n_cell = length(theta_vec),
        rho = NA_real_,
        mean_phase_rad = NA_real_,
        mean_phase_hour = NA_real_
      ))
    }

    theta_circular <- circular::circular(
      theta_vec,
      units = "radians",
      modulo = "2pi"
    )

    mean_phase_rad <- as.numeric(
      circular::mean.circular(theta_circular)
    )

    list(
      n_cell = length(theta_vec),
      rho = as.numeric(circular::rho.circular(theta_circular)),
      mean_phase_rad = wrap_phase(mean_phase_rad),
      mean_phase_hour = rad_to_hour(mean_phase_rad)
    )
  }

  p_to_label <- function(p_adj) {
    ifelse(
      is.na(p_adj),
      "ns",
      ifelse(
        p_adj < 1e-4,
        "****",
        ifelse(
          p_adj < 1e-3,
          "***",
          ifelse(
            p_adj < 1e-2,
            "**",
            ifelse(p_adj < 5e-2, "*", "ns")
          )
        )
      )
    )
  }

  idx <- match(res_phase$cell, meta$cell)

  if (any(is.na(idx))) {
    stop("Some cells in res_phase are not found in metadata.")
  }

  df_phase <- data.frame(
    cell = as.character(res_phase$cell),
    theta = wrap_phase(res_phase[[theta_by]]),
    phase_hour = rad_to_hour(res_phase[[theta_by]]),
    cell_type = as.character(meta[[celltype_by]][idx]),
    Group = as.character(meta[[group_by]][idx]),
    stringsAsFactors = FALSE
  )

  if (!is.null(sample_by)) {
    df_phase$Sample <- as.character(meta[[sample_by]][idx])
  }

  keep <- is.finite(df_phase$theta) &
    !is.na(df_phase$cell_type) &
    nzchar(df_phase$cell_type) &
    !is.na(df_phase$Group) &
    nzchar(df_phase$Group)

  if (!is.null(sample_by)) {
    keep <- keep &
      !is.na(df_phase$Sample) &
      nzchar(df_phase$Sample)
  }

  df_phase <- df_phase[keep, , drop = FALSE]

  if (nrow(df_phase) == 0) {
    stop("No finite phase values remain after filtering.")
  }

  if (is.null(celltype_levels)) {
    celltype_levels <- sort(unique(df_phase$cell_type))
  } else {
    celltype_levels <- as.character(celltype_levels)
  }

  if (is.null(group_levels)) {
    group_levels <- sort(unique(df_phase$Group))
  } else {
    group_levels <- as.character(group_levels)
  }

  if (length(celltype_levels) == 0) {
    stop("celltype_levels must contain at least one value.")
  }

  if (length(group_levels) == 0) {
    stop("group_levels must contain at least one value.")
  }

  df_phase$cell_type <- factor(
    df_phase$cell_type,
    levels = celltype_levels
  )

  df_phase$Group <- factor(
    df_phase$Group,
    levels = group_levels
  )

  df_phase <- df_phase[
    !is.na(df_phase$cell_type) &
      !is.na(df_phase$Group),
    ,
    drop = FALSE
  ]

  if (nrow(df_phase) == 0) {
    stop("No cells remain after applying celltype_levels and group_levels.")
  }

  if (isTRUE(Verbose)) {
    message("Use ", nrow(df_phase), " cells for phase-rho summary.")
    message(
      "Use ", length(celltype_levels), " cell type(s): ",
      paste(celltype_levels, collapse = ", ")
    )
    message(
      "Use ", length(group_levels), " group(s): ",
      paste(group_levels, collapse = ", ")
    )

    if (!is.null(sample_by)) {
      message("Use sample column: ", sample_by)
    }
  }

  rho_list <- list()
  k <- 1L

  for (ct in celltype_levels) {
    for (gp in group_levels) {
      theta_use <- df_phase$theta[
        as.character(df_phase$cell_type) == ct &
          as.character(df_phase$Group) == gp
      ]

      rho_use <- calc_rho(theta_use)

      rho_list[[k]] <- data.frame(
        cell_type = ct,
        Group = gp,
        n_cell = rho_use$n_cell,
        rho = rho_use$rho,
        mean_phase_rad = rho_use$mean_phase_rad,
        mean_phase_hour = rho_use$mean_phase_hour,
        stringsAsFactors = FALSE
      )

      k <- k + 1L
    }
  }

  rho_result <- do.call(rbind, rho_list)

  rho_result$cell_type <- factor(
    rho_result$cell_type,
    levels = celltype_levels
  )

  rho_result$Group <- factor(
    rho_result$Group,
    levels = group_levels
  )

  rho_order <- rho_result |>
    dplyr::group_by(.data$cell_type) |>
    dplyr::summarise(
      rho_mean = mean(.data$rho, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$rho_mean))

  rho_result <- rho_result |>
    dplyr::left_join(
      rho_order,
      by = "cell_type"
    )

  rho_sample <- NULL
  rho_sample_summary <- NULL
  sig_result <- NULL

  if (!is.null(sample_by)) {
    sample_list <- list()
    k <- 1L

    for (ct in celltype_levels) {
      for (gp in group_levels) {
        sample_levels_now <- sort(unique(
          as.character(df_phase$Sample[
            as.character(df_phase$cell_type) == ct &
              as.character(df_phase$Group) == gp
          ])
        ))

        sample_levels_now <- sample_levels_now[
          !is.na(sample_levels_now) & nzchar(sample_levels_now)
        ]

        if (length(sample_levels_now) == 0) {
          next
        }

        for (sm in sample_levels_now) {
          theta_use <- df_phase$theta[
            as.character(df_phase$cell_type) == ct &
              as.character(df_phase$Group) == gp &
              as.character(df_phase$Sample) == sm
          ]

          rho_use <- calc_rho(theta_use)

          sample_list[[k]] <- data.frame(
            cell_type = ct,
            Group = gp,
            Sample = sm,
            n_cell = rho_use$n_cell,
            rho = rho_use$rho,
            mean_phase_rad = rho_use$mean_phase_rad,
            mean_phase_hour = rho_use$mean_phase_hour,
            stringsAsFactors = FALSE
          )

          k <- k + 1L
        }
      }
    }

    if (length(sample_list) > 0) {
      rho_sample <- do.call(rbind, sample_list)

      rho_sample$cell_type <- factor(
        rho_sample$cell_type,
        levels = celltype_levels
      )

      rho_sample$Group <- factor(
        rho_sample$Group,
        levels = group_levels
      )

      rho_sample_summary <- rho_sample |>
        dplyr::group_by(.data$cell_type, .data$Group) |>
        dplyr::summarise(
          n_sample = sum(is.finite(.data$rho)),
          mean_n_cell = mean(.data$n_cell, na.rm = TRUE),
          rho_mean = mean(.data$rho, na.rm = TRUE),
          rho_median = stats::median(.data$rho, na.rm = TRUE),
          rho_sd = stats::sd(.data$rho, na.rm = TRUE),
          rho_se = .data$rho_sd / sqrt(pmax(.data$n_sample, 1)),
          rho_q025 = stats::quantile(.data$rho, 0.025, na.rm = TRUE),
          rho_q975 = stats::quantile(.data$rho, 0.975, na.rm = TRUE),
          .groups = "drop"
        )

      if (is.null(compare_groups)) {
        if (length(group_levels) >= 2) {
          compare_groups <- group_levels[seq_len(2)]
        }
      } else {
        compare_groups <- as.character(compare_groups)
      }

      if (!is.null(compare_groups)) {
        if (length(compare_groups) != 2) {
          stop("compare_groups must be NULL or a character vector of length two.")
        }

        if (!all(compare_groups %in% group_levels)) {
          stop("compare_groups must be present in group_levels.")
        }

        sig_list <- lapply(celltype_levels, function(ct) {
          df_now <- rho_sample[
            as.character(rho_sample$cell_type) == ct &
              as.character(rho_sample$Group) %in% compare_groups &
              is.finite(rho_sample$rho),
            ,
            drop = FALSE
          ]

          df_now$Group <- factor(
            as.character(df_now$Group),
            levels = compare_groups
          )

          n_1 <- sum(as.character(df_now$Group) == compare_groups[1])
          n_2 <- sum(as.character(df_now$Group) == compare_groups[2])

          p_val <- if (n_1 >= 2 && n_2 >= 2) {
            stats::wilcox.test(
              rho ~ Group,
              data = df_now,
              exact = FALSE
            )$p.value
          } else {
            NA_real_
          }

          data.frame(
            cell_type = ct,
            group1 = compare_groups[1],
            group2 = compare_groups[2],
            n_group1 = n_1,
            n_group2 = n_2,
            p_value = p_val,
            stringsAsFactors = FALSE
          )
        })

        sig_result <- do.call(rbind, sig_list)

        sig_result$p_adj <- stats::p.adjust(
          sig_result$p_value,
          method = p_adjust_method
        )

        sig_result$p_label <- p_to_label(sig_result$p_adj)
        sig_result$x_pos <- 1.5

        sig_result$cell_type <- factor(
          sig_result$cell_type,
          levels = celltype_levels
        )
      }
    }
  }

  out <- list(
    data = df_phase,
    rho_result = rho_result,
    rho_sample = rho_sample,
    rho_sample_summary = rho_sample_summary,
    sig_result = sig_result,
    params = list(
      celltype_by = celltype_by,
      group_by = group_by,
      sample_by = sample_by,
      theta_by = theta_by,
      celltype_levels = celltype_levels,
      group_levels = group_levels,
      compare_groups = compare_groups,
      min_cell = min_cell,
      p_adjust_method = p_adjust_method
    )
  )

  class(out) <- "scRhythmPhaseRho"

  out
}


#' Print a scRhythmPhaseRho object
#'
#' This S3 method prints a concise summary of a
#' \code{scRhythmPhaseRho} object returned by \code{run_phase_rho()}.
#'
#' @param x A \code{scRhythmPhaseRho} object.
#' @param ... Additional arguments. Currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.scRhythmPhaseRho <- function(x, ...) {
  if (!inherits(x, "scRhythmPhaseRho")) {
    stop("x must be a scRhythmPhaseRho object.")
  }

  cat("scRhythmPhaseRho object\n")

  cat("\nGroup-level rho:\n")
  print(x$rho_result)

  if (!is.null(x$rho_sample_summary)) {
    cat("\nSample-level rho summary:\n")
    print(x$rho_sample_summary)
  }

  if (!is.null(x$sig_result)) {
    cat("\nSample-level group comparison:\n")
    print(x$sig_result)
  } else {
    cat("\nSample-level group comparison: not available\n")
  }

  invisible(x)
}
