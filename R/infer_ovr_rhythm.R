
#' Infer one-vs-rest rhythmic markers
#'
#' For each target group, all remaining cells are pooled as Rest.
#' Target and Rest are fitted with infer_gene_param() using the same inputs.
#'
#' @param count Gene-by-cell raw count matrix.
#' @param metadata Cell-level metadata aligned by the \code{cell} column.
#' @param theta_pred Inferred cell phases in radians.
#' @param r_cell Cell-level rhythmic reliability values between 0 and 1.
#' @param group_by Metadata column containing the groups to compare.
#' @param ZT_by Metadata column containing zeitgeber-time bins.
#' @param target_groups Optional groups to test. By default, all groups are used.
#' @param q_cut FDR threshold for target-specific rhythmic markers.
#' @param amp_cut Minimum fitted target amplitude.
#' @param n_worker Number of parallel workers passed to \code{infer_gene_param()}.
#' @param setseed Random seed.
#' @param Verbose Logical; whether to print progress messages.
#' @param sample_by Optional metadata column identifying independent biological
#' samples for sample-aware inference.
#' @param ... Additional arguments passed to \code{infer_gene_param()}, such as
#' \code{B} and \code{sample_boot_by}.
#'
#' @return A list containing \code{ovr_result}, \code{ovr_markers},
#' \code{ovr_summary}, the name of \code{group_by}, and analysis settings.
#' The result table retains fitted Target and Rest Mesor, mean, amplitude, and
#' phase values so that model-based curves can be drawn with
#' \code{plot_ovr_curve()}.
#'
#' @export
infer_ovr_rhythm <- function(
    count,
    metadata,
    theta_pred,
    r_cell,
    group_by,
    ZT_by,
    target_groups = NULL,
    q_cut = 0.05,
    amp_cut = 0.25,
    n_worker = NULL,
    setseed = 1,
    Verbose = TRUE,
    sample_by = NULL,
    ...
) {
  if (!(is.matrix(count) || inherits(count, "Matrix"))) stop("count must be a matrix or sparse Matrix.")
  if (is.null(colnames(count))) stop("count must have cell names.")
  if (anyNA(colnames(count)) || anyDuplicated(colnames(count))) stop("Invalid or duplicated count cell IDs.")
  if (!is.data.frame(metadata)) stop("metadata must be a data.frame.")

  if (!is.character(group_by) || length(group_by) != 1L || is.na(group_by)) stop("Invalid group_by.")
  if (!is.character(ZT_by) || length(ZT_by) != 1L || is.na(ZT_by)) stop("Invalid ZT_by.")

  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!"cell" %in% names(meta)) meta$cell <- rownames(meta)
  meta$cell <- as.character(meta$cell)

  if (anyNA(meta$cell) || anyDuplicated(meta$cell) || any(!nzchar(meta$cell))) {
    stop("Invalid or duplicated metadata cell IDs.")
  }

  if (!group_by %in% names(meta)) stop("group_by is not found in metadata.")
  if (!ZT_by %in% names(meta)) stop("ZT_by is not found in metadata.")
  if (!is.null(sample_by) && !sample_by %in% names(meta)) stop("sample_by is not found in metadata.")
  if (group_by == ".ovr_group") stop("group_by cannot be '.ovr_group'.")

  id <- match(colnames(count), meta$cell)
  if (anyNA(id)) stop("Some count cells are absent from metadata.")
  meta <- meta[id, , drop = FALSE]

  if (length(theta_pred) != ncol(count) || length(r_cell) != ncol(count)) {
    stop("theta_pred and r_cell must match the number of count cells.")
  }

  if (!is.null(names(theta_pred))) {
    if (anyNA(names(theta_pred)) || anyDuplicated(names(theta_pred))) {
      stop("Invalid theta_pred cell names.")
    }
    phase_id <- match(colnames(count), names(theta_pred))
    if (anyNA(phase_id)) stop("Some count cells are absent from theta_pred.")
    theta_pred <- theta_pred[phase_id]
  }

  if (!is.null(names(r_cell))) {
    if (anyNA(names(r_cell)) || anyDuplicated(names(r_cell))) {
      stop("Invalid r_cell cell names.")
    }
    r_id <- match(colnames(count), names(r_cell))
    if (anyNA(r_id)) stop("Some count cells are absent from r_cell.")
    r_cell <- r_cell[r_id]
  }

  theta_pred <- as.numeric(theta_pred)
  r_cell <- as.numeric(r_cell)

  if (any(!is.finite(theta_pred)) || any(!is.finite(r_cell))) {
    stop("Non-finite theta_pred or r_cell values.")
  }

  if (any(r_cell < 0 | r_cell > 1)) stop("r_cell must be between 0 and 1.")
  if (anyNA(meta[[ZT_by]])) stop("ZT_by contains missing values.")
  if (anyNA(meta[[group_by]])) stop("group_by contains missing values.")

  if (!is.numeric(q_cut) || length(q_cut) != 1L ||
      !is.finite(q_cut) || q_cut <= 0 || q_cut >= 1) {
    stop("q_cut must be between 0 and 1.")
  }

  if (!is.numeric(amp_cut) || length(amp_cut) != 1L ||
      !is.finite(amp_cut) || amp_cut < 0) {
    stop("amp_cut must be non-negative.")
  }

  groups <- sort(unique(as.character(meta[[group_by]])))
  if (length(groups) < 2L) stop("OVR requires at least two groups.")
  if (any(!nzchar(groups))) stop("group_by contains empty group labels.")

  target_groups <- if (is.null(target_groups)) groups else unique(as.character(target_groups))
  if (!length(target_groups) || anyNA(target_groups)) stop("Invalid target_groups.")

  unknown <- setdiff(target_groups, groups)
  if (length(unknown)) stop("Unknown target_groups: ", paste(unknown, collapse = ", "))

  ans <- vector("list", length(target_groups))

  for (i in seq_along(target_groups)) {
    ct <- target_groups[i]
    if (Verbose) message(sprintf("[OVR %d/%d] %s vs Rest", i, length(target_groups), ct))

    meta$.ovr_group <- ifelse(as.character(meta[[group_by]]) == ct, "Target", "Rest")

    fit <- infer_gene_param(
      count = count,
      metadata = meta,
      theta_pred = theta_pred,
      r_cell = r_cell,
      strata_by = ".ovr_group",
      ZT_by = ZT_by,
      n_worker = n_worker,
      setseed = setseed,
      Verbose = FALSE,
      sample_by = sample_by,
      ...
    )

    p <- fit$param_long

    required_cols <- c(
      "gene", ".ovr_group", "Mesor", "Mean", "Amplitude", "Phase", "FDR"
    )
    if (!all(required_cols %in% names(p))) {
      stop("infer_gene_param() returned unexpected param_long columns.")
    }

    keep <- c("gene", "Mesor", "Mean", "Amplitude", "Phase", "FDR")
    tar <- p[p$.ovr_group == "Target", keep, drop = FALSE]
    rst <- p[p$.ovr_group == "Rest", keep, drop = FALSE]

    if (!nrow(tar) || !nrow(rst)) {
      stop("Target or Rest gene-level results are empty.")
    }

    if (anyDuplicated(tar$gene) || anyDuplicated(rst$gene)) {
      stop("Duplicated gene IDs in gene-level results.")
    }

    names(tar)[-1] <- paste0(names(tar)[-1], "_target")
    names(rst)[-1] <- paste0(names(rst)[-1], "_rest")

    x <- merge(tar, rst, by = "gene", all = FALSE, sort = FALSE)

    if (!nrow(x)) stop("No common genes between Target and Rest.")

    x[[group_by]] <- ct

    x$Delta_Mesor <- x$Mesor_target - x$Mesor_rest
    x$Delta_Mean <- x$Mean_target - x$Mean_rest
    x$Delta_Amp <- x$Amplitude_target - x$Amplitude_rest
    x$Delta_Phase <- ((x$Phase_target - x$Phase_rest + pi) %% (2 * pi)) - pi

    x$ovr_marker <- with(
      x,
      is.finite(FDR_target) &
        is.finite(Amplitude_target) &
        is.finite(FDR_rest) &
        FDR_target < q_cut &
        Amplitude_target >= amp_cut &
        FDR_rest >= q_cut
    )

    x$ovr_marker[is.na(x$ovr_marker)] <- FALSE

    ans[[i]] <- x[c(
      "gene", group_by,
      "Mesor_target", "Mesor_rest", "Delta_Mesor",
      "Mean_target", "Mean_rest", "Delta_Mean",
      "Amplitude_target", "Amplitude_rest", "Delta_Amp",
      "Phase_target", "Phase_rest", "Delta_Phase",
      "FDR_target", "FDR_rest", "ovr_marker"
    )]
  }

  result <- tibble::as_tibble(do.call(rbind, ans))
  markers <- result[result$ovr_marker, , drop = FALSE]

  summary <- data.frame(
    group = target_groups,
    n_marker = as.integer(table(factor(markers[[group_by]], levels = target_groups))),
    stringsAsFactors = FALSE
  )

  names(summary)[1] <- group_by
  summary <- tibble::as_tibble(summary)

  list(
    ovr_result = result,
    ovr_markers = markers,
    ovr_summary = summary,
    group_by = group_by,
    settings = list(
      ZT_by = ZT_by,
      sample_by = sample_by,
      groups = groups,
      target_groups = target_groups,
      q_cut = q_cut,
      amp_cut = amp_cut
    )
  )
}
