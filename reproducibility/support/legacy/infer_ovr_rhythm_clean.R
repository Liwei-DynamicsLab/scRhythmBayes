#' Infer one-vs-rest rhythmic markers
#'
#' For each target group, all remaining cells are pooled as Rest.
#' Target and Rest are fitted with infer_gene_param() using the same inputs.
#'
#' @export
infer_ovr_rhythm <- function(
    count,
    metadata,
    theta_pred,
    r_cell,
    group_by,
    target_groups = NULL,
    q_cut = 0.05,
    amp_cut = 0.25,
    n_worker = NULL,
    setseed = 1,
    Verbose = TRUE,
    sample_by = NULL,
    ...
) {
  if (!is.data.frame(metadata)) stop("metadata must be a data.frame.")
  if (!is.character(group_by) || length(group_by) != 1L || !group_by %in% names(metadata)) stop("group_by is not found in metadata.")

  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!"cell" %in% names(meta)) meta$cell <- rownames(meta)
  id <- match(colnames(count), meta$cell)
  if (anyNA(id)) stop("Some count cells are absent from metadata.")
  meta <- meta[id, , drop = FALSE]

  groups <- sort(unique(as.character(meta[[group_by]])))
  target_groups <- if (is.null(target_groups)) groups else unique(as.character(target_groups))
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
      n_worker = n_worker,
      setseed = setseed,
      Verbose = FALSE,
      sample_by = sample_by,
      ...
    )

    p <- fit$param_long

    tar <- p[p$.ovr_group == "Target", c("gene", "Amplitude", "Phase", "FDR"), drop = FALSE]
    rst <- p[p$.ovr_group == "Rest", c("gene", "Amplitude", "Phase", "FDR"), drop = FALSE]

    names(tar)[-1] <- paste0(names(tar)[-1], "_target")
    names(rst)[-1] <- paste0(names(rst)[-1], "_rest")

    x <- merge(tar, rst, by = "gene", all = FALSE, sort = FALSE)
    x[[group_by]] <- ct
    x$Delta_Amp <- x$Amplitude_target - x$Amplitude_rest
    x$Delta_Phase <- ((x$Phase_target - x$Phase_rest + pi) %% (2 * pi)) - pi
    x$ovr_marker <- with(
      x,
      is.finite(FDR_target) & is.finite(Amplitude_target) & is.finite(FDR_rest) &
      FDR_target < q_cut & Amplitude_target >= amp_cut & FDR_rest >= q_cut
    )

    ans[[i]] <- x[c(
      "gene", group_by,
      "Amplitude_target", "Amplitude_rest", "Delta_Amp",
      "Phase_target", "Phase_rest", "Delta_Phase",
      "FDR_target", "FDR_rest", "ovr_marker"
    )]
  }

  result <- tibble::as_tibble(do.call(rbind, ans))
  markers <- result[result$ovr_marker, , drop = FALSE]

  summary <- data.frame(group = target_groups, n_marker = as.integer(table(factor(markers[[group_by]], levels = target_groups))), stringsAsFactors = FALSE)
  names(summary)[1] <- group_by
  summary <- tibble::as_tibble(summary)

  list(
    ovr_result = result,
    ovr_markers = markers,
    ovr_summary = summary
  )
}
