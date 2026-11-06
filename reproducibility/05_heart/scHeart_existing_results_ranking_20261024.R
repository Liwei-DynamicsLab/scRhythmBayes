heart_root <- "/path/to/analysis_home/scRhythmBayes/data/scHeart_ACR_NCR_new_20261022"
`%>%` <- dplyr::`%>%`
heart_figure_dir <- file.path(
  heart_root,
  "figures_complete_noif_20261023"
)
heart_rho_path <- file.path(
  heart_figure_dir,
  "heart_run_phase_rho.rds"
)
heart_dtw_path <- file.path(
  heart_figure_dir,
  "heart_run_dtw_similarity.rds"
)
figure6_root <- "/path/to/source_archive/new figures /figure6"
review_dir <- file.path(
  figure6_root,
  "existing_results_review_20261024"
)
ranking_dir <- file.path(review_dir, "ranking")
qa_dir <- file.path(review_dir, "qa")

dir.create(ranking_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(heart_rho_path), file.exists(heart_dtw_path))

save_vector_output <- function(plot, prefix, width, height) {
  grDevices::cairo_pdf(
    paste0(prefix, ".pdf"),
    width = width,
    height = height,
    family = "Arial"
  )
  print(plot)
  grDevices::dev.off()
  svglite::svglite(
    paste0(prefix, ".svg"),
    width = width,
    height = height,
    system_fonts = list(sans = "Arial")
  )
  print(plot)
  grDevices::dev.off()
}

heart_rho <- readRDS(heart_rho_path)
heart_dtw <- readRDS(heart_dtw_path)

stopifnot(
  all(c("rho_sample_summary", "sig_result", "params") %in%
    names(heart_rho)),
  all(c("distance", "annotation", "params") %in%
    names(heart_dtw)),
  identical(heart_rho$params$compare_groups, c("NCR", "ACR")),
  identical(heart_rho$params$min_cell, 10L),
  identical(heart_dtw$params$group_levels, c("NCR", "ACR")),
  identical(heart_dtw$params$min_cell, 10L),
  identical(heart_dtw$params$normalize_distance, TRUE)
)

heart_rho_summary <- heart_rho$rho_sample_summary %>%
  tibble::as_tibble() %>%
  dplyr::select(
    cell_type, Group, n_sample, mean_n_cell,
    rho_mean, rho_median
  ) %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = c(
      n_sample, mean_n_cell, rho_mean, rho_median
    ),
    names_glue = "{.value}_{Group}"
  )

heart_rho_tests <- heart_rho$sig_result %>%
  tibble::as_tibble() %>%
  dplyr::transmute(
    cell_type,
    rho_test_group1 = group1,
    rho_test_group2 = group2,
    rho_test_n_NCR = n_group1,
    rho_test_n_ACR = n_group2,
    rho_p_value = p_value,
    rho_BH_FDR = p_adj,
    rho_tested = is.finite(p_value) & is.finite(p_adj)
  )

heart_dtw_pairs <- heart_dtw$annotation %>%
  tibble::as_tibble() %>%
  dplyr::select(cell_type, Group, combo, n_cell) %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = c(combo, n_cell),
    names_glue = "{.value}_{Group}"
  )

heart_dtw_pairs$paired_NCR_ACR_DTW_distance <- heart_dtw$distance[cbind(
  match(
    heart_dtw_pairs$combo_NCR,
    rownames(heart_dtw$distance)
  ),
  match(
    heart_dtw_pairs$combo_ACR,
    colnames(heart_dtw$distance)
  )
)]

heart_ranking <- heart_rho_summary %>%
  dplyr::left_join(heart_rho_tests, by = "cell_type") %>%
  dplyr::left_join(heart_dtw_pairs, by = "cell_type") %>%
  dplyr::mutate(
    rho_summary_statistic = "mean of saved biological-sample rho values",
    delta_rho_ACR_minus_NCR = rho_mean_ACR - rho_mean_NCR,
    abs_delta_rho = abs(delta_rho_ACR_minus_NCR),
    rho_direction = dplyr::case_when(
      delta_rho_ACR_minus_NCR > 0 ~ "ACR higher rho",
      delta_rho_ACR_minus_NCR < 0 ~ "NCR higher rho",
      TRUE ~ "equal mean rho"
    ),
    rho_test_status = dplyr::case_when(
      !rho_tested ~ "not tested; p-value and FDR are NA",
      rho_BH_FDR < 0.05 ~ "tested; BH FDR < 0.05",
      TRUE ~ "tested; BH FDR >= 0.05"
    ),
    rho_min_cell = heart_rho$params$min_cell,
    dtw_min_cell = heart_dtw$params$min_cell,
    dtw_distance_definition = paste0(
      "same-cell-type NCR versus ACR path-length-normalized distance from ",
      "run_dtw_similarity()$distance"
    ),
    passes_dtw_min_cell = n_cell_NCR >= dtw_min_cell &
      n_cell_ACR >= dtw_min_cell,
    rho_rank = dplyr::min_rank(dplyr::desc(abs_delta_rho)),
    dtw_rank = dplyr::min_rank(
      dplyr::desc(paired_NCR_ACR_DTW_distance)
    ),
    rank_sum = rho_rank + dtw_rank
  ) %>%
  dplyr::arrange(rank_sum, rho_rank, dtw_rank, cell_type) %>%
  dplyr::mutate(
    descriptive_rank = dplyr::row_number(),
    ranking_rule = paste0(
      "ascending rho_rank + dtw_rank; rank 1 within each metric is the ",
      "largest absolute mean-rho change or largest paired raw normalized ",
      "DTW distance; descriptive only"
    )
  ) %>%
  dplyr::select(
    descriptive_rank, rank_sum, rho_rank, dtw_rank, cell_type,
    rho_mean_NCR, rho_mean_ACR,
    delta_rho_ACR_minus_NCR, abs_delta_rho, rho_direction,
    n_sample_NCR, n_sample_ACR,
    mean_n_cell_NCR, mean_n_cell_ACR,
    rho_test_n_NCR, rho_test_n_ACR,
    rho_p_value, rho_BH_FDR, rho_tested, rho_test_status,
    rho_min_cell,
    paired_NCR_ACR_DTW_distance,
    n_cell_NCR, n_cell_ACR, dtw_min_cell, passes_dtw_min_cell,
    rho_summary_statistic, dtw_distance_definition, ranking_rule,
    combo_NCR, combo_ACR
  )

readr::write_csv(
  heart_ranking,
  file.path(review_dir, "scHeart_NCR_ACR_celltype_ranking.csv"),
  na = "unavailable"
)

heart_group_colours <- c(
  "NCR higher rho" = "#5F8EA3",
  "ACR higher rho" = "#A17494",
  "equal mean rho" = "#666666"
)

heart_ranking_plot <- heart_ranking %>%
  dplyr::mutate(
    rho_direction = factor(
      rho_direction,
      levels = c(
        "NCR higher rho",
        "ACR higher rho",
        "equal mean rho"
      )
    )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = abs_delta_rho,
      y = paired_NCR_ACR_DTW_distance,
      colour = rho_direction,
      label = cell_type
    )
  ) +
  ggplot2::geom_point(size = 3.2) +
  ggrepel::geom_text_repel(
    size = 3.1,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = Inf,
    segment.colour = NA,
    max.overlaps = Inf,
    seed = 1,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = heart_group_colours,
    drop = TRUE,
    name = NULL
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0.05, 0.22))
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.08, 0.20))
  ) +
  ggplot2::labs(
    x = "|Mean rho(ACR) - mean rho(NCR)|",
    y = "Paired NCR-ACR DTW distance"
  ) +
  ggplot2::theme_classic(base_size = 10, base_family = "Arial") +
  ggplot2::theme(
    axis.text = ggplot2::element_text(colour = "black"),
    legend.position = "top",
    legend.justification = "center",
    plot.margin = ggplot2::margin(8, 18, 8, 8)
  )

ggplot2::ggsave(
  file.path(ranking_dir, "scHeart_NCR_ACR_celltype_ranking_scatter.png"),
  heart_ranking_plot,
  width = 6.6,
  height = 5.4,
  dpi = 600,
  bg = "white"
)
save_vector_output(
  heart_ranking_plot,
  file.path(ranking_dir, "scHeart_NCR_ACR_celltype_ranking_scatter"),
  6.6,
  5.4
)

heart_input_audit <- tibble::tibble(
  object = c("heart rho result", "heart DTW result"),
  path = c(heart_rho_path, heart_dtw_path)
) %>%
  dplyr::mutate(
    exists = file.exists(path),
    size_bytes = file.info(path)$size,
    modified = as.character(file.info(path)$mtime),
    md5 = unname(tools::md5sum(path))
  )

readr::write_csv(
  heart_input_audit,
  file.path(qa_dir, "scHeart_existing_input_files.csv")
)

writeLines(
  c(
    "scHeart saved-result ranking",
    "",
    "No infer_gene_param(), infer_cell_phase(), run_phase_rho() or run_dtw_similarity() call is made by this script.",
    "Rho change is ACR minus NCR using the group means already saved in rho_sample_summary.",
    "DTW is the same-cell-type NCR-versus-ACR raw normalized distance from the saved distance matrix, not the heatmap similarity colour.",
    "Ranking is descriptive Borda ordering by rho rank plus DTW rank; significance does not enter the rank.",
    "Rows with unavailable tests retain NA p-values and FDR and are never relabelled as ns."
  ),
  file.path(review_dir, "scHeart_METHODS_AND_LIMITS.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(qa_dir, "sessionInfo_scHeart_existing_results.txt")
)

generated_files <- tibble::tibble(
  path = list.files(review_dir, recursive = TRUE, full.names = TRUE)
) %>%
  dplyr::filter(!startsWith(basename(path), "._")) %>%
  dplyr::mutate(
    relative_path = substring(path, nchar(review_dir) + 2),
    size_bytes = file.info(path)$size
  ) %>%
  dplyr::select(relative_path, size_bytes) %>%
  dplyr::arrange(relative_path)

readr::write_csv(
  generated_files,
  file.path(review_dir, "generated_files_manifest.csv")
)
