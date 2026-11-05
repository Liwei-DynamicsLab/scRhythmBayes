suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
load("data/srb_example_data.rda")

count <- srb_example_data$Count
meta <- srb_example_data$meta
design <- unique(meta[c("sample", "condition", "ZT", "replicate")])

stopifnot(
  identical(colnames(count), meta$cell),
  nrow(design) == 36L,
  all(with(design, table(condition, ZT)) == 3L),
  all(table(meta$sample) >= 50L & table(meta$sample) <= 100L),
  identical(sort(unique(srb_example_data$truth$gene$module_truth[1:96])), paste0("M", 0:7))
)

timing <- numeric()
timing["cell_phase"] <- system.time({
  res_phase <- infer_cell_phase(
    count = count,
    metadata = meta,
    ZT_by = "ZT",
    strata_by = "celltype",
    Species = "mouse",
    n_grid = 48,
    n_iter_outer = 3,
    n_iter_seed = 2,
    Min_n = 100,
    setseed = 1,
    Verbose = FALSE
  )
})[["elapsed"]]

phase_error <- abs(
  ((res_phase$theta_pred - meta$Theta_true + pi) %% (2 * pi)) - pi
) * 24 / (2 * pi)

stopifnot(
  all(is.finite(res_phase$theta_pred)),
  all(is.finite(res_phase$r_cell)),
  all(res_phase$r_cell >= 0 & res_phase$r_cell <= 1),
  attr(res_phase, "n_seed_retained") >= 20L,
  stats::median(phase_error) < 1
)

use <- meta$celltype == "1"
timing["sample_gene"] <- system.time({
  res_gene <- infer_gene_param(
    count = count[, use, drop = FALSE],
    metadata = meta[use, , drop = FALSE],
    theta_pred = res_phase$theta_pred[use],
    r_cell = res_phase$r_cell[use],
    strata_by = "condition",
    ZT_by = "ZT",
    B = 5,
    n_worker = 1,
    setseed = 1,
    Verbose = FALSE,
    sample_by = "sample",
    sample_boot_by = "ZT"
  )
})[["elapsed"]]

stopifnot(
  identical(res_gene$inference$mode, "restricted_null_cluster_multiplier_score"),
  all(res_gene$inference$sample_counts == 18L),
  all(is.finite(res_gene$param_long$P_value)),
  all(is.finite(res_gene$param_long$FDR)),
  all(res_gene$fit_diagnostics$Sample_status == "ok_cluster_multiplier_score")
)

timing["module"] <- system.time({
  res_module <- infer_rhythm_module(
    res_gene$param_wide,
    compare_by = "condition",
    group1 = "A",
    group2 = "B",
    Verbose = FALSE
  )
})[["elapsed"]]

module_truth <- merge(
  res_module$module_result[c("gene", "seed_module", "module_map")],
  srb_example_data$truth$gene[c("gene", "module_truth", "role")],
  by = "gene"
)
module_truth <- module_truth[module_truth$role == "module", ]

stopifnot(
  identical(sort(unique(module_truth$seed_module)), 0:7),
  identical(sort(unique(module_truth$module_map)), 0:7)
)

timing["ovr"] <- system.time({
  res_ovr <- infer_ovr_rhythm(
    count = count,
    metadata = meta,
    theta_pred = res_phase$theta_pred,
    r_cell = res_phase$r_cell,
    group_by = "celltype",
    ZT_by = "ZT",
    q_cut = 0.05,
    amp_cut = 0.25,
    n_worker = 1,
    setseed = 1,
    Verbose = FALSE,
    sample_by = "sample",
    sample_boot_by = "ZT",
    B = 3
  )
})[["elapsed"]]

marker_truth <- merge(
  res_ovr$ovr_markers[c("gene", "celltype")],
  srb_example_data$truth$gene[c("gene", "role", "target_celltype")],
  by = "gene",
  all.x = TRUE
)

stopifnot(
  all(res_ovr$ovr_summary$n_marker > 0L),
  any(marker_truth$role == "ovr_marker" & marker_truth$celltype == marker_truth$target_celltype)
)

genes <- head(res_ovr$ovr_markers$gene[res_ovr$ovr_markers$celltype == "1"], 4)
ovr_plot <- plot_ovr_curve(res_ovr, genes = genes, target_group = "1", ncol = 2)
curve_data <- attr(ovr_plot, "curve_data")
stopifnot(inherits(ovr_plot, "ggplot"), all(is.finite(curve_data$fitted_mean)))

dir.create("data-raw/generated", recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    timing_seconds = timing,
    median_phase_error_hours = stats::median(phase_error),
    phase_error_p90_hours = unname(stats::quantile(phase_error, 0.9)),
    n_seed_retained = attr(res_phase, "n_seed_retained"),
    module_seed_accuracy = mean(paste0("M", module_truth$seed_module) == module_truth$module_truth),
    module_map_accuracy = mean(paste0("M", module_truth$module_map) == module_truth$module_truth),
    ovr_summary = res_ovr$ovr_summary,
    ovr_true_marker_count = sum(
      marker_truth$role == "ovr_marker" & marker_truth$celltype == marker_truth$target_celltype
    )
  ),
  "data-raw/generated/e2e_validation.rds"
)

print(timing)
print(res_ovr$ovr_summary)
