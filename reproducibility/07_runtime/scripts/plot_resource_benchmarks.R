suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(tidyr)
})

source("/path/to/scRhythmBayes/reproducibility/07_runtime/scripts/panel_alignment_compat.R")

root <- "/path/to/scRhythmBayes/reproducibility/07_runtime"
data_dir <- file.path(root, "data")
metric_dir <- file.path(root, "metrics")
figure_dir <- file.path(root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

method_order <- c("scRhythmBayes", "scRITMO", "Tempo", "CoPhaser", "CYCLOPS", "reCAT", "Cyclum", "scPrisma")
method_colors <- c(
  scRhythmBayes = "#A64078", scRITMO = "#6E87B0", Tempo = "#B2C1D6",
  CoPhaser = "#D3CDDA", CYCLOPS = "#C0D87B", reCAT = "#C7B5CC",
  Cyclum = "#FBE69E", scPrisma = "#B99BB5"
)
gene_colors <- c(scRhythmBayes = "#C8D5B9", dryR = "#F8DDDA")
cpu_fill <- "#B65A45"

theme_resource <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 24, hjust = 0.5, margin = ggplot2::margin(b = 16)),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold", size = 20),
      axis.text = ggplot2::element_text(colour = "#20242A", size = 18),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      legend.text = ggplot2::element_text(size = 16),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "#E8E8E8", linewidth = 0.3),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

save_plot_set <- function(plot, stem, width, height) {
  svg_file <- file.path(figure_dir, paste0(stem, ".svg"))
  pdf_file <- file.path(figure_dir, paste0(stem, ".pdf"))
  svglite::svglite(svg_file, width = width, height = height, bg = "white")
  print(plot)
  grDevices::dev.off()
  native_pdf <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(native_pdf, width = width, height = height, family = "sans", bg = "white")
  print(plot)
  grDevices::dev.off()
  unlink(native_pdf)
  pdf_tmp_dir <- tempfile(pattern = "svg-pdf-")
  dir.create(pdf_tmp_dir)
  soffice <- Sys.which("soffice")
  if (!nzchar(soffice)) stop("LibreOffice is required for the auditable SVG-to-PDF export")
  status <- system2(soffice, c("--headless", "--convert-to", "pdf", "--outdir", shQuote(pdf_tmp_dir), shQuote(svg_file)), stdout = FALSE, stderr = FALSE)
  converted_pdf <- file.path(pdf_tmp_dir, paste0(stem, ".pdf"))
  if (status != 0 || !file.exists(converted_pdf) || !file.copy(converted_pdf, pdf_file, overwrite = TRUE)) {
    stop("SVG-to-PDF conversion failed for ", stem)
  }
  unlink(pdf_tmp_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, device = ragg::agg_png, bg = "white")
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height, dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
  organized_dir <- file.path(figure_dir, if (endsWith(stem, "_3panel")) "combined_panels" else "single_figures")
  dir.create(organized_dir, recursive = TRUE, showWarnings = FALSE)
  exports <- paste0(stem, c(".pdf", ".svg", ".png", ".tiff"))
  stopifnot(all(file.copy(file.path(figure_dir, exports), file.path(organized_dir, exports), overwrite = TRUE)))
}

cell <- readr::read_csv(file.path(data_dir, "cellphase_resource_current_evidence.csv"), show_col_types = FALSE, na = c("", "NA", " ")) %>%
  dplyr::mutate(method = factor(method, levels = rev(method_order)))

cell_runtime_max <- max(cell$runtime_minutes, na.rm = TRUE)
p_cell_runtime <- ggplot2::ggplot(cell, ggplot2::aes(y = method)) +
  ggplot2::geom_col(
    data = dplyr::filter(cell, is.finite(runtime_minutes)),
    ggplot2::aes(x = runtime_minutes, fill = method), width = 0.68, colour = "#4B4B4B", linewidth = 0.25, show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = dplyr::filter(cell, is.finite(runtime_minutes)),
    ggplot2::aes(x = runtime_minutes + cell_runtime_max * 0.012, label = sprintf("%.2f", runtime_minutes)),
    hjust = 0, size = 18 / ggplot2::.pt, family = "sans"
  ) +
  ggplot2::geom_text(
    data = dplyr::filter(cell, !is.finite(runtime_minutes)),
    ggplot2::aes(x = 0.5, label = "Not recorded"),
    hjust = 0, size = 2.7, colour = "#8A8A8A", family = "sans", fontface = "italic"
  ) +
  ggplot2::scale_fill_manual(values = method_colors) +
  ggplot2::scale_y_discrete(limits = rev(method_order), drop = FALSE) +
  ggplot2::scale_x_continuous(limits = c(0, cell_runtime_max * 1.2), expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::labs(
    title = "Runtime",
    x = "Runtime (min)", y = NULL
  ) + theme_resource() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())

if (any(is.finite(cell$peak_rss_gb))) {
  memory_max <- max(cell$peak_rss_gb, na.rm = TRUE)
  p_cell_memory <- ggplot2::ggplot(cell, ggplot2::aes(y = method)) +
    ggplot2::geom_col(
      data = dplyr::filter(cell, is.finite(peak_rss_gb)),
      ggplot2::aes(x = peak_rss_gb, fill = method), width = 0.68, colour = "#4B4B4B", linewidth = 0.25, show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(cell, is.finite(peak_rss_gb)),
      ggplot2::aes(x = peak_rss_gb + memory_max * 0.012, label = sprintf("%.2f", peak_rss_gb)),
      hjust = 0, size = 18 / ggplot2::.pt, family = "sans"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(cell, !is.finite(peak_rss_gb)),
      ggplot2::aes(x = memory_max * 0.03, label = "Not recorded"),
      hjust = 0, size = 2.7, colour = "#8A8A8A", family = "sans", fontface = "italic"
    ) +
    ggplot2::scale_fill_manual(values = method_colors) +
    ggplot2::scale_y_discrete(limits = rev(method_order), drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, memory_max * 1.2), expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::labs(title = "Peak memory", x = "Peak aggregate RSS (GB)", y = NULL) +
    theme_resource() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
} else {
  p_cell_memory <- ggplot2::ggplot(cell, ggplot2::aes(y = method)) +
    ggplot2::geom_text(ggplot2::aes(x = 0.04, label = "Not recorded"), hjust = 0, size = 2.7, colour = "#8A8A8A", family = "sans", fontface = "italic") +
    ggplot2::scale_y_discrete(limits = rev(method_order), drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = NULL) +
    ggplot2::labs(title = "Peak memory", x = NULL, y = NULL) +
    theme_resource() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}

cell_core <- cell %>%
  dplyr::transmute(method, cores = mean_effective_cpu_cores)
core_max <- max(cell_core$cores, na.rm = TRUE)
p_cell_cores <- ggplot2::ggplot(cell_core, ggplot2::aes(y = method, x = cores)) +
  ggplot2::geom_col(
    data = dplyr::filter(cell_core, is.finite(cores)),
    width = 0.68, fill = cpu_fill, colour = "#4B4B4B", linewidth = 0.25
  ) +
  ggplot2::geom_text(
    data = dplyr::filter(cell_core, is.finite(cores)),
    ggplot2::aes(x = cores + core_max * 0.03, label = sprintf("%.2f", cores)),
    hjust = 0, size = 18 / ggplot2::.pt, family = "sans"
  ) +
  ggplot2::scale_y_discrete(limits = rev(method_order), drop = FALSE) +
  ggplot2::scale_x_continuous(limits = c(0, core_max * 1.25), breaks = 0:6, expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::labs(title = "CPU cores used", x = "CPU cores", y = NULL) +
  theme_resource() + ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position = "none"
  )

cell_combined <- p_cell_runtime / p_cell_memory / p_cell_cores + patchwork::plot_layout(heights = c(1, 1, 1))
require_patchwork_panel_alignment(
  plot = cell_combined,
  manifest_path = file.path(root, "qa", "CellPhase_resource_benchmark_3panel.panel-layout.json"),
  report_path = file.path(root, "qa", "CellPhase_resource_benchmark_3panel.panel-alignment.json"),
  width_in = 10.5,
  height_in = 18.6,
  panel_ids = c("A", "B", "C"),
  column_groups = list(c("A", "B", "C")),
  audit_script = "/path/to/external_plot_helpers/audit_panel_alignment.py",
  overlay_svg = file.path(root, "qa", "CellPhase_resource_benchmark_3panel.panel-alignment.svg")
)
save_plot_set(p_cell_runtime, "CellPhase_runtime_barplot", 10.5, 6.2)
save_plot_set(p_cell_memory, "CellPhase_peak_memory_barplot", 10.5, 6.2)
save_plot_set(p_cell_cores, "CellPhase_core_usage_barplot", 10.5, 6.2)
save_plot_set(cell_combined, "CellPhase_resource_benchmark_3panel", 10.5, 18.6)

read_metric <- function(file, method, input_level, workload_note) {
  x <- jsonlite::read_json(file, simplifyVector = TRUE)
  tibble::tibble(
    method = method,
    input_level = input_level,
    workload_note = workload_note,
    status = x$status,
    wall_seconds = x$wall_seconds,
    wall_minutes = x$wall_minutes,
    peak_process_tree_rss_mb = x$peak_process_tree_rss_mb,
    peak_process_tree_rss_gb = x$peak_process_tree_rss_gb,
    configured_cores = x$configured_cores,
    mean_effective_cpu_cores = x$mean_effective_cpu_cores,
    peak_effective_cpu_cores = x$peak_effective_cpu_cores,
    sample_interval_seconds = x$sample_interval_seconds,
    n_samples = x$n_samples,
    started_utc = x$started_utc,
    finished_utc = x$finished_utc,
    rss_definition = x$rss_definition,
    cpu_definition = x$cpu_definition,
    hardware = paste(x$hardware$machdep.cpu.brand_string, paste0(round(as.numeric(x$hardware$hw.memsize) / 1024^3), " GB RAM"), x$hardware$platform, sep = "; ")
  )
}

gene <- dplyr::bind_rows(
  read_metric(
    file.path(metric_dir, "genelevel_scRhythmBayes.json"), "scRhythmBayes", "single-cell counts: 4000 genes x 3600 cells",
    "two strata; 4 workers; sample-aware bootstrap: 100 genes/stratum, B=499"
  ),
  read_metric(
    file.path(metric_dir, "genelevel_dryR.json"), "dryR", "sample-pseudobulk counts: 4000 genes x 36 samples",
    "native dryseq pipeline; four configured cores"
  )
) %>%
  dplyr::mutate(method = factor(method, levels = c("scRhythmBayes", "dryR")))
stopifnot(all(gene$status == "complete"))
readr::write_csv(gene, file.path(data_dir, "genelevel_resource_measured.csv"), na = "")

p_gene_runtime_offset <- max(gene$wall_minutes) * 0.035
p_gene_runtime <- ggplot2::ggplot(gene, ggplot2::aes(x = method, y = wall_minutes, fill = method)) +
  ggplot2::geom_col(width = 0.62, colour = "#4B4B4B", linewidth = 0.3, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(y = wall_minutes + p_gene_runtime_offset, label = sprintf("%.2f", wall_minutes)), vjust = 0.5, size = 18 / ggplot2::.pt, family = "sans") +
  ggplot2::scale_fill_manual(values = gene_colors) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::labs(title = "Runtime", x = NULL, y = "Runtime (min)") +
  theme_resource() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())

p_gene_memory_offset <- max(gene$peak_process_tree_rss_gb) * 0.045
p_gene_memory <- ggplot2::ggplot(gene, ggplot2::aes(x = method, y = peak_process_tree_rss_gb, fill = method)) +
  ggplot2::geom_col(width = 0.62, colour = "#4B4B4B", linewidth = 0.3, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(y = peak_process_tree_rss_gb + p_gene_memory_offset, label = sprintf("%.2f", peak_process_tree_rss_gb)), vjust = 0.5, size = 18 / ggplot2::.pt, family = "sans") +
  ggplot2::scale_fill_manual(values = gene_colors) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::labs(title = "Peak memory", x = NULL, y = "Peak aggregate RSS (GB)") +
  theme_resource() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())

gene_core <- gene %>%
  dplyr::transmute(method, cores = mean_effective_cpu_cores)
p_gene_cores <- ggplot2::ggplot(gene_core, ggplot2::aes(x = method, y = cores)) +
  ggplot2::geom_col(width = 0.62, fill = cpu_fill, colour = "#4B4B4B", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(y = cores + 0.18, label = sprintf("%.2f", cores)),
    vjust = 0.5, size = 18 / ggplot2::.pt, family = "sans"
  ) +
  ggplot2::scale_y_continuous(limits = c(0, max(gene_core$cores) * 1.18), breaks = 0:5, expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::labs(title = "CPU cores used", x = NULL, y = "CPU cores") +
  theme_resource() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())

gene_combined <- p_gene_runtime | p_gene_memory | p_gene_cores
require_patchwork_panel_alignment(
  plot = gene_combined,
  manifest_path = file.path(root, "qa", "GeneLevel_resource_benchmark_3panel.panel-layout.json"),
  report_path = file.path(root, "qa", "GeneLevel_resource_benchmark_3panel.panel-alignment.json"),
  width_in = 22.8,
  height_in = 7.2,
  panel_ids = c("A", "B", "C"),
  row_groups = list(c("A", "B", "C")),
  audit_script = "/path/to/external_plot_helpers/audit_panel_alignment.py",
  overlay_svg = file.path(root, "qa", "GeneLevel_resource_benchmark_3panel.panel-alignment.svg")
)
save_plot_set(p_gene_runtime, "GeneLevel_runtime_barplot", 6.8, 6.4)
save_plot_set(p_gene_memory, "GeneLevel_peak_memory_barplot", 6.8, 6.4)
save_plot_set(p_gene_cores, "GeneLevel_core_usage_barplot", 8.4, 6.4)
save_plot_set(gene_combined, "GeneLevel_resource_benchmark_3panel", 22.8, 7.2)

writeLines(capture.output(sessionInfo()), file.path(root, "R_sessionInfo.txt"))
