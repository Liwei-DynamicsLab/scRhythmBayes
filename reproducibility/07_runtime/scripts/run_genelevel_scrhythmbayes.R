suppressPackageStartupMessages({
  library(qs)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
set.seed(1)

input_file <- "outputs/gene_benchmark/algorithm_results/benchmark_input.qs"
source_file <- "../R/infer_gene_param.R"
output_file <- "07_runtime/results/scRhythmBayes_profile_result.qs"
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_file), file.exists(source_file))
benchmark_input <- qs::qread(input_file)
counts <- benchmark_input$counts
meta <- benchmark_input$metadata
stopifnot(nrow(counts) == 4000L, nrow(meta) == ncol(counts))

source(source_file)
result <- infer_gene_param(
  count = counts,
  metadata = meta,
  theta_pred = meta$theta_obs,
  r_cell = meta$p,
  strata_by = "group",
  ZT_by = "ZT_bin_num",
  sample_by = "sample",
  sample_boot_by = "ZT_bin",
  n_worker = 4,
  setseed = 1,
  Verbose = TRUE
)
qs::qsave(result, output_file)
cat("Saved:", output_file, "\n")
