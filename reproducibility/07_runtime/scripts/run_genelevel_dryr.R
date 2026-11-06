suppressPackageStartupMessages({
  library(qs)
  library(dryR)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
set.seed(1)

input_file <- "outputs/gene_benchmark/algorithm_results/benchmark_input.qs"
output_file <- "07_runtime/results/dryR_profile_result.qs"
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_file))
benchmark_input <- qs::qread(input_file)
count_mat <- benchmark_input$pseudobulk_counts
sample_info <- benchmark_input$sample_info
stopifnot(nrow(count_mat) == 4000L, ncol(count_mat) == nrow(sample_info))

result <- dryR::dryseq(
  countData = count_mat,
  group = sample_info$group,
  time = sample_info$ZT_bin_num,
  period = 24,
  n.cores = 4
)
qs::qsave(result, output_file)
cat("Saved:", output_file, "\n")
