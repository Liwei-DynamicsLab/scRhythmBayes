# Generated from the primary formal CYCLOPS section; integration appendix excluded.
options(device = 'pdf')

# ---- Rmd code chunk ----
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 8, fig.height = 5)

library(Matrix)
library(matrixStats)
library(torch)
library(ggplot2)

input_dir <- "outputs/simulation/output/scRITMO_extended_simulation/sigma_gradient_1p5_2p5_3p5_4p5_5p5_6p5_empirical_fixed_v1/R_input"
output_dir <- "07_runtime/runs/CYCLOPS"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

SEED <- 1L
N_MODELS <- 80L
MAX_STEPS <- 2050L
LR <- 0.001
BLUNT_PERCENT <- 0.975
SEED_MIN_CV <- 0.14
SEED_MAX_CV <- 0.70
SEED_TOP_M <- 10000L
EIGEN_TOTAL_VAR <- 0.85
EIGEN_MIN_VAR <- 0.06
MAX_EIGEN <- 30L
ALIGN_P_CUTOFF <- 0.05
ALIGN_SHIFT_ITERATIONS <- 192L

set.seed(SEED)
torch_manual_seed(SEED)

wrap_hour <- function(x) (x + 12) %% 24 - 12


# ---- Rmd code chunk ----
count <- as(
  readMM(file.path(input_dir, "count_gene_by_cell.mtx")),
  "CsparseMatrix"
)

genes <- read.csv(
  file.path(input_dir, "genes.csv"),
  stringsAsFactors = FALSE
)

stopifnot("gene" %in% names(genes))
stopifnot(nrow(count) == nrow(genes))
stopifnot(!anyDuplicated(genes$gene))
stopifnot(all(is.finite(count@x)))
stopifnot(all(count@x >= 0))
stopifnot(all(count@x == round(count@x)))

rownames(count) <- genes$gene

cat("Raw count matrix:", nrow(count), "genes x", ncol(count), "cells\n")


# ---- Rmd code chunk ----
library_size <- Matrix::colSums(count)

stopifnot(all(is.finite(library_size)))
stopifnot(all(library_size > 0))

expr_norm <- count %*% Diagonal(x = 1e4 / library_size)
expr_log <- expr_norm
expr_log@x <- log1p(expr_log@x)

expr_matrix <- as.matrix(expr_log)
metadata_index <- read.csv(
  file.path(
    input_dir,
    "metadata.csv"
  ),
  stringsAsFactors = FALSE
)

metadata_index <- metadata_index[
  ,
  c(
    "cell",
    "sample_name"
  )
]

stopifnot(
  nrow(metadata_index) ==
    ncol(count)
)

stopifnot(
  !anyDuplicated(
    metadata_index$cell
  )
)

colnames(count) <- metadata_index$cell
sample_levels <- unique(
  metadata_index$sample_name
)

expr_seed_filter <- sapply(
  sample_levels,
  function(sample_id) {
    rowMeans(
      expr_matrix[
        ,
        metadata_index$sample_name ==
          sample_id,
        drop = FALSE
      ]
    )
  }
)

rownames(expr_seed_filter) <- rownames(
  expr_matrix
)

colnames(expr_seed_filter) <- sample_levels

cat(
  "Cell-level expression:",
  nrow(expr_matrix),
  "genes x",
  ncol(expr_matrix),
  "cells\n"
)

cat(
  "Seed-filter pseudo-bulk expression:",
  nrow(expr_seed_filter),
  "genes x",
  ncol(expr_seed_filter),
  "samples\n"
)

stopifnot(all(is.finite(expr_matrix)))
stopifnot(nrow(expr_matrix) == nrow(count))
stopifnot(ncol(expr_matrix) == ncol(count))

cat("Normalized expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "cells\n")


# ---- Rmd code chunk ----
canonical_gene <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y == "BMAL1"] <- "ARNTL"
  y
}

seed_genes_base <- c(
  "CDKN1A", "MMP14", "CLOCK", "PER2", "PER1", "ARNTL", "POR", "WEE1",
  "SLC16A1", "PER3", "NR1D1", "DBP", "NPAS2", "TSC22D3", "TEF", "TENC1",
  "HLF", "NFIL3", "USP2", "NR1D2", "FMO2", "TSPAN4", "LEO1", "TMEM57",
  "STK35", "LONRF3", "GRAMD4", "COQ10B", "DTX4", "C1ORF51"
)

mouse_alignment_genes <- c(
  "ARNTL", "CLOCK", "NPAS2", "NR1D1", "BHLHE41", "NR1D2", "DBP",
  "CIART", "PER1", "PER3", "TEF", "HLF", "CRY2", "PER2", "CRY1",
  "RORC", "NFIL3"
)

seed_genes_list <- unique(
  c(
    seed_genes_base,
    mouse_alignment_genes
  )
)
mouse_alignment_acrophases <- c(
  ARNTL = 0,
  CLOCK = 0.079,
  NPAS2 = 0.151,
  NR1D1 = 2.296,
  BHLHE41 = 2.909,
  NR1D2 = 2.987,
  DBP = 2.991,
  CIART = 3.008,
  PER1 = 3.122,
  PER3 = 3.306,
  TEF = 3.314,
  HLF = 3.426,
  CRY2 = 3.501,
  PER2 = 3.887,
  CRY1 = 4.995,
  RORC = 5.050,
  NFIL3 = 6.008
)

select_seed_genes <- function(
  expression_data,
  filter_expression_data,
  gene_symbols,
  seed_gene_list,
  top_m = 10000,
  min_cv = 0.14,
  max_cv = 0.70
) {
  stopifnot(
    nrow(expression_data) ==
      length(gene_symbols)
  )

  stopifnot(
    nrow(filter_expression_data) ==
      length(gene_symbols)
  )

  gene_key <- canonical_gene(
    gene_symbols
  )

  filter_means_all <- rowMeans(
    filter_expression_data
  )

  top_n <- min(
    top_m,
    nrow(filter_expression_data)
  )

  top_idx <- order(
    filter_means_all,
    decreasing = TRUE
  )[
    seq_len(top_n)
  ]

  top_mask <- rep(
    FALSE,
    nrow(filter_expression_data)
  )

  top_mask[
    top_idx
  ] <- TRUE

  seed_mask <- (
    gene_key %in%
      seed_gene_list
  ) & top_mask

  if (!any(seed_mask)) {
    stop(
      "No CYCLOPS seed genes were found."
    )
  }

  candidate_filter <- filter_expression_data[
    seed_mask,
    ,
    drop = FALSE
  ]

  candidate_model <- expression_data[
    seed_mask,
    ,
    drop = FALSE
  ]

  candidate_symbols <- gene_symbols[
    seed_mask
  ]

  candidate_key <- gene_key[
    seed_mask
  ]

  means <- rowMeans(
    candidate_filter
  )

  sds <- matrixStats::rowSds(
    candidate_filter
  )

  cvs <- sds / means

  keep <- (
    is.finite(cvs) &
      is.finite(means) &
      means > 0 &
      cvs >= min_cv &
      cvs <= max_cv
  )

  qc <- data.frame(
    gene = candidate_symbols,
    canonical_gene = candidate_key,
    mean_expression = means,
    sd_expression = sds,
    CV = cvs,
    retained = keep,
    stringsAsFactors = FALSE
  )

  list(
    data = candidate_model[
      keep,
      ,
      drop = FALSE
    ],
    genes = candidate_symbols[
      keep
    ],
    canonical_genes = candidate_key[
      keep
    ],
    qc = qc
  )
}

blunt_outliers <- function(
  data,
  percentile = 0.975
) {
  stopifnot(
    percentile > 0.5,
    percentile < 1
  )

  data_blunted <- data

  qc <- data.frame(
    gene = rownames(data),
    nonzero_fraction = numeric(nrow(data)),
    lower_bound = numeric(nrow(data)),
    upper_bound = numeric(nrow(data)),
    sparse_fallback = logical(nrow(data)),
    mean_before = numeric(nrow(data)),
    mean_after = numeric(nrow(data)),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(data))) {
    x <- as.numeric(
      data[i, ]
    )

    lower_bound <- as.numeric(
      quantile(
        x,
        probs = 1 - percentile,
        names = FALSE,
        type = 7
      )
    )

    upper_bound <- as.numeric(
      quantile(
        x,
        probs = percentile,
        names = FALSE,
        type = 7
      )
    )

    sparse_fallback <- FALSE

    if (
      !is.finite(upper_bound) ||
      upper_bound <= 0
    ) {
      positive_values <- x[
        is.finite(x) &
          x > 0
      ]

      if (length(positive_values) == 0) {
        stop(
          paste0(
            "Seed gene ",
            rownames(data)[i],
            " has no positive expression values."
          )
        )
      }

      lower_bound <- 0

      upper_bound <- as.numeric(
        quantile(
          positive_values,
          probs = percentile,
          names = FALSE,
          type = 7
        )
      )

      sparse_fallback <- TRUE
    }

    data_blunted[i, ] <- pmin(
      pmax(
        x,
        lower_bound
      ),
      upper_bound
    )

    qc$nonzero_fraction[i] <- mean(
      x > 0
    )

    qc$lower_bound[i] <- lower_bound

    qc$upper_bound[i] <- upper_bound

    qc$sparse_fallback[i] <- sparse_fallback

    qc$mean_before[i] <- mean(
      x
    )

    qc$mean_after[i] <- mean(
      data_blunted[i, ]
    )
  }

  list(
    data = data_blunted,
    qc = qc
  )
}

mean_normalize <- function(data) {
  means <- rowMeans(data)

  if (any(!is.finite(means)) || any(means <= 0)) {
    stop("Mean normalization encountered a non-positive or non-finite gene mean.")
  }

  sweep(
    sweep(data, 1, means, FUN = "-"),
    1,
    means,
    FUN = "/"
  )
}

perform_svd_cyclops <- function(
  data,
  total_var = 0.85,
  min_individual_var = 0.06,
  max_dims = 30
) {
  s <- svd(data, nu = 0)

  var_fraction <- s$d^2 / sum(s$d^2)
  cumulative_var <- cumsum(var_fraction)

  max_possible <- min(
    max_dims,
    length(s$d),
    ncol(data)
  )

  eligible <- which(
    seq_along(var_fraction) <= max_possible &
      var_fraction >= min_individual_var
  )

  if (length(eligible) == 0) {
    stop("No eigengene satisfied the minimum individual variance criterion.")
  }

  n_by_individual <- max(eligible)

  n_by_total <- which(cumulative_var >= total_var)[1]

  if (is.na(n_by_total)) {
    n_by_total <- max_possible
  }

  n_dims <- min(
    n_by_individual,
    n_by_total,
    max_possible
  )

  if (n_dims < 2) {
    n_dims <- min(2L, max_possible)
    warning(
      "Fewer than two eigengenes met both formal criteria; retaining two dimensions because the circular autoencoder requires a usable multivariate input."
    )
  }

  eigengenes <- 10 * t(
    s$v[, seq_len(n_dims), drop = FALSE]
  )

  variance_table <- data.frame(
    dimension = seq_len(length(var_fraction)),
    variance_fraction = var_fraction,
    cumulative_variance = cumulative_var
  )

  list(
    eigengenes = eigengenes,
    n_dims = n_dims,
    variance_table = variance_table,
    retained_variance = cumulative_var[n_dims]
  )
}

layer_circular_node <- nn_module(
  "CircularNode",
  forward = function(x) {
    nnf_normalize(
      x,
      p = 2,
      dim = 2,
      eps = 1e-12
    )
  }
)

cyclops_module <- nn_module(
  "CYCLOPS_Model",
  initialize = function(n_in) {
    self$n_in <- n_in
    self$l1 <- nn_linear(n_in, 2)
    self$circular <- layer_circular_node()
    self$l2 <- nn_linear(2, n_in)
  },
  forward = function(x) {
    pre_circle <- self$l1(x)
    bottleneck <- self$circular(pre_circle)
    reconstruction <- self$l2(bottleneck)

    list(
      reconstruction = reconstruction,
      bottleneck = bottleneck,
      pre_circle = pre_circle
    )
  }
)

train_cyclops_ensemble <- function(
  eigengenes,
  n_models = 80,
  max_steps = 2050,
  lr = 0.001,
  seed = 1
) {
  x_tensor <- torch_tensor(
    t(eigengenes),
    dtype = torch_float32()
  )

  n_dims <- nrow(eigengenes)
  model_losses <- numeric(n_models)
  best_loss <- Inf
  best_model <- NULL
  best_model_index <- NA_integer_

  cat(
    sprintf(
      "Training CYCLOPS ensemble: %d models x %d steps\n",
      n_models,
      max_steps
    )
  )

  for (i in seq_len(n_models)) {
    torch_manual_seed(seed + i - 1L)

    model <- cyclops_module(
      n_in = n_dims
    )

    optimizer <- optim_adam(
      model$parameters,
      lr = lr,
      betas = c(0.9, 0.999)
    )

    model$train()

    for (step in seq_len(max_steps)) {
      optimizer$zero_grad()

      output <- model(x_tensor)

      loss <- nnf_mse_loss(
        output$reconstruction,
        x_tensor
      )

      loss$backward()
      optimizer$step()
    }

    model$eval()

    final_loss <- with_no_grad({
      output <- model(x_tensor)

      as.numeric(
        nnf_mse_loss(
          output$reconstruction,
          x_tensor
        )$item()
      )
    })

    model_losses[i] <- final_loss

    if (final_loss < best_loss) {
      best_loss <- final_loss
      best_model <- model
      best_model_index <- i
    }

    if (i == 1L || i %% 10L == 0L || i == n_models) {
      cat(
        sprintf(
          "Model %d/%d - final loss: %.6f - best loss: %.6f\n",
          i,
          n_models,
          final_loss,
          best_loss
        )
      )
    }
  }

  list(
    model = best_model,
    best_loss = best_loss,
    best_model_index = best_model_index,
    model_losses = model_losses
  )
}

get_cyclops_phase <- function(model, eigengenes) {
  x_tensor <- torch_tensor(
    t(eigengenes),
    dtype = torch_float32()
  )

  model$eval()

  output <- with_no_grad({
    model(x_tensor)
  })

  bottleneck <- as.matrix(output$bottleneck)
  pre_circle <- as.matrix(output$pre_circle)

  phase <- atan2(
    bottleneck[, 2],
    bottleneck[, 1]
  ) %% (2 * pi)

  magnitude <- sqrt(
    pre_circle[, 1]^2 +
      pre_circle[, 2]^2
  )

  list(
    phase = phase,
    bottleneck = bottleneck,
    magnitude = magnitude
  )
}

cosinor_fit <- function(expression_data, phases) {
  stopifnot(ncol(expression_data) == length(phases))

  design <- cbind(
    intercept = 1,
    sin_phase = sin(phases),
    cos_phase = cos(phases)
  )

  n <- nrow(design)
  p <- ncol(design)
  df_model <- p - 1
  df_resid <- n - p

  results <- vector(
    "list",
    nrow(expression_data)
  )

  for (i in seq_len(nrow(expression_data))) {
    y <- as.numeric(
      expression_data[i, ]
    )

    fit <- lm.fit(
      x = design,
      y = y
    )

    coefficients <- fit$coefficients
    residuals <- fit$residuals

    rss <- sum(residuals^2)
    tss <- sum(
      (y - mean(y))^2
    )

    if (!is.finite(tss) || tss <= 0 || !is.finite(rss)) {
      r_squared <- 0
      p_value <- 1
    } else {
      r_squared <- max(
        0,
        1 - rss / tss
      )

      ms_model <- max(
        0,
        (tss - rss) / df_model
      )

      ms_resid <- rss / df_resid

      if (ms_resid <= 0) {
        p_value <- 0
      } else {
        f_stat <- ms_model / ms_resid
        p_value <- pf(
          f_stat,
          df_model,
          df_resid,
          lower.tail = FALSE
        )
      }
    }

    sin_coef <- coefficients[2]
    cos_coef <- coefficients[3]

    results[[i]] <- data.frame(
      Gene_Symbols = rownames(expression_data)[i],
      Sine_Coefficient = sin_coef,
      Cosine_Coefficient = cos_coef,
      Fit_Average = coefficients[1],
      P_Value = p_value,
      R_Squared = r_squared,
      Acrophase = atan2(
        sin_coef,
        cos_coef
      ) %% (2 * pi),
      Amplitude = sqrt(
        sin_coef^2 +
          cos_coef^2
      ),
      stringsAsFactors = FALSE
    )
  }

  output <- do.call(
    rbind,
    results
  )

  output$BHQ_Statistic <- p.adjust(
    output$P_Value,
    method = "BH"
  )

  output
}

align_cyclops_mouse_clock <- function(
  predicted_phases,
  expression_data,
  gene_symbols,
  marker_acrophases,
  p_cutoff = 0.05,
  shift_iterations = 192
) {
  gene_key <- canonical_gene(
    gene_symbols
  )

  available_markers <- intersect(
    names(marker_acrophases),
    unique(gene_key)
  )

  if (length(available_markers) < 2) {
    stop("Fewer than two CYCLOPS mouse alignment genes are available.")
  }

  marker_idx <- vapply(
    available_markers,
    function(g) {
      which(gene_key == g)[1]
    },
    integer(1)
  )

  marker_data <- expression_data[
    marker_idx,
    ,
    drop = FALSE
  ]

  rownames(marker_data) <- available_markers

  marker_fit <- cosinor_fit(
    marker_data,
    predicted_phases
  )

  marker_fit$canonical_gene <- marker_fit$Gene_Symbols
  marker_fit$expected_acrophase <- marker_acrophases[
    marker_fit$canonical_gene
  ]

  marker_fit$used_for_alignment <- (
    is.finite(marker_fit$P_Value) &
      marker_fit$P_Value < p_cutoff &
      is.finite(marker_fit$Acrophase)
  )

  used <- marker_fit[
    marker_fit$used_for_alignment,
    ,
    drop = FALSE
  ]

  if (nrow(used) < 2) {
    stop(
      paste0(
        "Fewer than two native core-clock markers passed the CYCLOPS alignment p-value cutoff. ",
        "Inspect CYCLOPS_alignment_markers.csv; no truth-based or ZT-based fallback is applied."
      )
    )
  }

  shifts <- seq(
    0,
    2 * pi,
    length.out = shift_iterations + 1L
  )

  shifts <- shifts[
    -length(shifts)
  ]

  directions <- c(
    1,
    -1
  )

  best_error <- Inf
  best_shift <- NA_real_
  best_direction <- NA_real_

  observed <- used$Acrophase
  expected <- used$expected_acrophase

  for (direction in directions) {
    for (shift in shifts) {
      shifted_observed <- (
        direction * observed +
          shift
      ) %% (2 * pi)

      delta <- atan2(
        sin(
          shifted_observed -
            expected
        ),
        cos(
          shifted_observed -
            expected
        )
      )

      error <- mean(
        delta^2
      )

      if (error < best_error) {
        best_error <- error
        best_shift <- shift
        best_direction <- direction
      }
    }
  }

  aligned_phases <- (
    best_direction *
      predicted_phases +
      best_shift
  ) %% (2 * pi)

  list(
    phase = aligned_phases,
    best_error = best_error,
    best_shift = best_shift,
    best_direction = best_direction,
    marker_fit = marker_fit
  )
}


# ---- Rmd code chunk ----
seed_result <- select_seed_genes(
  expression_data = expr_matrix,
  filter_expression_data = expr_seed_filter,
  gene_symbols = rownames(expr_matrix),
  seed_gene_list = seed_genes_list,
  top_m = SEED_TOP_M,
  min_cv = SEED_MIN_CV,
  max_cv = SEED_MAX_CV
)

write.csv(
  seed_result$qc,
  file.path(
    output_dir,
    "CYCLOPS_seed_gene_qc.csv"
  ),
  row.names = FALSE
)

print(
  seed_result$qc[
    order(
      seed_result$qc$CV
    ),
  ]
)

cat(
  "Candidate CYCLOPS seed genes:",
  nrow(seed_result$qc),
  "\n"
)

cat(
  "Retained CYCLOPS seed genes:",
  nrow(seed_result$data),
  "\n"
)

if (
  nrow(seed_result$data) < 2
) {
  stop(
    paste0(
      "Fewer than two seed genes passed CYCLOPS seed filtering. ",
      "Inspect CYCLOPS_seed_gene_qc.csv."
    )
  )
}


# ---- Rmd code chunk ----
seed_blunted <- blunt_outliers(
  seed_result$data,
  percentile = BLUNT_PERCENT
)

seed_normalized <- mean_normalize(
  seed_blunted$data
)

stopifnot(
  all(
    is.finite(
      seed_normalized
    )
  )
)

cat(
  "CYCLOPS preprocessing matrix:",
  nrow(seed_normalized),
  "seed genes x",
  ncol(seed_normalized),
  "cells\n"
)


# ---- Rmd code chunk ----
svd_result <- perform_svd_cyclops(
  data = seed_normalized,
  total_var = EIGEN_TOTAL_VAR,
  min_individual_var = EIGEN_MIN_VAR,
  max_dims = MAX_EIGEN
)

eigengenes <- svd_result$eigengenes

write.csv(
  svd_result$variance_table,
  file.path(
    output_dir,
    "CYCLOPS_eigengene_variance.csv"
  ),
  row.names = FALSE
)

cat(
  "Retained eigengenes:",
  svd_result$n_dims,
  "\n"
)

cat(
  "Cumulative variance retained:",
  round(
    svd_result$retained_variance,
    4
  ),
  "\n"
)

print(
  head(
    svd_result$variance_table,
    10
  )
)


# ---- Rmd code chunk ----
variance_plot_data <- svd_result$variance_table[
  seq_len(
    min(
      15,
      nrow(
        svd_result$variance_table
      )
    )
  ),
]

p_eigen <- ggplot(
  variance_plot_data,
  aes(
    x = dimension,
    y = variance_fraction
  )
) +
  geom_col() +
  geom_hline(
    yintercept = EIGEN_MIN_VAR,
    linetype = "dashed"
  ) +
  labs(
    x = "SVD dimension",
    y = "Variance fraction",
    title = "CYCLOPS eigengene variance"
  ) +
  theme_classic()

p_eigen

ggsave(
  file.path(
    output_dir,
    "CYCLOPS_eigengene_variance.png"
  ),
  p_eigen,
  width = 6,
  height = 4,
  dpi = 300,
  bg = "white"
)


# ---- Rmd code chunk ----
training_result <- train_cyclops_ensemble(
  eigengenes = eigengenes,
  n_models = N_MODELS,
  max_steps = MAX_STEPS,
  lr = LR,
  seed = SEED
)

cyclops_model <- training_result$model

model_loss_table <- data.frame(
  model = seq_len(N_MODELS),
  final_reconstruction_loss = training_result$model_losses
)

write.csv(
  model_loss_table,
  file.path(
    output_dir,
    "CYCLOPS_model_losses.csv"
  ),
  row.names = FALSE
)

cat(
  "Best model:",
  training_result$best_model_index,
  "\n"
)

cat(
  "Best reconstruction loss:",
  training_result$best_loss,
  "\n"
)


# ---- Rmd code chunk ----
p_loss <- ggplot(
  model_loss_table,
  aes(
    x = model,
    y = final_reconstruction_loss
  )
) +
  geom_point() +
  geom_hline(
    yintercept = training_result$best_loss,
    linetype = "dashed"
  ) +
  labs(
    x = "Ensemble model",
    y = "Final reconstruction loss",
    title = "CYCLOPS ensemble training"
  ) +
  theme_classic()

p_loss

ggsave(
  file.path(
    output_dir,
    "CYCLOPS_model_losses.png"
  ),
  p_loss,
  width = 6,
  height = 4,
  dpi = 300,
  bg = "white"
)


# ---- Rmd code chunk ----
raw_prediction <- get_cyclops_phase(
  cyclops_model,
  eigengenes
)

raw_phase <- raw_prediction$phase
raw_magnitude <- raw_prediction$magnitude

stopifnot(
  length(raw_phase) == ncol(expr_matrix)
)

stopifnot(
  all(is.finite(raw_phase))
)

stopifnot(
  all(raw_phase >= 0 & raw_phase < 2 * pi)
)

cat(
  "Raw phase range:",
  range(raw_phase),
  "\n"
)

cat(
  "Raw magnitude range:",
  range(raw_magnitude),
  "\n"
)


# ---- Rmd code chunk ----
raw_phase_df <- data.frame(
  phase = raw_phase
)

p_raw <- ggplot(
  raw_phase_df,
  aes(x = phase)
) +
  geom_histogram(
    bins = 100,
    boundary = 0
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      2 * pi,
      length.out = 5
    ),
    labels = c(
      "0",
      expression(pi / 2),
      expression(pi),
      expression(3 * pi / 2),
      expression(2 * pi)
    )
  ) +
  labs(
    x = "Raw CYCLOPS phase (rad)",
    y = "Cells",
    title = "CYCLOPS raw phase"
  ) +
  theme_classic()

p_raw

ggsave(
  file.path(
    output_dir,
    "CYCLOPS_raw_phase_histogram.png"
  ),
  p_raw,
  width = 6,
  height = 4,
  dpi = 300,
  bg = "white"
)


# ---- Rmd code chunk ----
alignment_result <- align_cyclops_mouse_clock(
  predicted_phases = raw_phase,
  expression_data = expr_matrix,
  gene_symbols = rownames(expr_matrix),
  marker_acrophases = mouse_alignment_acrophases,
  p_cutoff = ALIGN_P_CUTOFF,
  shift_iterations = ALIGN_SHIFT_ITERATIONS
)

aligned_phase <- alignment_result$phase

write.csv(
  alignment_result$marker_fit,
  file.path(
    output_dir,
    "CYCLOPS_alignment_markers.csv"
  ),
  row.names = FALSE
)

alignment_summary <- data.frame(
  metric = c(
    "best_direction",
    "best_shift_rad",
    "best_shift_hour",
    "alignment_mean_squared_circular_error",
    "available_alignment_markers",
    "significant_alignment_markers"
  ),
  value = c(
    alignment_result$best_direction,
    alignment_result$best_shift,
    alignment_result$best_shift * 24 / (2 * pi),
    alignment_result$best_error,
    nrow(alignment_result$marker_fit),
    sum(
      alignment_result$marker_fit$used_for_alignment
    )
  )
)

write.csv(
  alignment_summary,
  file.path(
    output_dir,
    "CYCLOPS_alignment_summary.csv"
  ),
  row.names = FALSE
)

alignment_summary

alignment_result$marker_fit


# ---- Rmd code chunk ----
metadata <- read.csv(
  file.path(
    input_dir,
    "metadata.csv"
  ),
  stringsAsFactors = FALSE
)

meta_cols <- c(
  "cell",
  "sample_name",
  "ZT",
  "celltype",
  "sigma_bio_hour",
  "T_true",
  "Theta_true"
)

stopifnot(
  all(
    meta_cols %in%
      names(metadata)
  )
)

stopifnot(
  nrow(metadata) ==
    length(aligned_phase)
)

stopifnot(
  !anyDuplicated(
    metadata$cell
  )
)


# ---- Rmd code chunk ----
CYCLOPS_cell_phase <- metadata[
  ,
  meta_cols
]

CYCLOPS_cell_phase$raw_phase_rad <- raw_phase

CYCLOPS_cell_phase$raw_phase_hour <- (
  raw_phase *
    24 /
    (2 * pi)
)

CYCLOPS_cell_phase$phase_magnitude <- raw_magnitude

CYCLOPS_cell_phase$pred_rad <- aligned_phase

CYCLOPS_cell_phase$pred_hour <- (
  aligned_phase *
    24 /
    (2 * pi)
)

CYCLOPS_cell_phase$wrapped_error_hour <- wrap_hour(
  CYCLOPS_cell_phase$pred_hour -
    CYCLOPS_cell_phase$T_true
)

CYCLOPS_cell_phase$abs_error_hour <- abs(
  CYCLOPS_cell_phase$wrapped_error_hour
)

stopifnot(
  all(
    is.finite(
      CYCLOPS_cell_phase$pred_hour
    )
  )
)

stopifnot(
  all(
    CYCLOPS_cell_phase$pred_hour >= 0 &
      CYCLOPS_cell_phase$pred_hour < 24
  )
)

write.csv(
  CYCLOPS_cell_phase,
  file.path(
    output_dir,
    "CYCLOPS_cell_phase.csv"
  ),
  row.names = FALSE
)

head(
  CYCLOPS_cell_phase
)


# ---- Rmd code chunk ----
groups <- c(
  "overall",
  "CT1",
  "CT2",
  "CT3",
  "CT4",
  "CT5",
  "CT6"
)

CYCLOPS_summary <- do.call(
  rbind,
  lapply(
    groups,
    function(group) {
      if (group == "overall") {
        x <- CYCLOPS_cell_phase
      } else {
        x <- CYCLOPS_cell_phase[
          CYCLOPS_cell_phase$celltype ==
            group,
          ,
          drop = FALSE
        ]
      }

      error <- x$wrapped_error_hour

      data.frame(
        group = group,
        n_cells = nrow(x),
        circular_MAE_hour = mean(
          abs(error)
        ),
        circular_MedAE_hour = median(
          abs(error)
        ),
        circular_RMSE_hour = sqrt(
          mean(
            error^2
          )
        )
      )
    }
  )
)

write.csv(
  CYCLOPS_summary,
  file.path(
    output_dir,
    "CYCLOPS_circular_error_summary.csv"
  ),
  row.names = FALSE
)

CYCLOPS_summary


# ---- Rmd code chunk ----
df_cyclops <- data.frame(
  true = CYCLOPS_cell_phase$T_true %% 24,
  pred = CYCLOPS_cell_phase$pred_hour %% 24
)

color_point_cyclops <- "#D6A84B"
color_diag <- "gray55"
color_band2 <- "#F7DDE3"
color_band4 <- "#FCECEF"

x_seq <- seq(
  -6,
  30,
  length.out = 1000
)

ribbon_data <- expand.grid(
  x = x_seq,
  offset = c(
    -24,
    0,
    24
  )
)

ribbon_data$y <- (
  ribbon_data$x +
    ribbon_data$offset
)

ribbon_data$y_low2 <- (
  ribbon_data$y - 2
)

ribbon_data$y_high2 <- (
  ribbon_data$y + 2
)

ribbon_data$y_low4 <- (
  ribbon_data$y - 4
)

ribbon_data$y_high4 <- (
  ribbon_data$y + 4
)

p_cyclops <- ggplot(
  df_cyclops,
  aes(
    x = true,
    y = pred
  )
) +
  geom_ribbon(
    data = ribbon_data,
    aes(
      x = x,
      ymin = y_low4,
      ymax = y_low2,
      group = offset
    ),
    fill = color_band4,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = ribbon_data,
    aes(
      x = x,
      ymin = y_low2,
      ymax = y_high2,
      group = offset
    ),
    fill = color_band2,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = ribbon_data,
    aes(
      x = x,
      ymin = y_high2,
      ymax = y_high4,
      group = offset
    ),
    fill = color_band4,
    inherit.aes = FALSE
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = color_diag,
    linewidth = 0.45
  ) +
  geom_abline(
    intercept = 24,
    slope = 1,
    linetype = "dashed",
    color = color_diag,
    linewidth = 0.45
  ) +
  geom_abline(
    intercept = -24,
    slope = 1,
    linetype = "dashed",
    color = color_diag,
    linewidth = 0.45
  ) +
  geom_point(
    size = 0.65,
    color = color_point_cyclops,
    alpha = 0.60
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      24,
      2
    ),
    expand = c(
      0,
      0
    )
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      24,
      2
    ),
    expand = c(
      0,
      0
    )
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(
      0,
      24
    ),
    ylim = c(
      0,
      24
    ),
    expand = FALSE,
    clip = "on"
  ) +
  labs(
    x = "True Time (h)",
    y = "Predicted Time (h)",
    title = "CYCLOPS"
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "serif"
    ),
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 13,
      face = "bold"
    ),
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),
    panel.grid = element_blank(),
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  )

p_cyclops

ggsave(
  file.path(
    output_dir,
    "FigureC_CYCLOPS.png"
  ),
  p_cyclops,
  width = 4,
  height = 4,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "FigureC_CYCLOPS.pdf"
  ),
  p_cyclops,
  width = 4,
  height = 4,
  device = cairo_pdf,
  bg = "white"
)


# ---- Rmd code chunk ----
run_parameters <- data.frame(
  parameter = c(
    "seed",
    "n_models",
    "max_steps",
    "learning_rate",
    "blunt_percent",
    "seed_min_CV",
    "seed_max_CV",
    "seed_top_m",
    "eigen_total_variance",
    "eigen_min_individual_variance",
    "max_eigengenes",
    "alignment_p_cutoff",
    "alignment_shift_iterations",
    "retained_seed_genes",
    "retained_eigengenes",
    "best_model_index",
    "best_reconstruction_loss"
  ),
  value = c(
    SEED,
    N_MODELS,
    MAX_STEPS,
    LR,
    BLUNT_PERCENT,
    SEED_MIN_CV,
    SEED_MAX_CV,
    SEED_TOP_M,
    EIGEN_TOTAL_VAR,
    EIGEN_MIN_VAR,
    MAX_EIGEN,
    ALIGN_P_CUTOFF,
    ALIGN_SHIFT_ITERATIONS,
    nrow(seed_result$data),
    svd_result$n_dims,
    training_result$best_model_index,
    training_result$best_loss
  )
)

write.csv(
  run_parameters,
  file.path(
    output_dir,
    "CYCLOPS_run_parameters.csv"
  ),
  row.names = FALSE
)

run_parameters


# ---- Rmd code chunk ----
cat(
  "Cells:",
  nrow(CYCLOPS_cell_phase),
  "\n"
)

cat(
  "Candidate seed genes:",
  nrow(seed_result$qc),
  "\n"
)

cat(
  "Retained seed genes:",
  nrow(seed_result$data),
  "\n"
)

cat(
  "Retained eigengenes:",
  svd_result$n_dims,
  "\n"
)

cat(
  "Significant alignment markers:",
  sum(
    alignment_result$marker_fit$used_for_alignment
  ),
  "\n"
)

cat(
  "Circular MAE:",
  CYCLOPS_summary$circular_MAE_hour[
    CYCLOPS_summary$group == "overall"
  ],
  "\n"
)

cat(
  "Circular MedAE:",
  CYCLOPS_summary$circular_MedAE_hour[
    CYCLOPS_summary$group == "overall"
  ],
  "\n"
)

cat(
  "Circular RMSE:",
  CYCLOPS_summary$circular_RMSE_hour[
    CYCLOPS_summary$group == "overall"
  ],
  "\n"
)

cat(
  "Result:",
  file.path(
    output_dir,
    "CYCLOPS_cell_phase.csv"
  ),
  "\n"
)

