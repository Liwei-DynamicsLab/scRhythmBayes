#' Infer gene-level rhythmic parameters
#'
#' This function estimates gene-level rhythmic parameters within each metadata
#' stratum using fixed cell phase and cell-level rhythmic probability.
#'
#' @param count A gene-by-cell raw count matrix. Rows are genes and columns are
#' cells. It can be a base R matrix or a sparse \code{Matrix} object.
#' @param metadata A data.frame containing cell-level metadata.
#' @param theta_pred A numeric vector of inferred cell phases in radians.
#' @param r_cell A numeric vector of cell-level posterior rhythmic
#' probabilities.
#' @param strata_by A character string specifying the metadata column used to
#' define strata for gene-level rhythmic fitting.
#' @param B Number of bootstrap iterations. Defaults to \code{100}.
#' @param n_worker Number of workers used for parallel computation. If
#' \code{NULL}, \code{max(1, future::availableCores() - 1)} is used. Defaults to
#' \code{NULL}.
#' @param setseed Random seed. Defaults to \code{1}.
#' @param Verbose Logical; if \code{TRUE}, progress messages are printed.
#'
#' @return A list containing gene-level rhythmic parameter estimates. The main
#' components are \code{param_long}, a long-format parameter table, and
#' \code{param_wide}, a wide-format table for downstream rhythm-state module
#' inference.
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
#' cell_use <- srb_example_data$meta$celltype == "1"
#'
#' res_param <- infer_gene_param(
#'   count = srb_example_data$Count[, cell_use, drop = FALSE],
#'   metadata = srb_example_data$meta[cell_use, , drop = FALSE],
#'   theta_pred = res_phase$theta_pred[cell_use],
#'   r_cell = res_phase$r_cell[cell_use],
#'   strata_by = "condition",
#'   B = 5,
#'   n_worker = 1,
#'   setseed = 1,
#'   Verbose = FALSE
#' )
#'
#' head(res_param$param_long)
#' head(res_param$param_wide)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom stats plogis p.adjust quantile sd median
#' @export
infer_gene_param <- function(
    count,
    metadata,
    theta_pred,
    r_cell,
    strata_by,
    B = 100,
    n_worker = NULL,
    setseed = 1,
    Verbose = TRUE
) {
  t_start <- Sys.time()

  if (!(is.matrix(count) || inherits(count, "Matrix"))) {
    stop("count must be a matrix or Matrix-derived sparse matrix, with genes in rows and cells in columns.")
  }

  if (is.null(rownames(count))) {
    stop("count must have rownames as gene symbols.")
  }

  if (is.null(colnames(count))) {
    colnames(count) <- paste0("Cell_", seq_len(ncol(count)))
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame.")
  }

  if (!is.character(strata_by) || length(strata_by) != 1) {
    stop("strata_by must be a single character string specifying a column name in metadata.")
  }

  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!"cell" %in% colnames(meta)) {
    if (!is.null(rownames(meta))) {
      meta$cell <- rownames(meta)
    } else {
      stop("metadata must contain a 'cell' column or rownames matching colnames(count).")
    }
  }

  meta$cell <- as.character(meta$cell)

  if (anyDuplicated(meta$cell)) {
    stop("metadata contains duplicated cell identifiers.")
  }

  if (!strata_by %in% colnames(meta)) {
    stop("strata_by is not found in metadata.")
  }

  idx <- match(colnames(count), meta$cell)

  if (any(is.na(idx))) {
    stop("Some cells in count are not found in metadata.")
  }

  meta <- meta[idx, , drop = FALSE]

  if (!all(meta$cell == colnames(count))) {
    stop("Cell order mismatch between count and metadata after matching.")
  }

  if (length(theta_pred) != ncol(count)) {
    stop("Length of theta_pred must equal ncol(count).")
  }

  if (length(r_cell) != ncol(count)) {
    stop("Length of r_cell must equal ncol(count).")
  }

  if (!is.numeric(theta_pred)) {
    theta_pred <- as.numeric(theta_pred)
  }

  if (!is.numeric(r_cell)) {
    r_cell <- as.numeric(r_cell)
  }

  stratum <- as.character(meta[[strata_by]])

  if (any(!is.finite(theta_pred))) {
    stop("theta_pred contains NA/Inf. Please provide valid numeric phase values.")
  }

  if (any(!is.finite(r_cell))) {
    stop("r_cell contains NA/Inf. Please provide valid numeric rhythmic probabilities.")
  }

  if (any(is.na(stratum))) {
    stop("The selected strata_by column contains NA. Please provide valid stratum labels.")
  }

  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    stop("B must be a positive integer.")
  }

  B <- as.integer(B)

  if (is.null(n_worker)) {
    n_worker <- max(1L, future::availableCores() - 1L)
  } else {
    if (!is.numeric(n_worker) || length(n_worker) != 1 || n_worker < 1) {
      stop("n_worker must be NULL or a positive integer.")
    }

    n_worker <- as.integer(n_worker)
  }

  set.seed(setseed)

  if (inherits(count, "Matrix")) {
    X <- methods::as(count, "dgCMatrix")
    X@x <- as.double(X@x)
    X@x <- round(X@x)
    X@x[X@x < 0] <- 0
  } else {
    storage.mode(count) <- "double"
    count <- round(count)
    count[count < 0] <- 0
    X <- Matrix::Matrix(count, sparse = TRUE)
    X <- methods::as(X, "dgCMatrix")
  }

  cell <- colnames(X)

  meta_fit <- tibble::tibble(
    cell = cell,
    stratum = stratum,
    theta_pred = as.numeric(theta_pred),
    r_cell = as.numeric(r_cell)
  )

  meta_fit <- as.data.frame(meta_fit)
  rownames(meta_fit) <- meta_fit$cell

  meta_fit$nCount_RNA <- as.numeric(Matrix::colSums(X))
  meta_fit$nFeature_RNA <- as.numeric(Matrix::colSums(X > 0))

  gene_nnz <- Matrix::rowSums(X > 0)

  gene_log_gm <- rep(NA_real_, nrow(X))
  row_log_mean <- tapply(log(X@x), X@i + 1L, mean)
  gene_log_gm[as.integer(names(row_log_mean))] <- as.numeric(row_log_mean)
  names(gene_log_gm) <- rownames(X)

  cell_sf <- vapply(seq_len(ncol(X)), function(k) {
    idx_start <- X@p[k] + 1L
    idx_end <- X@p[k + 1L]

    if (idx_start > idx_end) {
      return(NA_real_)
    }

    row_idx <- X@i[idx_start:idx_end] + 1L
    val_x <- X@x[idx_start:idx_end]

    lgm_k <- gene_log_gm[row_idx]
    ok <- is.finite(lgm_k) & is.finite(val_x) & val_x > 0

    if (!any(ok)) {
      return(NA_real_)
    }

    exp(stats::median(log(val_x[ok]) - lgm_k[ok], na.rm = TRUE))
  }, numeric(1))

  names(cell_sf) <- colnames(X)

  if (all(!is.finite(cell_sf))) {
    cell_sf[] <- 1
  } else {
    cell_sf <- cell_sf / exp(mean(log(cell_sf[is.finite(cell_sf)]), na.rm = TRUE))
    cell_sf[!is.finite(cell_sf)] <- 1
    cell_sf[cell_sf <= 0] <- 1
  }

  meta_fit$size_factor_deseq <- as.numeric(cell_sf[meta_fit$cell])
  meta_fit$log_size_factor <- log(meta_fit$size_factor_deseq)

  meta_fit <- meta_fit %>%
    dplyr::mutate(
      p_cell = pmin(pmax(r_cell, 0.01), 0.99)
    )

  strata_use <- sort(unique(meta_fit$stratum))
  genes_all <- rownames(X)
  genes_fit <- genes_all[gene_nnz > 0]
  n_groups <- length(strata_use)

  if (Verbose) {
    message(
      sprintf(
        "[infer_gene_param] %d strata, %d cells, %d genes, %d worker%s",
        n_groups, ncol(X), nrow(X), n_worker,
        ifelse(n_worker > 1, "s", "")
      )
    )
  }

  options(future.globals.maxSize = Inf)
  on.exit(future::plan(future::sequential), add = TRUE)

  future::plan(
    strategy = future::multisession,
    workers = n_worker
  )

  furrr_opts <- furrr::furrr_options(
    seed = TRUE,
    scheduling = 1,
    packages = c("Matrix", "dplyr", "tibble", "tidyr", "stats", "methods")
  )

  fit_one_gene_group <- function(y_vec, cos_vec, sin_vec, p_vec, log_sf_vec) {
    nnz <- sum(y_vec > 0)
    y_sum <- sum(y_vec)
    feasible <- (nnz >= 2) & (y_sum > 0)

    out <- rep(NA_real_, 14)
    names(out) <- c(
      "Mesor", "Amplitude", "Phase", "Alpha", "Dropout_rate",
      "LRT_stat", "P_value", "Mean",
      "beta1_est", "beta2_est",
      "log_alpha_est", "logit_pi_est", "conv_code", "feasible"
    )

    y_norm <- y_vec / exp(log_sf_vec)
    z <- log(y_norm + 0.1)

    design_mat <- cbind(1, cos_vec, sin_vec)

    fit_lm <- try(
      stats::lm.wfit(
        x = design_mat,
        y = z,
        w = p_vec
      ),
      silent = TRUE
    )

    ok_lm <- !inherits(fit_lm, "try-error")

    coef_init <- rep(0, 3)

    if (ok_lm) {
      coef_init <- fit_lm$coefficients
    }

    coef_init[!is.finite(coef_init)] <- 0

    init_M <- coef_init[1]
    init_beta1 <- coef_init[2]
    init_beta2 <- coef_init[3]

    mu0 <- mean(y_norm)
    v0 <- stats::var(y_norm)
    alpha0 <- pmax((v0 - mu0) / pmax(mu0^2, 1e-8), 1e-3)
    init_log_alpha <- log(alpha0)

    zero_rate <- mean(y_vec == 0)
    init_logit_pi <- stats::qlogis(pmin(pmax(zero_rate * 0.95, 0.01), 0.95))

    nll_full <- function(par) {
      M <- par[1]
      beta1 <- par[2]
      beta2 <- par[3]
      log_alpha <- par[4]
      logit_pi <- par[5]

      alpha <- exp(log_alpha)
      pi_drop <- plogis(logit_pi)

      eta_rhythm <- M + beta1 * cos_vec + beta2 * sin_vec + log_sf_vec
      eta_null <- M + log_sf_vec

      mu_rhythm <- exp(eta_rhythm)
      mu_null <- exp(eta_null)

      prob_zero_rhythm <- pi_drop + (1 - pi_drop) *
        stats::dnbinom(0, mu = mu_rhythm, size = 1 / alpha)
      prob_nonzero_rhythm <- (1 - pi_drop) *
        stats::dnbinom(y_vec, mu = mu_rhythm, size = 1 / alpha)

      lik_rhythm <- prob_nonzero_rhythm
      lik_rhythm[y_vec == 0] <- prob_zero_rhythm[y_vec == 0]

      prob_zero_null <- pi_drop + (1 - pi_drop) *
        stats::dnbinom(0, mu = mu_null, size = 1 / alpha)
      prob_nonzero_null <- (1 - pi_drop) *
        stats::dnbinom(y_vec, mu = mu_null, size = 1 / alpha)

      lik_null <- prob_nonzero_null
      lik_null[y_vec == 0] <- prob_zero_null[y_vec == 0]

      lik_mix <- p_vec * lik_rhythm + (1 - p_vec) * lik_null
      lik_mix[lik_mix < 1e-12] <- 1e-12

      -sum(log(lik_mix))
    }

    nll_null <- function(par) {
      M <- par[1]
      log_alpha <- par[2]
      logit_pi <- par[3]

      alpha <- exp(log_alpha)
      pi_drop <- plogis(logit_pi)

      eta_null <- M + log_sf_vec
      mu_null <- exp(eta_null)

      prob_zero_null <- pi_drop + (1 - pi_drop) *
        stats::dnbinom(0, mu = mu_null, size = 1 / alpha)
      prob_nonzero_null <- (1 - pi_drop) *
        stats::dnbinom(y_vec, mu = mu_null, size = 1 / alpha)

      lik_null <- prob_nonzero_null
      lik_null[y_vec == 0] <- prob_zero_null[y_vec == 0]
      lik_null[lik_null < 1e-12] <- 1e-12

      -sum(log(lik_null))
    }

    fit_full <- try(
      stats::optim(
        par = c(init_M, init_beta1, init_beta2, init_log_alpha, init_logit_pi),
        fn = nll_full,
        method = "L-BFGS-B",
        lower = c(-20, -10, -10, -10, -8),
        upper = c(20, 10, 10, 10, 8),
        control = list(maxit = 150, factr = 1e7)
      ),
      silent = TRUE
    )

    ok_full <- !inherits(fit_full, "try-error")

    null_init <- c(init_M, init_log_alpha, init_logit_pi)

    if (ok_full) {
      null_init <- c(
        fit_full$par[1],
        fit_full$par[4],
        fit_full$par[5]
      )
    }

    fit_null <- try(
      stats::optim(
        par = null_init,
        fn = nll_null,
        method = "L-BFGS-B",
        lower = c(-20, -10, -8),
        upper = c(20, 10, 8),
        control = list(maxit = 120, factr = 1e7)
      ),
      silent = TRUE
    )

    ok_null <- !inherits(fit_null, "try-error")
    ok_both <- feasible & ok_full & ok_null

    if (!ok_both) {
      out["feasible"] <- as.numeric(feasible)
      return(out)
    }

    M_est <- fit_full$par[1]
    beta1_est <- fit_full$par[2]
    beta2_est <- fit_full$par[3]
    log_alpha_est <- fit_full$par[4]
    logit_pi_est <- fit_full$par[5]
    conv_code <- fit_full$convergence

    lrt_stat <- 2 * (fit_null$value - fit_full$value)
    lrt_stat <- pmax(lrt_stat, 0)

    Mesor <- M_est
    Amplitude <- sqrt(beta1_est^2 + beta2_est^2)
    Phase <- atan2(beta2_est, beta1_est) %% (2 * pi)
    Alpha <- exp(log_alpha_est)
    Dropout_rate <- plogis(logit_pi_est)
    P_value <- stats::pchisq(lrt_stat, df = 2, lower.tail = FALSE)
    Mean <- exp(M_est)

    out[] <- c(
      Mesor, Amplitude, Phase, Alpha, Dropout_rate,
      lrt_stat, P_value, Mean,
      beta1_est, beta2_est,
      log_alpha_est, logit_pi_est, conv_code, 1
    )

    out
  }

  stratum_results <- vector("list", length(strata_use))
  names(stratum_results) <- strata_use

  for (grp in strata_use) {
    stratum_start <- Sys.time()

    idx <- which(meta_fit$stratum == grp)
    meta_g <- meta_fit[idx, , drop = FALSE]
    Xg <- X[genes_fit, idx, drop = FALSE]

    if (Verbose) {
      message(
        sprintf(
          "[stratum %d/%d] %s: %d cells",
          match(grp, strata_use), length(strata_use), grp, ncol(Xg)
        )
      )
    }

    theta_vec <- meta_g$theta_pred
    cos_vec <- cos(theta_vec)
    sin_vec <- sin(theta_vec)
    p_vec <- meta_g$p_cell
    log_sf_vec <- meta_g$log_size_factor

    gene_ids <- seq_len(nrow(Xg))

    res_list <- furrr::future_map(
      .x = gene_ids,
      .f = function(j) {
        y_vec <- as.numeric(Xg[j, ])

        fit_one_gene_group(
          y_vec = y_vec,
          cos_vec = cos_vec,
          sin_vec = sin_vec,
          p_vec = p_vec,
          log_sf_vec = log_sf_vec
        )
      },
      .options = furrr_opts,
      .progress = Verbose
    )

    res_mat <- do.call(rbind, res_list)

    res_df <- as.data.frame(res_mat)
    res_df$gene <- rownames(Xg)
    res_df$stratum <- grp
    res_df$FDR <- stats::p.adjust(res_df$P_value, method = "BH")

    res_df <- res_df %>%
      dplyr::select(
        gene, stratum,
        Mesor, Amplitude, Phase, Alpha, Dropout_rate,
        LRT_stat, P_value, FDR, Mean,
        beta1_est, beta2_est, log_alpha_est, logit_pi_est, conv_code
      )

    stratum_results[[grp]] <- res_df

    if (Verbose) {
      message(
        sprintf(
          "[stratum %d/%d] done in %.2f min",
          match(grp, strata_use), length(strata_use),
          as.numeric(difftime(Sys.time(), stratum_start, units = "mins"))
        )
      )
    }
  }

  param_long <- dplyr::bind_rows(stratum_results)

  param_wide <- param_long %>%
    dplyr::select(
      gene, stratum,
      Mesor, Mean, Amplitude, Phase, Alpha, Dropout_rate, FDR
    ) %>%
    tidyr::pivot_wider(
      names_from = stratum,
      values_from = c(Mesor, Mean, Amplitude, Phase, Alpha, Dropout_rate, FDR)
    )

  if (n_groups == 2) {
    grp1 <- strata_use[1]
    grp2 <- strata_use[2]

    param_wide <- param_wide %>%
      dplyr::mutate(
        Delta_Mesor = .data[[paste0("Mesor_", grp2)]] -
          .data[[paste0("Mesor_", grp1)]],
        Delta_Mean = .data[[paste0("Mean_", grp2)]] -
          .data[[paste0("Mean_", grp1)]],
        Delta_Amp = .data[[paste0("Amplitude_", grp2)]] -
          .data[[paste0("Amplitude_", grp1)]],
        Delta_Phase = (
          (.data[[paste0("Phase_", grp2)]] -
             .data[[paste0("Phase_", grp1)]] + pi) %% (2 * pi)
        ) - pi
      )
  }

  gene_rank_df <- param_long %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      min_FDR = if (all(is.na(FDR))) NA_real_ else min(FDR, na.rm = TRUE),
      max_Amp = if (all(is.na(Amplitude))) NA_real_ else max(Amplitude, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(min_FDR), is.finite(max_Amp)) %>%
    dplyr::arrange(min_FDR, dplyr::desc(max_Amp))

  genes_boot <- gene_rank_df %>%
    dplyr::slice_head(n = 100) %>%
    dplyr::pull(gene)

  if (length(genes_boot) == 0) {
    genes_boot <- character(0)
  }

  bootstrap_one_gene_group <- function(
    y0, cos0, sin0, p0, logsf0,
    point_mesor, point_mean, point_amp, point_phase, B
  ) {
    boot_mesor <- rep(NA_real_, B)
    boot_mean <- rep(NA_real_, B)
    boot_amp <- rep(NA_real_, B)
    boot_phase <- rep(NA_real_, B)

    n_cell <- length(y0)

    for (b in seq_len(B)) {
      idb <- sample.int(n_cell, n_cell, replace = TRUE)

      fit_b <- fit_one_gene_group(
        y_vec = y0[idb],
        cos_vec = cos0[idb],
        sin_vec = sin0[idb],
        p_vec = p0[idb],
        log_sf_vec = logsf0[idb]
      )

      boot_mesor[b] <- fit_b["Mesor"]
      boot_mean[b] <- fit_b["Mean"]
      boot_amp[b] <- fit_b["Amplitude"]
      boot_phase[b] <- fit_b["Phase"]
    }

    valid_boot <- is.finite(boot_mesor) &
      is.finite(boot_mean) &
      is.finite(boot_amp) &
      is.finite(boot_phase)

    if (!any(valid_boot)) {
      return(data.frame(
        Mesor_est = point_mesor,
        Mesor_CI_low = NA_real_,
        Mesor_CI_high = NA_real_,
        Mean_est = point_mean,
        Mean_CI_low = NA_real_,
        Mean_CI_high = NA_real_,
        Amplitude_est = point_amp,
        Amp_CI_low = NA_real_,
        Amp_CI_high = NA_real_,
        Phase_est = point_phase,
        Phase_CI_low = NA_real_,
        Phase_CI_high = NA_real_,
        Bootstrap_valid_fraction = 0,
        stringsAsFactors = FALSE
      ))
    }

    mesor_ci <- stats::quantile(
      boot_mesor[valid_boot],
      probs = c(0.025, 0.975),
      na.rm = TRUE
    )

    mean_ci <- stats::quantile(
      boot_mean[valid_boot],
      probs = c(0.025, 0.975),
      na.rm = TRUE
    )

    amp_ci <- stats::quantile(
      boot_amp[valid_boot],
      probs = c(0.025, 0.975),
      na.rm = TRUE
    )

    phase_diff <- ((boot_phase[valid_boot] - point_phase + pi) %% (2 * pi)) - pi

    phase_diff_ci <- stats::quantile(
      phase_diff[is.finite(phase_diff)],
      probs = c(0.025, 0.975),
      na.rm = TRUE
    )

    phase_ci <- (point_phase + phase_diff_ci) %% (2 * pi)

    data.frame(
      Mesor_est = point_mesor,
      Mesor_CI_low = as.numeric(mesor_ci[1]),
      Mesor_CI_high = as.numeric(mesor_ci[2]),
      Mean_est = point_mean,
      Mean_CI_low = as.numeric(mean_ci[1]),
      Mean_CI_high = as.numeric(mean_ci[2]),
      Amplitude_est = point_amp,
      Amp_CI_low = as.numeric(amp_ci[1]),
      Amp_CI_high = as.numeric(amp_ci[2]),
      Phase_est = point_phase,
      Phase_CI_low = as.numeric(phase_ci[1]),
      Phase_CI_high = as.numeric(phase_ci[2]),
      Bootstrap_valid_fraction = mean(valid_boot),
      stringsAsFactors = FALSE
    )
  }

  boot_all <- list()

  if (length(genes_boot) > 0) {
    for (grp in strata_use) {
      boot_start <- Sys.time()

      idx <- which(meta_fit$stratum == grp)
      meta_g <- meta_fit[idx, , drop = FALSE]
      Xg <- X[genes_boot, idx, drop = FALSE]

      if (Verbose) {
        message(
          sprintf(
            "[bootstrap %d/%d] %s: %d genes, B = %d",
            match(grp, strata_use), length(strata_use), grp,
            length(genes_boot), B
          )
        )
      }

      theta0 <- meta_g$theta_pred
      cos0 <- cos(theta0)
      sin0 <- sin(theta0)
      p0 <- meta_g$p_cell
      logsf0 <- meta_g$log_size_factor

      point_g <- stratum_results[[grp]]
      point_g <- point_g[match(genes_boot, point_g$gene), ]

      gene_ids <- seq_along(genes_boot)

      boot_gene_list <- furrr::future_map(
        .x = gene_ids,
        .f = function(j) {
          gene_j <- genes_boot[j]

          y0 <- as.numeric(Xg[j, ])
          point_mesor <- point_g$Mesor[j]
          point_mean <- point_g$Mean[j]
          point_amp <- point_g$Amplitude[j]
          point_phase <- point_g$Phase[j]

          boot_df_j <- bootstrap_one_gene_group(
            y0 = y0,
            cos0 = cos0,
            sin0 = sin0,
            p0 = p0,
            logsf0 = logsf0,
            point_mesor = point_mesor,
            point_mean = point_mean,
            point_amp = point_amp,
            point_phase = point_phase,
            B = B
          )

          boot_df_j$gene <- gene_j
          boot_df_j$stratum <- grp

          boot_df_j %>%
            dplyr::select(
              gene, stratum,
              Mesor_est, Mesor_CI_low, Mesor_CI_high,
              Mean_est, Mean_CI_low, Mean_CI_high,
              Amplitude_est, Amp_CI_low, Amp_CI_high,
              Phase_est, Phase_CI_low, Phase_CI_high,
              Bootstrap_valid_fraction
            )
        },
        .options = furrr_opts,
        .progress = Verbose
      )

      boot_df <- dplyr::bind_rows(boot_gene_list)
      boot_all[[grp]] <- boot_df

      if (Verbose) {
        message(
          sprintf(
            "[bootstrap %d/%d] done in %.2f min",
            match(grp, strata_use), length(strata_use),
            as.numeric(difftime(Sys.time(), boot_start, units = "mins"))
          )
        )
      }
    }

    bootstrap_res <- dplyr::bind_rows(boot_all)
  } else {
    bootstrap_res <- tibble::tibble(
      gene = character(0),
      stratum = character(0),
      Mesor_est = numeric(0),
      Mesor_CI_low = numeric(0),
      Mesor_CI_high = numeric(0),
      Mean_est = numeric(0),
      Mean_CI_low = numeric(0),
      Mean_CI_high = numeric(0),
      Amplitude_est = numeric(0),
      Amp_CI_low = numeric(0),
      Amp_CI_high = numeric(0),
      Phase_est = numeric(0),
      Phase_CI_low = numeric(0),
      Phase_CI_high = numeric(0),
      Bootstrap_valid_fraction = numeric(0)
    )
  }

  if ("stratum" %in% colnames(param_long)) {
    colnames(param_long)[colnames(param_long) == "stratum"] <- strata_by
  }

  if ("stratum" %in% colnames(bootstrap_res)) {
    colnames(bootstrap_res)[colnames(bootstrap_res) == "stratum"] <- strata_by
  }

  if (Verbose) {
    message(
      sprintf(
        "[infer_gene_param] finished in %.2f min",
        as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      )
    )
  }

  out <- list(
    param_long = param_long,
    param_wide = param_wide,
    bootstrap_res = bootstrap_res
  )

  return(out)
}
