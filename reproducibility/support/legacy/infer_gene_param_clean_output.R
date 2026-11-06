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
#' @param sample_by Optional metadata column identifying independent biological
#' samples. NULL preserves the original cell-level branch, including its RNG,
#' output columns, LRT P values and top-100 cell bootstrap.
#' @param sample_B Number of whole-sample pairs bootstrap refits used for
#' confidence intervals of the top 100 ranked genes in sample-aware mode
#' (integer >= 50, default 499). This bootstrap does not determine the
#' sample-aware P value; B is used only by the legacy branch.
#' @param sample_score_B Number of biological-sample Rademacher multiplier draws
#' used for the sample-aware score-test P value (integer >= 199, default 9999).
#' The default provides a Monte Carlo P-value resolution of approximately
#' 1e-4 when all multiplier draws are valid, which is more suitable for
#' transcriptome-scale multiple-testing analyses.
#' @param sample_boot_by Optional metadata column used to stratify the
#' whole-sample pairs bootstrap for confidence intervals, e.g. "ZT". It must be
#' constant within sample within fitting stratum, with at least two independent
#' samples per resampling stratum. It is not used to generate score-test
#' multipliers. Never use theta_pred.
#' @param sample_min Minimum independent samples per fitting stratum (>= 6).
#' This is a computational guard, not a claim that six samples suffice.
#' @param sample_min_valid Minimum converged whole-sample bootstrap fraction for
#' confidence intervals (default 0.8).
#'
#' @details In sample-aware mode, the original cell-level ZINB mixture working
#' likelihood is retained for point estimation. Rhythmic significance instead
#' uses a restricted-null biological-sample cluster multiplier score test for
#' the joint null beta1 = beta2 = 0. Cell-level score contributions are summed
#' within biological samples, nuisance-score contributions are projected out
#' using the restricted-null Hessian, and Rademacher multipliers are applied to
#' the biological-sample score vectors to obtain an empirical null distribution.
#' LRT_stat and P_value_LRT_naive retain the original cell-independence working-
#' likelihood comparison for diagnostics only. Whole-sample pairs bootstrap is
#' run only after ranking, for the top 100 genes, to obtain parameter confidence
#' intervals; it does not determine the sample-aware P value. FDR is BH of the sample-aware P values within each
#' fitting stratum. Invalid score tests return NA without falling back to the
#' naive LRT. Point estimates remain cell-weighted. Phase, responsibility and
#' size factors are held fixed during sample-aware inference. No phase-uncertainty
#' propagation, expression-phase circularity correction, between-stratum paired
#' inference, or module-definition changes are made. The cluster multiplier test
#' is an inferential candidate that should be validated by null simulations for
#' the intended sample size and design before use as a calibrated final test.
#'
#' @return A list containing gene-level rhythmic parameter estimates.
#' \code{param_long} contains the main user-facing results;
#' \code{param_wide} provides the wide-format parameter table for downstream
#' rhythm-state analysis; \code{fit_diagnostics} stores numerical and
#' sample-aware fitting diagnostics separately; \code{bootstrap_res} stores
#' confidence intervals for the top-ranked genes. In sample-aware mode,
#' \code{inference} records the biological-sample inference settings.
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
    Verbose = TRUE,
    sample_by = NULL,
    sample_B = 499L,
    sample_score_B = 9999L,
    sample_boot_by = NULL,
    sample_min = 6L,
    sample_min_valid = 0.8
) {
  t_start <- Sys.time()
  # Local binding makes the uploaded function usable via source(), without attach().
  `%>%` <- magrittr::`%>%`
  sample_mode <- !is.null(sample_by)
  if (sample_mode) {
    if (!is.character(sample_by) || length(sample_by) != 1L ||
        is.na(sample_by) || !nzchar(sample_by)) stop("sample_by must be a column name or NULL.")
    if (!is.numeric(sample_B) || length(sample_B) != 1L ||
        !is.finite(sample_B) || sample_B < 50 || sample_B != floor(sample_B))
      stop("sample_B must be an integer >= 50.")
    if (!is.numeric(sample_score_B) || length(sample_score_B) != 1L ||
        !is.finite(sample_score_B) || sample_score_B < 199 ||
        sample_score_B != floor(sample_score_B))
      stop("sample_score_B must be an integer >= 199.")
    if (!is.numeric(sample_min) || length(sample_min) != 1L ||
        !is.finite(sample_min) || sample_min < 6 || sample_min != floor(sample_min))
      stop("sample_min must be an integer >= 6.")
    if (!is.numeric(sample_min_valid) || length(sample_min_valid) != 1L ||
        !is.finite(sample_min_valid) || sample_min_valid <= 0 || sample_min_valid > 1)
      stop("sample_min_valid must be in (0, 1].")
    if (!is.null(sample_boot_by) && (!is.character(sample_boot_by) ||
        length(sample_boot_by) != 1L || is.na(sample_boot_by)))
      stop("sample_boot_by must be a column name or NULL.")
    sample_B <- as.integer(sample_B)
    sample_score_B <- as.integer(sample_score_B)
  } else if (!is.null(sample_boot_by)) {
    stop("sample_boot_by requires sample_by.")
  }

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
  if (sample_mode) {
    if (!sample_by %in% names(meta)) stop("sample_by is not found in metadata.")
    sample_id <- as.character(meta[[sample_by]])
    if (anyNA(sample_id) || any(!nzchar(trimws(sample_id))))
      stop("sample_by contains missing or empty sample identifiers.")
    boot_stratum <- rep("all", nrow(meta))
    if (!is.null(sample_boot_by)) {
      if (!sample_boot_by %in% names(meta)) stop("sample_boot_by is not found in metadata.")
      boot_stratum <- as.character(meta[[sample_boot_by]])
      if (anyNA(boot_stratum) || any(!nzchar(trimws(boot_stratum))))
        stop("sample_boot_by contains missing or empty labels.")
    }
    for (g in unique(stratum[!is.na(stratum)])) {
      ii <- which(stratum == g)
      blocks <- split(ii, sample_id[ii])
      if (length(blocks) < sample_min)
        stop(sprintf("Stratum '%s' has %d independent samples; need at least %d.",
                     g, length(blocks), sample_min))
      labels <- vapply(blocks, function(i) {
        v <- unique(boot_stratum[i])
        if (length(v) != 1L) stop("sample_boot_by must be constant within each sample and fitting stratum.")
        v
      }, character(1))
      if (any(table(labels) < 2L))
        stop("Each sample_boot_by level needs at least two independent samples per fitting stratum.")
      if (qr(cbind(1, cos(theta_pred[ii]), sin(theta_pred[ii])))$rank < 3L)
        stop(sprintf("Stratum '%s' has a rank-deficient phase design.", g))
    }
  }

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

  if (sample_mode) {
    meta_fit$sample_id <- sample_id
    meta_fit$boot_stratum <- boot_stratum
  }
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
  if (sample_mode && length(genes_fit) == 0L) stop("No genes with positive counts to fit.")

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
    if (sample_mode && qr(cbind(1, cos_vec, sin_vec))$rank < 3L) {
      out <- setNames(rep(NA_real_, 14), c(
        "Mesor", "Amplitude", "Phase", "Alpha", "Dropout_rate", "LRT_stat",
        "P_value", "Mean", "beta1_est", "beta2_est", "log_alpha_est",
        "logit_pi_est", "conv_code", "feasible"))
      out["feasible"] <- 0
      return(out)
    }
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
    # Only the opt-in branch tightens numerical checks. Retry the full model
    # from the nested null to avoid a spurious negative likelihood improvement.
    if (sample_mode && ok_null && is.finite(fit_null$value)) {
      retry <- try(stats::optim(
        c(fit_null$par[1], 0, 0, fit_null$par[2:3]), nll_full,
        method = "L-BFGS-B", lower = c(-20, -10, -10, -10, -8),
        upper = c(20, 10, 10, 10, 8),
        control = list(maxit = 300, factr = 1e7)), silent = TRUE)
      if (!inherits(retry, "try-error") && retry$convergence == 0L &&
          is.finite(retry$value) && (!ok_full || fit_full$convergence != 0L ||
          !is.finite(fit_full$value) || retry$value < fit_full$value)) {
        fit_full <- retry
        ok_full <- TRUE
      }
    }
    ok_both <- feasible & ok_full & ok_null
    if (sample_mode && ok_both) {
      ok_both <- fit_full$convergence == 0L && fit_null$convergence == 0L &&
        all(is.finite(c(fit_full$par, fit_null$par, fit_full$value, fit_null$value))) &&
        fit_full$value <= fit_null$value + 1e-6 &&
        all(abs(fit_full$par[2:3]) < 9.999) && abs(fit_full$par[1]) < 19.999
    }

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

  # Resample biological units with replacement, not their individual cells.
  # Index into seq_along(pool): sample(pool) is unsafe for singleton numeric pools.
  fit_bootstrap_full <- function(y_vec, cos_vec, sin_vec, p_vec, log_sf_vec, init_par) {
    out <- setNames(rep(NA_real_, 9), c("Mesor", "Mean", "Amplitude", "Phase", "beta1_est", "beta2_est", "log_alpha_est", "logit_pi_est", "conv_code"))
    if (sum(y_vec > 0) < 2L || sum(y_vec) <= 0) return(out)
    if (qr(cbind(1, cos_vec, sin_vec))$rank < 3L) return(out)
    if (length(init_par) != 5L || any(!is.finite(init_par))) return(out)
    
    nll_full <- function(par) {
      M <- par[1]
      beta1 <- par[2]
      beta2 <- par[3]
      log_alpha <- par[4]
      logit_pi <- par[5]
      
      alpha <- exp(log_alpha)
      pi_drop <- stats::plogis(logit_pi)
      
      eta_rhythm <- M + beta1 * cos_vec + beta2 * sin_vec + log_sf_vec
      eta_null <- M + log_sf_vec
      
      mu_rhythm <- exp(eta_rhythm)
      mu_null <- exp(eta_null)
      
      prob_zero_rhythm <- pi_drop + (1 - pi_drop) * stats::dnbinom(0, mu = mu_rhythm, size = 1 / alpha)
      prob_nonzero_rhythm <- (1 - pi_drop) * stats::dnbinom(y_vec, mu = mu_rhythm, size = 1 / alpha)
      
      lik_rhythm <- prob_nonzero_rhythm
      lik_rhythm[y_vec == 0] <- prob_zero_rhythm[y_vec == 0]
      
      prob_zero_null <- pi_drop + (1 - pi_drop) * stats::dnbinom(0, mu = mu_null, size = 1 / alpha)
      prob_nonzero_null <- (1 - pi_drop) * stats::dnbinom(y_vec, mu = mu_null, size = 1 / alpha)
      
      lik_null <- prob_nonzero_null
      lik_null[y_vec == 0] <- prob_zero_null[y_vec == 0]
      
      lik_mix <- p_vec * lik_rhythm + (1 - p_vec) * lik_null
      lik_mix[lik_mix < 1e-12] <- 1e-12
      
      -sum(log(lik_mix))
    }
    
    fit_full <- try(
      stats::optim(
        par = init_par,
        fn = nll_full,
        method = "L-BFGS-B",
        lower = c(-20, -10, -10, -10, -8),
        upper = c(20, 10, 10, 10, 8),
        control = list(maxit = 80, factr = 1e7)
      ),
      silent = TRUE
    )
    
    if (inherits(fit_full, "try-error")) return(out)
    if (fit_full$convergence != 0L) return(out)
    if (!is.finite(fit_full$value) || any(!is.finite(fit_full$par))) return(out)
    
    M_est <- fit_full$par[1]
    beta1_est <- fit_full$par[2]
    beta2_est <- fit_full$par[3]
    log_alpha_est <- fit_full$par[4]
    logit_pi_est <- fit_full$par[5]
    
    out[] <- c(
      M_est,
      exp(M_est),
      sqrt(beta1_est^2 + beta2_est^2),
      atan2(beta2_est, beta1_est) %% (2 * pi),
      beta1_est,
      beta2_est,
      log_alpha_est,
      logit_pi_est,
      fit_full$convergence
    )
    
    out
  }
  
  make_sample_draws <- function(sample0, boot0) {
    blocks <- split(seq_along(sample0), sample0)
    labels <- vapply(blocks, function(i) boot0[i[1]], character(1))
    pools <- split(seq_along(blocks), labels)
    lapply(seq_len(sample_B), function(b) {
      chosen <- unlist(lapply(pools, function(pool)
        pool[sample.int(length(pool), length(pool), replace = TRUE)]), use.names = FALSE)
      unlist(blocks[chosen], use.names = FALSE)
    })
  }

  cluster_multiplier_score_test <- function(y, co, si, pr, sf, sample_id,
                                            sample_levels, wild_weights) {
    out <- list(
      P_value = NA_real_, statistic = NA_real_, n_valid = 0L,
      valid_fraction = 0, status = "score_test_failed"
    )

    if (sum(y > 0) < 2L || sum(y) <= 0) {
      out$status <- "score_test_infeasible"
      return(out)
    }

    y_norm <- y / exp(sf)
    init_M <- log(mean(y_norm) + 1e-8)
    mu0 <- mean(y_norm)
    v0 <- stats::var(y_norm)
    alpha0 <- pmax((v0 - mu0) / pmax(mu0^2, 1e-8), 1e-3)
    init_log_alpha <- log(alpha0)
    zero_rate <- mean(y == 0)
    init_logit_pi <- stats::qlogis(pmin(pmax(zero_rate * 0.95, 0.01), 0.95))

    nll_null <- function(par) {
      M <- par[1]
      log_alpha <- par[2]
      logit_pi <- par[3]
      alpha <- exp(log_alpha)
      pi_drop <- stats::plogis(logit_pi)
      mu <- exp(M + sf)
      prob_zero <- pi_drop + (1 - pi_drop) *
        stats::dnbinom(0, mu = mu, size = 1 / alpha)
      prob_nonzero <- (1 - pi_drop) *
        stats::dnbinom(y, mu = mu, size = 1 / alpha)
      lik <- prob_nonzero
      lik[y == 0] <- prob_zero[y == 0]
      lik[lik < 1e-12] <- 1e-12
      -sum(log(lik))
    }

    fit_null <- try(
      stats::optim(
        par = c(init_M, init_log_alpha, init_logit_pi),
        fn = nll_null,
        method = "L-BFGS-B",
        lower = c(-20, -10, -8),
        upper = c(20, 10, 8),
        control = list(maxit = 150, factr = 1e7)
      ),
      silent = TRUE
    )

    if (inherits(fit_null, "try-error") || fit_null$convergence != 0L ||
        !is.finite(fit_null$value) || any(!is.finite(fit_null$par))) {
      out$status <- "null_fit_failed"
      return(out)
    }

    par0 <- c(fit_null$par[1], 0, 0, fit_null$par[2], fit_null$par[3])

    loglik_vec <- function(par) {
      M <- par[1]
      beta1 <- par[2]
      beta2 <- par[3]
      log_alpha <- par[4]
      logit_pi <- par[5]
      alpha <- exp(log_alpha)
      pi_drop <- stats::plogis(logit_pi)

      eta_rhythm <- M + beta1 * co + beta2 * si + sf
      eta_null <- M + sf
      mu_rhythm <- exp(eta_rhythm)
      mu_null <- exp(eta_null)

      prob_zero_rhythm <- pi_drop + (1 - pi_drop) *
        stats::dnbinom(0, mu = mu_rhythm, size = 1 / alpha)
      prob_nonzero_rhythm <- (1 - pi_drop) *
        stats::dnbinom(y, mu = mu_rhythm, size = 1 / alpha)
      lik_rhythm <- prob_nonzero_rhythm
      lik_rhythm[y == 0] <- prob_zero_rhythm[y == 0]

      prob_zero_null <- pi_drop + (1 - pi_drop) *
        stats::dnbinom(0, mu = mu_null, size = 1 / alpha)
      prob_nonzero_null <- (1 - pi_drop) *
        stats::dnbinom(y, mu = mu_null, size = 1 / alpha)
      lik_null <- prob_nonzero_null
      lik_null[y == 0] <- prob_zero_null[y == 0]

      lik_mix <- pr * lik_rhythm + (1 - pr) * lik_null
      lik_mix[lik_mix < 1e-12] <- 1e-12
      log(lik_mix)
    }

    score <- matrix(NA_real_, nrow = length(y), ncol = 5L)
    eps <- 1e-5 * pmax(1, abs(par0))

    for (k in seq_len(5L)) {
      p_plus <- par0
      p_minus <- par0
      p_plus[k] <- p_plus[k] + eps[k]
      p_minus[k] <- p_minus[k] - eps[k]
      ll_plus <- loglik_vec(p_plus)
      ll_minus <- loglik_vec(p_minus)
      score[, k] <- (ll_plus - ll_minus) / (2 * eps[k])
    }

    if (any(!is.finite(score))) {
      out$status <- "nonfinite_score"
      return(out)
    }

    nll_full0 <- function(par) -sum(loglik_vec(par))
    H <- try(stats::optimHess(par0, nll_full0), silent = TRUE)

    if (inherits(H, "try-error") || !is.matrix(H) ||
        !all(dim(H) == c(5L, 5L)) || any(!is.finite(H))) {
      out$status <- "hessian_failed"
      return(out)
    }

    beta_idx <- c(2L, 3L)
    nuisance_idx <- c(1L, 4L, 5L)
    H_nn <- H[nuisance_idx, nuisance_idx, drop = FALSE]
    H_bn <- H[beta_idx, nuisance_idx, drop = FALSE]
    H_nn_inv <- try(solve(H_nn), silent = TRUE)

    if (inherits(H_nn_inv, "try-error") || any(!is.finite(H_nn_inv))) {
      out$status <- "nuisance_hessian_singular"
      return(out)
    }

    projection <- H_bn %*% H_nn_inv
    efficient_cell_score <- score[, beta_idx, drop = FALSE] -
      score[, nuisance_idx, drop = FALSE] %*% t(projection)

    G <- length(sample_levels)
    if (G < sample_min || ncol(wild_weights) != G) {
      out$status <- "sample_mapping_failed"
      return(out)
    }

    psi <- matrix(0, nrow = G, ncol = 2L)
    for (s in seq_len(G)) {
      ii <- which(sample_id == sample_levels[s])
      if (length(ii) == 0L) {
        out$status <- "sample_mapping_failed"
        return(out)
      }
      psi[s, ] <- colSums(efficient_cell_score[ii, , drop = FALSE])
    }

    if (any(!is.finite(psi))) {
      out$status <- "nonfinite_cluster_score"
      return(out)
    }

    U_obs <- colSums(psi)
    # Restricted-null cluster sandwich meat. Do not center the cluster scores:
    # under H0 their expectation is zero, and the observed total score is the
    # quantity being tested. The finite-G multiplier distribution is empirical.
    V <- crossprod(psi) * G / (G - 1)

    if (any(!is.finite(V))) {
      out$status <- "nonfinite_cluster_covariance"
      return(out)
    }

    ev <- try(eigen(V, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
    if (inherits(ev, "try-error") || length(ev) != 2L || any(!is.finite(ev)) ||
        max(ev) <= 0 || min(ev) <= max(ev) * 1e-8) {
      out$status <- "singular_cluster_covariance"
      return(out)
    }

    Vinv <- try(solve(V), silent = TRUE)
    if (inherits(Vinv, "try-error") || any(!is.finite(Vinv))) {
      out$status <- "cluster_covariance_inverse_failed"
      return(out)
    }

    T_obs <- as.numeric(crossprod(U_obs, Vinv %*% U_obs))
    if (!is.finite(T_obs) || T_obs < 0) {
      out$status <- "invalid_observed_statistic"
      return(out)
    }

    U_star <- wild_weights %*% psi
    T_star <- rowSums((U_star %*% Vinv) * U_star)
    valid <- is.finite(T_star) & T_star >= 0
    n_valid <- sum(valid)

    if (n_valid < ceiling(0.8 * nrow(wild_weights))) {
      out$status <- "insufficient_valid_multiplier_draws"
      out$n_valid <- n_valid
      out$valid_fraction <- mean(valid)
      return(out)
    }

    out$P_value <- (1 + sum(T_star[valid] >= T_obs)) / (n_valid + 1)
    out$statistic <- T_obs
    out$n_valid <- n_valid
    out$valid_fraction <- mean(valid)
    out$status <- "ok_cluster_multiplier_score"
    out
  }

  sample_fit <- function(y, co, si, pr, sf, n_samples,
                             sample_id, sample_levels, wild_weights) {
    fit <- fit_one_gene_group(y, co, si, pr, sf)
    naive_p <- unname(fit["P_value"])
    fit["P_value"] <- NA_real_

    score_test <- cluster_multiplier_score_test(
      y = y, co = co, si = si, pr = pr, sf = sf,
      sample_id = sample_id,
      sample_levels = sample_levels,
      wild_weights = wild_weights
    )
    fit["P_value"] <- score_test$P_value

    list(
      fit = fit,
      diagnostics = data.frame(
        P_value_LRT_naive = naive_p,
        Sample_score_stat = score_test$statistic,
        N_samples = n_samples,
        Sample_score_valid = score_test$n_valid,
        Sample_score_valid_fraction = score_test$valid_fraction,
        Sample_status = score_test$status
      )
    )
  }

  stratum_results <- vector("list", length(strata_use))
  names(stratum_results) <- strata_use
  diagnostic_results <- vector("list", length(strata_use))
  names(diagnostic_results) <- strata_use

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

    if (sample_mode) {
      sample_levels <- unique(meta_g$sample_id)
      n_samples <- length(sample_levels)
      set.seed(setseed + 10000L * match(grp, strata_use))
      wild_weights <- matrix(
        sample(c(-1, 1), sample_score_B * n_samples, replace = TRUE),
        nrow = sample_score_B,
        ncol = n_samples
      )
    } else {
      sample_levels <- NULL
      n_samples <- NA_integer_
      wild_weights <- NULL
    }
    gene_ids <- seq_len(nrow(Xg))

    res_list <- furrr::future_map(
      .x = gene_ids,
      .f = function(j) {
        y_vec <- as.numeric(Xg[j, ])

        if (sample_mode) {
          return(sample_fit(
            y = y_vec,
            co = cos_vec,
            si = sin_vec,
            pr = p_vec,
            sf = log_sf_vec,
            n_samples = n_samples,
            sample_id = meta_g$sample_id,
            sample_levels = sample_levels,
            wild_weights = wild_weights
          ))
        }
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

    res_mat <- do.call(rbind, if (sample_mode) lapply(res_list, `[[`, "fit") else res_list)

    res_df <- as.data.frame(res_mat)
    res_df$gene <- rownames(Xg)
    res_df$stratum <- grp
    res_df$FDR <- stats::p.adjust(res_df$P_value, method = "BH", n = nrow(res_df))

    fit_diag <- res_df %>%
      dplyr::select(
        gene, stratum, feasible, conv_code, LRT_stat,
        beta1_est, beta2_est, log_alpha_est, logit_pi_est
      )

    if (sample_mode) {
      sample_diag <- dplyr::bind_rows(lapply(res_list, `[[`, "diagnostics"))
      fit_diag <- dplyr::bind_cols(fit_diag, sample_diag)
    }

    diagnostic_results[[grp]] <- fit_diag

    res_df <- res_df %>%
      dplyr::select(
        gene, stratum,
        Mesor, Mean, Amplitude, Phase, Alpha, Dropout_rate,
        P_value, FDR
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
  fit_diagnostics <- dplyr::bind_rows(diagnostic_results)

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

  bootstrap_sample_gene_group <- function(
    y0, cos0, sin0, p0, logsf0, draws,
    point_mesor, point_mean, point_amp, point_phase,
    point_beta1, point_beta2, point_log_alpha, point_logit_pi
  ) {
    fields <- c("Mesor", "Mean", "Amplitude", "Phase")
    boot <- matrix(NA_real_, length(draws), length(fields), dimnames = list(NULL, fields))
    init_par <- c(point_mesor, point_beta1, point_beta2, point_log_alpha, point_logit_pi)

    if (all(is.finite(init_par))) {
      for (b in seq_along(draws)) {
        i <- draws[[b]]
        z <- fit_bootstrap_full(
          y_vec = y0[i], cos_vec = cos0[i], sin_vec = sin0[i], p_vec = p0[i],
          log_sf_vec = logsf0[i], init_par = init_par
        )
        if (is.finite(z["conv_code"]) && z["conv_code"] == 0L) boot[b, ] <- z[fields]
      }
    }

    valid <- apply(boot, 1, function(z) all(is.finite(z)))
    n_valid <- sum(valid)
    fraction <- mean(valid)
    enough <- n_valid >= max(30L, ceiling(sample_min_valid * length(draws)))

    ci <- function(nm) {
      if (!enough) return(c(NA_real_, NA_real_))
      v <- boot[valid, nm]
      if (nm == "Phase") {
        v <- ((v - point_phase + pi) %% (2 * pi)) - pi
        return((point_phase + as.numeric(stats::quantile(v, c(0.025, 0.975), na.rm = TRUE))) %% (2 * pi))
      }
      as.numeric(stats::quantile(v, c(0.025, 0.975), na.rm = TRUE))
    }

    mc <- ci("Mesor")
    nc <- ci("Mean")
    ac <- ci("Amplitude")
    pc <- ci("Phase")

    data.frame(
      Mesor_est = point_mesor,
      Mesor_CI_low = mc[1], Mesor_CI_high = mc[2],
      Mean_est = point_mean,
      Mean_CI_low = nc[1], Mean_CI_high = nc[2],
      Amplitude_est = point_amp,
      Amp_CI_low = ac[1], Amp_CI_high = ac[2],
      Phase_est = point_phase,
      Phase_CI_low = pc[1], Phase_CI_high = pc[2],
      Bootstrap_valid_fraction = fraction,
      stringsAsFactors = FALSE
    )
  }

  boot_all <- list()

  if (sample_mode && length(genes_boot) > 0) {
    for (grp in strata_use) {
      boot_start <- Sys.time()
      idx <- which(meta_fit$stratum == grp)
      meta_g <- meta_fit[idx, , drop = FALSE]
      Xg <- X[genes_boot, idx, drop = FALSE]

      if (Verbose) {
        message(sprintf(
          "[sample bootstrap %d/%d] %s: %d genes, B = %d",
          match(grp, strata_use), length(strata_use), grp, length(genes_boot), sample_B
        ))
      }

      theta0 <- meta_g$theta_pred
      cos0 <- cos(theta0)
      sin0 <- sin(theta0)
      p0 <- meta_g$p_cell
      logsf0 <- meta_g$log_size_factor
      point_g <- stratum_results[[grp]]
      point_g <- point_g[match(genes_boot, point_g$gene), ]
      set.seed(setseed + 20000L * match(grp, strata_use))
      draws <- make_sample_draws(meta_g$sample_id, meta_g$boot_stratum)

      boot_gene_list <- furrr::future_map(
        .x = seq_along(genes_boot),
        .f = function(j) {
          boot_df_j <- bootstrap_sample_gene_group(
            y0 = as.numeric(Xg[j, ]),
            cos0 = cos0,
            sin0 = sin0,
            p0 = p0,
            logsf0 = logsf0,
            draws = draws,
            point_mesor = point_g$Mesor[j],
            point_mean = point_g$Mean[j],
            point_amp = point_g$Amplitude[j],
            point_phase = point_g$Phase[j],
            point_beta1 = point_g$beta1_est[j],
            point_beta2 = point_g$beta2_est[j],
            point_log_alpha = point_g$log_alpha_est[j],
            point_logit_pi = point_g$logit_pi_est[j]
          )
          boot_df_j$gene <- genes_boot[j]
          boot_df_j$stratum <- grp
          boot_df_j %>% dplyr::select(
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

      boot_all[[grp]] <- dplyr::bind_rows(boot_gene_list)
      if (Verbose) {
        message(sprintf(
          "[sample bootstrap %d/%d] done in %.2f min",
          match(grp, strata_use), length(strata_use),
          as.numeric(difftime(Sys.time(), boot_start, units = "mins"))
        ))
      }
    }
    bootstrap_res <- dplyr::bind_rows(boot_all)
  } else if (!sample_mode && length(genes_boot) > 0) {
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

  if ("stratum" %in% colnames(fit_diagnostics)) {
    colnames(fit_diagnostics)[colnames(fit_diagnostics) == "stratum"] <- strata_by
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
    fit_diagnostics = fit_diagnostics,
    bootstrap_res = bootstrap_res
  )

  if (sample_mode) {
    out$inference <- list(
      mode = "restricted_null_cluster_multiplier_score",
      sample_by = sample_by,
      sample_boot_by = sample_boot_by,
      sample_B = sample_B,
      sample_score_B = sample_score_B,
      nominal_score_p_floor = 1 / (sample_score_B + 1),
      sample_min = sample_min,
      sample_min_valid = sample_min_valid,
      P_value = paste0(
        "restricted-null biological-sample cluster multiplier score test ",
        "with Rademacher weights"
      ),
      LRT = "original cell-independence working-likelihood diagnostic only",
      CI = "whole-sample pairs bootstrap for top 100 ranked genes only",
      fixed_inputs = c("theta_pred", "r_cell", "size_factors"),
      sample_counts = vapply(
        split(meta_fit$sample_id, meta_fit$stratum),
        function(x) length(unique(x)), integer(1)
      )
    )
  }
  return(out)
}
