#' Infer rhythm-state modules
#'
#' This function assigns genes to rhythm-state modules based on two-group
#' gene-level rhythmic parameter estimates.
#'
#' @param param_wide A wide-format data.frame returned by
#' \code{infer_gene_param()}.
#' @param compare_by Comparison variable name. Defaults to
#' \code{"condition"}.
#' @param group1 The first group used for comparison.
#' @param group2 The second group used for comparison.
#' @param tau_q FDR threshold for significant rhythmicity. Defaults to
#' \code{0.05}.
#' @param tau_A Minimum rhythmic amplitude threshold. Defaults to \code{0.25}.
#' @param tau_dA Log2 amplitude-ratio threshold. Defaults to \code{1.00}.
#' @param tau_dphi Phase-shift threshold in radians. Defaults to \code{pi / 6}.
#' @param out_dir Output directory for result tables. If \code{NULL}, no files
#' are written. Defaults to \code{NULL}.
#' @param Verbose Logical; if \code{TRUE}, progress messages are printed.
#'
#' @return A list containing rhythm-state module assignment results and module
#' summaries. The main components are \code{module_result} and
#' \code{module_summary}.
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
#' res_module <- infer_rhythm_module(
#'   param_wide = res_param$param_wide,
#'   compare_by = "condition",
#'   group1 = "A",
#'   group2 = "B",
#'   Verbose = FALSE
#' )
#'
#' head(res_module$module_result)
#' }
#'
#' @importFrom magrittr %>%
#' @export
infer_rhythm_module <- function(
    param_wide,
    compare_by = "condition",
    group1,
    group2,
    tau_q = 0.05,
    tau_A = 0.25,
    tau_dA = 1.00,
    tau_dphi = pi / 6,
    out_dir = NULL,
    Verbose = TRUE
) {
  t_start <- Sys.time()

  ## fixed internal parameters
  eps_amp <- 1e-8
  alpha0 <- 1
  kappa0 <- 0.5
  lambda_diag <- 0.10
  eps_diag <- 1e-6

  if (!is.data.frame(param_wide)) {
    stop("param_wide must be a data.frame or tibble.")
  }

  if (!"gene" %in% colnames(param_wide)) {
    stop("param_wide must contain a 'gene' column.")
  }

  if (!is.character(compare_by) || length(compare_by) != 1) {
    stop("compare_by must be a single character string.")
  }

  if (!is.character(group1) || length(group1) != 1) {
    stop("group1 must be a single character string.")
  }

  if (!is.character(group2) || length(group2) != 1) {
    stop("group2 must be a single character string.")
  }

  if (identical(group1, group2)) {
    stop("group1 and group2 must be different.")
  }

  if (!is.numeric(tau_q) || length(tau_q) != 1 || tau_q <= 0 || tau_q >= 1) {
    stop("tau_q must be a single numeric value between 0 and 1.")
  }

  if (!is.numeric(tau_A) || length(tau_A) != 1 || tau_A < 0) {
    stop("tau_A must be a single non-negative numeric value.")
  }

  if (!is.numeric(tau_dA) || length(tau_dA) != 1 || tau_dA <= 0) {
    stop("tau_dA must be a single positive numeric value.")
  }

  if (!is.numeric(tau_dphi) || length(tau_dphi) != 1 || tau_dphi <= 0) {
    stop("tau_dphi must be a single positive numeric value.")
  }

  grp1 <- group1
  grp2 <- group2

  mean1_col <- paste0("Mean_", grp1)
  mean2_col <- paste0("Mean_", grp2)
  fdr1_col <- paste0("FDR_", grp1)
  fdr2_col <- paste0("FDR_", grp2)
  amp1_col <- paste0("Amplitude_", grp1)
  amp2_col <- paste0("Amplitude_", grp2)
  phase1_col <- paste0("Phase_", grp1)
  phase2_col <- paste0("Phase_", grp2)

  required_cols <- c(
    mean1_col, mean2_col,
    fdr1_col, fdr2_col,
    amp1_col, amp2_col,
    phase1_col, phase2_col
  )

  miss_cols <- setdiff(required_cols, colnames(param_wide))

  if (length(miss_cols) > 0) {
    stop(
      "param_wide is missing required columns: ",
      paste(miss_cols, collapse = ", ")
    )
  }

  if (Verbose) {
    message(sprintf(
      "[infer_rhythm_module] %s comparison: %s vs %s",
      compare_by, grp1, grp2
    ))
  }

  if (!is.null(out_dir)) {
    dir.create(
      file.path(out_dir, "module_assignment"),
      showWarnings = FALSE,
      recursive = TRUE
    )
  }

  wrap_pi <- function(x) {
    ((x + pi) %% (2 * pi)) - pi
  }

  softmax_row <- function(log_mat) {
    row_max <- apply(log_mat, 1, max)
    exp_mat <- exp(log_mat - row_max)
    exp_mat / rowSums(exp_mat)
  }

  top2_margin <- function(prob_mat) {
    t(apply(prob_mat, 1, function(z) sort(z, decreasing = TRUE)[1:2])) |>
      as.data.frame() |>
      stats::setNames(c("p1", "p2")) |>
      dplyr::mutate(margin = p1 - p2) |>
      dplyr::pull(margin)
  }

  module_df <- param_wide %>%
    dplyr::transmute(
      gene = gene,
      !!fdr1_col := .data[[fdr1_col]],
      !!fdr2_col := .data[[fdr2_col]],
      !!amp1_col := .data[[amp1_col]],
      !!amp2_col := .data[[amp2_col]],
      !!phase1_col := .data[[phase1_col]],
      !!phase2_col := .data[[phase2_col]]
    ) %>%
    dplyr::mutate(
      !!fdr1_col := pmin(pmax(as.numeric(.data[[fdr1_col]]), 1e-300), 1),
      !!fdr2_col := pmin(pmax(as.numeric(.data[[fdr2_col]]), 1e-300), 1),
      !!amp1_col := as.numeric(.data[[amp1_col]]),
      !!amp2_col := as.numeric(.data[[amp2_col]]),
      !!phase1_col := as.numeric(.data[[phase1_col]]),
      !!phase2_col := as.numeric(.data[[phase2_col]]),
      delta_phi_rad = wrap_pi(.data[[phase2_col]] - .data[[phase1_col]]),
      u = cos(delta_phi_rad),
      v = sin(delta_phi_rad),
      rA = log2((.data[[amp2_col]] + eps_amp) / (.data[[amp1_col]] + eps_amp))
    ) %>%
    dplyr::mutate(
      score_amp = abs(rA) / tau_dA,
      score_phi = abs(delta_phi_rad) / tau_dphi,
      seed_module = dplyr::case_when(
        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] < tau_q &
          .data[[amp1_col]] < tau_A & .data[[amp2_col]] < tau_A ~ 0L,

        .data[[fdr1_col]] >= tau_q & .data[[fdr2_col]] >= tau_q ~ 0L,

        .data[[fdr1_col]] >= tau_q & .data[[fdr2_col]] < tau_q ~ 2L,
        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] >= tau_q ~ 3L,

        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] < tau_q &
          rA > tau_dA &
          abs(delta_phi_rad) < tau_dphi ~ 4L,

        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] < tau_q &
          rA < -tau_dA &
          abs(delta_phi_rad) < tau_dphi ~ 5L,

        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] < tau_q &
          abs(rA) < tau_dA &
          delta_phi_rad < -tau_dphi ~ 6L,

        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] < tau_q &
          abs(rA) < tau_dA &
          delta_phi_rad > tau_dphi ~ 7L,

        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] < tau_q &
          abs(rA) >= tau_dA &
          abs(delta_phi_rad) >= tau_dphi &
          score_amp >= score_phi &
          rA > 0 ~ 4L,

        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] < tau_q &
          abs(rA) >= tau_dA &
          abs(delta_phi_rad) >= tau_dphi &
          score_amp >= score_phi &
          rA < 0 ~ 5L,

        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] < tau_q &
          abs(rA) >= tau_dA &
          abs(delta_phi_rad) >= tau_dphi &
          score_amp < score_phi &
          delta_phi_rad < 0 ~ 6L,

        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] < tau_q &
          abs(rA) >= tau_dA &
          abs(delta_phi_rad) >= tau_dphi &
          score_amp < score_phi &
          delta_phi_rad > 0 ~ 7L,

        .data[[fdr1_col]] < tau_q & .data[[fdr2_col]] < tau_q &
          abs(rA) < tau_dA &
          abs(delta_phi_rad) < tau_dphi ~ 1L,

        TRUE ~ NA_integer_
      )
    )

  feat_cols <- c(
    fdr1_col,
    fdr2_col,
    amp1_col,
    amp2_col,
    "rA",
    "u",
    "v"
  )

  X_raw <- module_df %>%
    dplyr::select(dplyr::all_of(feat_cols)) %>%
    as.matrix()

  X_center <- colMeans(X_raw, na.rm = TRUE)
  X_scale <- apply(X_raw, 2, stats::sd, na.rm = TRUE)
  X_scale <- pmax(X_scale, 1e-8)

  X <- scale(X_raw, center = X_center, scale = X_scale) %>%
    as.matrix()

  X[!is.finite(X)] <- 0

  m0 <- colMeans(X)
  S0 <- diag(
    apply(X, 2, stats::var),
    nrow = ncol(X),
    ncol = ncol(X)
  )

  S0 <- S0 + diag(eps_diag, ncol(X))

  module_ids <- 0:7

  fit_one_module <- function(k) {
    idx <- which(module_df$seed_module == k)
    Xk <- X[idx, , drop = FALSE]
    nk <- nrow(Xk)
    p <- ncol(X)

    xbar_raw <- colSums(Xk) / pmax(nk, 1)
    use_seed <- as.numeric(nk > 0)
    xbar_k <- use_seed * xbar_raw + (1 - use_seed) * m0

    XCk <- sweep(Xk, 2, xbar_k, "-")
    Wk <- crossprod(XCk)

    kappa_star <- kappa0 + nk
    m_star <- (kappa0 * m0 + nk * xbar_k) / kappa_star

    S_star <- S0 +
      Wk +
      (kappa0 * nk / kappa_star) * tcrossprod(xbar_k - m0)

    Sigma_hat <- S_star / pmax(nk + kappa0 + p + 2, 1)

    Sigma_hat <- (1 - lambda_diag) * Sigma_hat +
      lambda_diag * diag(diag(Sigma_hat), nrow = p, ncol = p) +
      diag(eps_diag, p)

    tibble::tibble(
      module = k,
      n_seed = nk,
      mu = list(m_star),
      Sigma = list(Sigma_hat)
    )
  }

  if (Verbose) {
    message("[infer_rhythm_module] fitting module distributions")
  }

  fit_tbl <- purrr::map_dfr(module_ids, fit_one_module)

  pi_tbl <- fit_tbl %>%
    dplyr::transmute(
      module,
      n_seed,
      pi_k = (n_seed + alpha0) / sum(n_seed + alpha0)
    )

  fit_tbl <- fit_tbl %>%
    dplyr::left_join(pi_tbl, by = c("module", "n_seed"))

  score_mat <- purrr::map_dfc(module_ids, function(k) {
    fit_k <- fit_tbl %>%
      dplyr::filter(module == k)

    mu_k <- fit_k$mu[[1]]
    Sigma_k <- fit_k$Sigma[[1]]
    pi_k <- fit_k$pi_k[[1]]

    tibble::tibble(
      !!paste0("post_M", k) :=
        log(pi_k) +
        mvtnorm::dmvnorm(
          X,
          mean = mu_k,
          sigma = Sigma_k,
          log = TRUE
        )
    )
  }) %>%
    as.matrix()

  post_mat <- softmax_row(score_mat)
  colnames(post_mat) <- paste0("post_M", module_ids)

  module_map <- max.col(post_mat, ties.method = "first") - 1L
  post_max <- apply(post_mat, 1, max)
  margin <- top2_margin(post_mat)

  module_key <- tibble::tibble(
    module_id = 0:7,
    module_name = c(
      "Non-rhythmic",
      "Stable Rhythm",
      "Rhythm Gain",
      "Rhythm Loss",
      "Amplitude-enhanced rhythmic",
      "Amplitude-attenuated rhythmic",
      "Phase-advanced rhythmic",
      "Phase-delayed rhythmic"
    )
  )

  module_result <- module_df %>%
    dplyr::bind_cols(as.data.frame(post_mat)) %>%
    dplyr::mutate(
      compare_by = compare_by,
      group1 = grp1,
      group2 = grp2,
      module_map = module_map,
      post_max = post_max,
      margin = margin
    ) %>%
    dplyr::left_join(
      module_key,
      by = c("module_map" = "module_id")
    ) %>%
    dplyr::mutate(
      !!amp1_col := dplyr::case_when(
        module_map %in% c(0L, 2L) ~ NA_real_,
        TRUE ~ .data[[amp1_col]]
      ),
      !!amp2_col := dplyr::case_when(
        module_map %in% c(0L, 3L) ~ NA_real_,
        TRUE ~ .data[[amp2_col]]
      ),
      !!phase1_col := dplyr::case_when(
        module_map %in% c(0L, 2L) ~ NA_real_,
        TRUE ~ .data[[phase1_col]]
      ),
      !!phase2_col := dplyr::case_when(
        module_map %in% c(0L, 3L) ~ NA_real_,
        TRUE ~ .data[[phase2_col]]
      )
    ) %>%
    dplyr::mutate(
      delta_phi_rad = dplyr::case_when(
        module_map %in% c(0L, 2L, 3L) ~ NA_real_,
        TRUE ~ delta_phi_rad
      ),
      rA = dplyr::case_when(
        module_map %in% c(0L, 2L, 3L) ~ NA_real_,
        TRUE ~ rA
      ),
      u = dplyr::case_when(
        module_map %in% c(0L, 2L, 3L) ~ NA_real_,
        TRUE ~ u
      ),
      v = dplyr::case_when(
        module_map %in% c(0L, 2L, 3L) ~ NA_real_,
        TRUE ~ v
      )
    ) %>%
    dplyr::select(
      gene,
      compare_by,
      group1,
      group2,
      dplyr::all_of(fdr1_col),
      dplyr::all_of(fdr2_col),
      dplyr::all_of(amp1_col),
      dplyr::all_of(amp2_col),
      dplyr::all_of(phase1_col),
      dplyr::all_of(phase2_col),
      delta_phi_rad,
      rA,
      u,
      v,
      seed_module,
      module_map,
      module_name,
      dplyr::starts_with("post_M"),
      post_max,
      margin
    )

  module_summary <- module_result %>%
    dplyr::count(
      compare_by,
      group1,
      group2,
      module_map,
      module_name,
      name = "n_gene"
    ) %>%
    dplyr::mutate(
      freq = n_gene / sum(n_gene)
    )

  if (!is.null(out_dir)) {
    readr::write_tsv(
      module_result,
      file.path(
        out_dir,
        "module_assignment",
        "module_assignment_results.tsv"
      )
    )

    readr::write_tsv(
      module_summary,
      file.path(
        out_dir,
        "module_assignment",
        "module_assignment_summary.tsv"
      )
    )
  }

  if (Verbose) {
    message(
      sprintf(
        "[infer_rhythm_module] finished in %.2f min",
        as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      )
    )
  }

  list(
    module_result = module_result,
    module_summary = module_summary
  )
}
