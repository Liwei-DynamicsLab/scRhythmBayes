#' Infer cell-level circadian phase
#'
#' This function infers cell-level circadian phase from a gene-by-cell count
#' matrix and a cell-level metadata table.
#'
#' @param count A gene-by-cell raw count matrix. \code{count} must be either a
#' base R matrix or a sparse \code{Matrix} object, with genes in rows and cells
#' in columns.
#' @param metadata A data.frame containing cell-level metadata.
#' @param ZT_by A character string specifying the column in \code{metadata}
#' that contains sample-level Zeitgeber time in hours.
#' @param strata_by A character string specifying the column in \code{metadata}
#' used as the grouping variable for cell-level phase inference.
#' @param seedgene An optional character vector of seed genes used for phase
#' inference. If \code{NULL}, built-in default seed genes are used.
#' @param Species Species used to select the built-in default seed-gene set when
#' \code{seedgene = NULL}. One of \code{"auto"}, \code{"human"}, or
#' \code{"mouse"}.
#' @param n_grid Number of discrete phase grid points used to approximate the
#' circular latent phase space.
#' @param n_iter_outer Number of outer iterations for phase inference.
#' @param n_iter_seed Number of seed-gene fitting iterations within each outer
#' iteration.
#' @param Kappa Initial concentration parameter for the von Mises prior.
#' @param p_mix Initial mixture weight between the structured von Mises
#' component and the uniform background component.
#' @param Min_n Minimum number of cells required in a stratum for stratum-level
#' parameter updating.
#' @param setseed Integer random seed used for reproducibility.
#' @param Verbose Logical; if \code{TRUE}, progress messages are printed.
#' @param min_seed Minimum fraction of supplied cells with nonzero counts for
#' a matched seed gene. Defaults to 0.01.
#' @param family Count distribution used for gene fitting and phase likelihood.
#' One of \code{"NB"} (default) or \code{"Poisson"}.
#'
#' @return A tibble with one row per cell. The returned columns include
#' \code{cell}, \code{ZT_input}, \code{stratum}, \code{theta_pred},
#' \code{ZT_pred}, \code{Rbar}, and \code{r_cell}. Attributes seed_genes,
#' seed_diagnostics (all matched genes before filtering), n_seed_candidate,
#' n_seed_matched, n_seed_retained, min_seed and family
#' record the exact seed selection and count distribution.
#'
#' @details
#' The inferred cell phase may depend on the input analysis scope. Running
#' \code{infer_cell_phase()} on a subset of cells is not guaranteed to produce
#' identical results to running it on the full dataset and then subsetting,
#' because seed-gene fitting and stratum-level prior updating are performed on
#' the supplied data matrix.
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
#'   Verbose = TRUE
#' )
#'
#' head(res_phase)
#' }
#'
#' @importFrom magrittr %>%
#' @export
infer_cell_phase <- function(
    count,
    metadata,
    ZT_by,
    strata_by,
    seedgene = NULL,
    Species = c("auto", "human", "mouse"),
    n_grid = 60,
    n_iter_outer = 5,
    n_iter_seed = 3,
    Kappa = 5,
    p_mix = 0.85,
    Min_n = 200,
    setseed = 1,
    Verbose = TRUE,
    min_seed = 0.01,
    family = c("NB", "Poisson")
) {
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

  if (!is.character(ZT_by) || length(ZT_by) != 1) {
    stop("ZT_by must be a single character string specifying a column name in metadata.")
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

  if (!ZT_by %in% colnames(meta)) {
    stop("ZT_by is not found in metadata.")
  }

  if (!strata_by %in% colnames(meta)) {
    stop("strata_by is not found in metadata.")
  }

  idx <- match(colnames(count), meta$cell)

  if (any(is.na(idx))) {
    stop("Some cells in count are not found in metadata.")
  }

  meta_use <- meta[idx, , drop = FALSE]

  if (!all(meta_use$cell == colnames(count))) {
    stop("Cell order mismatch between count and metadata after matching.")
  }

  zt <- meta_use[[ZT_by]]

  if (is.factor(zt)) {
    zt <- as.character(zt)
  }

  if (is.character(zt)) {
    zt_num <- suppressWarnings(as.numeric(zt))
    if (any(is.na(zt_num))) {
      stop("ZT_by column must be numeric or coercible to numeric hours.")
    }
    zt <- zt_num
  }

  if (!is.numeric(zt)) {
    stop("ZT_by column must be numeric.")
  }

  zt <- as.numeric(zt)
  zt <- ((zt %% 24) + 24) %% 24

  stratum <- meta_use[[strata_by]]
  stratum <- as.character(stratum)

  if (any(is.na(stratum))) {
    stop("strata_by column contains NA values.")
  }

  .infer_cell_phase_core(
    Count = count,
    ZT = zt,
    Stratum = stratum,
    seedgene = seedgene,
    Species = Species,
    SeedMode = if (is.null(seedgene)) "Default" else "Custom",
    Mode = "Full",
    n_grid = n_grid,
    n_iter_outer = n_iter_outer,
    n_iter_seed = n_iter_seed,
    Kappa = Kappa,
    p_mix = p_mix,
    Min_n = Min_n,
    setseed = setseed,
    Verbose = Verbose,
    min_seed = min_seed,
    family = match.arg(family)
  )
}


#' Internal core routine for cell-phase inference
#'
#' Internal function used for benchmarking and ablation analysis.
#' Not intended as the primary user-facing interface.
#'
#' @noRd
.infer_cell_phase_core <- function(
    Count,
    ZT,
    Stratum,
    seedgene = NULL,
    Species = c("auto", "human", "mouse"),
    SeedMode = c("Default", "Custom", "Random"),
    Mode = c("Full", "NoPrior", "NoStratum", "NoMix", "NoWeight"),
    n_grid = 60,
    n_iter_outer = 5,
    n_iter_seed = 3,
    Kappa = 5,
    p_mix = 0.85,
    Min_n = 200,
    setseed = 1,
    Verbose = TRUE,
    min_seed = 0.01,
    family = c("NB", "Poisson")
) {
  family <- match.arg(family)
  stopifnot(is.numeric(min_seed), length(min_seed) == 1L,
            is.finite(min_seed), min_seed >= 0, min_seed <= 1)
  Species <- match.arg(Species)
  SeedMode <- match.arg(SeedMode)
  Mode <- match.arg(Mode)

  count <- Count
  zt <- ZT
  stratum <- Stratum

  default_mouse <- c(
    "Acacb", "Acadm", "Acly", "Parp1", "Cyp2d9", "Cyp7b1", "Arntl", "Clock",
    "Cry2", "Npas2", "Per1", "Per2", "Per3", "Acot1", "Acot2", "Acot4",
    "Cyp4a10", "Cyp4a14", "Ak3", "Rgs16", "Ccnd1", "Ddc", "H1f2", "Rell1",
    "Sbk1", "Smarcd1", "Tgfbr1", "Alas1", "Igfbp1", "Raet1c", "Rora", "Sult1e1",
    "Fas", "Rxra", "Cirbp", "Fus", "Hsp90aa1", "Hsph1", "Noct", "Cpt1a",
    "Pklr", "Slc2a2", "Slc37a4", "Ucp2", "Aacs", "Acot3", "Apoa1", "Cyp7a1",
    "Elovl3", "Elovl5", "Elovl6", "Fabp5", "Fads2", "Fasn", "Fgf15", "Insig1",
    "Insig2", "Ldlr", "Lipg", "Lpl", "Nr0b2", "Nr1h2", "Nr1h3", "Pctp",
    "Scap", "Srebf1", "Srebf2", "Apoc3", "Avpr1a", "Fkbp4", "Hspd1", "Tubb5",
    "Dbp", "Hlf", "Nfil3", "Tef", "Bhlhe40", "Por", "Fkbp5", "Usp2",
    "Amd1", "Gck", "Hmgcr", "Lpin1", "Saa4", "Wee1", "Nampt", "Aco2",
    "Akr7a5", "Aldh2", "Bhmt", "Cps1", "Eef1d", "Sord", "Bhlhe41", "Nr1d1",
    "Nr1d2", "Rorc", "G6pc", "Pck1", "Ppara", "Ppard", "Pparg", "Id2"
  )

  default_human <- c(
    "ACACB", "ACADM", "ACLY", "PARP1", "CYP2D6", "CYP7B1", "BMAL1", "CLOCK",
    "CRY2", "NPAS2", "PER1", "PER2", "PER3", "ACOT1", "ACOT2", "ACOT4",
    "CYP4A11", "CYP4A22", "AK3", "RGS16", "CCND1", "DDC", "H1-2", "RELL1",
    "SBK1", "SMARCD1", "TGFBR1", "ALAS1", "IGFBP1", "RAET1E", "RORA", "SULT1E1",
    "FAS", "RXRA", "CIRBP", "FUS", "HSP90AA1", "HSPH1", "NOCT", "CPT1A",
    "PKLR", "SLC2A2", "SLC37A4", "UCP2", "AACS", "APOA1", "CYP7A1", "ELOVL3",
    "ELOVL5", "ELOVL6", "FABP5", "FADS2", "FASN", "FGF19", "INSIG1", "INSIG2",
    "LDLR", "LIPG", "LPL", "NR0B2", "NR1H2", "NR1H3", "PCTP", "SCAP",
    "SREBF1", "SREBF2", "APOC3", "AVPR1A", "FKBP4", "HSPD1", "TUBB", "DBP",
    "HLF", "NFIL3", "TEF", "BHLHE40", "POR", "FKBP5", "USP2", "AMD1",
    "GCK", "HMGCR", "LPIN1", "SAA4", "WEE1", "NAMPT", "ACO2", "AKR7A2",
    "ALDH2", "BHMT", "CPS1", "EEF1D", "SORD", "BHLHE41", "NR1D1", "NR1D2",
    "RORC", "G6PC1", "PCK1", "PPARA", "PPARD", "PPARG", "ID2"
  )

  if (!(is.matrix(count) || inherits(count, "Matrix"))) {
    stop("Count must be a matrix or Matrix-derived sparse matrix, with genes x cells.")
  }

  if (is.null(rownames(count))) {
    stop("Count must have rownames as gene symbols.")
  }

  if (is.null(colnames(count))) {
    colnames(count) <- paste0("Cell_", seq_len(ncol(count)))
  }

  if (length(zt) != ncol(count)) {
    stop("Length of ZT must equal ncol(Count).")
  }

  if (length(stratum) != ncol(count)) {
    stop("Length of Stratum must equal ncol(Count).")
  }

  zt <- as.numeric(zt)
  stratum <- as.character(stratum)

  if (any(is.na(zt))) {
    stop("ZT contains NA. Please provide valid numeric time in hours.")
  }

  if (any(is.na(stratum))) {
    stop("Stratum contains NA. Please provide valid stratum labels.")
  }

  # Canonicalize aliases within each species, retaining original matrix row names.
  seed_aliases_mouse <- c(
    Hist1h1c = "H1f2", Ccrn4l = "Noct", Bmal1 = "Arntl"
  )
  seed_aliases_human <- c(
    ARNTL = "BMAL1", HIST1H1C = "H1-2", CCRN4L = "NOCT", G6PC = "G6PC1"
  )
  canonical_seed <- function(x, species) {
    aliases <- switch(species, mouse = seed_aliases_mouse, human = seed_aliases_human)
    x <- tolower(x)
    idx <- match(x, tolower(names(aliases)))
    x[!is.na(idx)] <- tolower(unname(aliases[idx[!is.na(idx)]]))
    x
  }
  match_seeds <- function(x, species = Species_use) {
    x <- unique(canonical_seed(x, species))
    rn <- rownames(count)
    rn_can <- canonical_seed(rn, species)
    idx <- match(x, tolower(rn))
    missing <- is.na(idx)
    idx[missing] <- match(x[missing], rn_can)
    unique(rn[idx[!is.na(idx)]])
  }

  # Auto remains a match-count heuristic; specify Species for reproducible benchmarks.
  if (Species == "auto") {
    n_mouse <- length(match_seeds(default_mouse, species = "mouse"))
    n_human <- length(match_seeds(default_human, species = "human"))
    Species_use <- if (n_human > n_mouse) "human" else "mouse"
  } else {
    Species_use <- Species
  }

  default_seeds <- switch(
    Species_use,
    human = default_human,
    mouse = default_mouse
  )

  if (SeedMode == "Custom") {
    if (is.null(seedgene)) {
      stop("seedgene must be provided when SeedMode = 'Custom'.")
    }
    seeds_use_input <- unique(as.character(seedgene))
  } else if (SeedMode == "Default") {
    seeds_use_input <- default_seeds
  } else {
    set.seed(setseed)
    n_seed_target <- length(match_seeds(default_seeds))
    if (n_seed_target == 0) {
      stop("No default seed genes found in Count; cannot determine random seed number.")
    }
    gene_pool <- rownames(count)
    if (length(gene_pool) < n_seed_target) {
      stop("Not enough genes in Count to sample random seed genes.")
    }
    seeds_use_input <- sample(gene_pool, size = n_seed_target, replace = FALSE)
  }

  use_prior <- TRUE
  use_stratum <- TRUE
  use_mix <- TRUE
  use_weight <- TRUE

  if (Mode == "NoPrior") use_prior <- FALSE
  if (Mode == "NoStratum") use_stratum <- FALSE
  if (Mode == "NoMix") use_mix <- FALSE
  if (Mode == "NoWeight") use_weight <- FALSE

  cell <- colnames(count)
  lib_size <- Matrix::colSums(count)
  lib_size <- pmax(1, as.numeric(lib_size))
  psi <- (2 * pi / 24) * zt

  meta_g <- tibble::tibble(
    cell = cell,
    ZT_input = zt,
    psi = psi,
    stratum = stratum,
    lib_size = lib_size
  )

  if (!use_stratum) {
    meta_g$stratum <- rep("all", nrow(meta_g))
  }

  X_raw <- count[, meta_g$cell, drop = FALSE]

  if (!all(colnames(X_raw) == meta_g$cell)) {
    stop("Cell order mismatch between Count and internal metadata.")
  }

  seed_matched <- match_seeds(seeds_use_input)

  if (length(seed_matched) == 0) {
    stop("No seed genes found in rownames(Count).")
  }

  seed_detect <- Matrix::rowMeans(X_raw[seed_matched, , drop = FALSE] > 0)
  seed_in <- seed_matched[seed_detect >= min_seed]

  if (length(seed_in) == 0) {
    stop("No matched seed genes pass min_seed detectability threshold.")
  }

  if (Verbose) {
    message(sprintf(
      "[SEED] candidate=%d | matched=%d | retained=%d | min_seed=%.3f",
      length(unique(seeds_use_input)),
      length(seed_matched),
      length(seed_in),
      min_seed
    ))
  }

  Y_seed <- X_raw[seed_in, , drop = FALSE]

  if (!all(colnames(Y_seed) == meta_g$cell)) {
    stop("Cell order mismatch after seed extraction.")
  }

  J <- nrow(Y_seed)
  N <- ncol(Y_seed)

  log_s <- log(pmax(1, meta_g$lib_size))
  exp_s <- exp(log_s)

  levs <- sort(unique(meta_g$stratum))

  theta_grid <- seq(0, 2 * pi, length.out = n_grid + 1)
  theta_grid <- theta_grid[-length(theta_grid)]
  M <- length(theta_grid)
  cth <- cos(theta_grid)
  sth <- sin(theta_grid)

  Yg <- Y_seed

  delta <- stats::setNames(rep(0, length(levs)), levs)
  kappa <- stats::setNames(rep(Kappa, length(levs)), levs)
  p_mix_vec <- stats::setNames(rep(p_mix, length(levs)), levs)

  min_n_ct <- Min_n
  kappa_cap <- 50
  kappa_floor <- 0.3
  p_floor <- 0.05
  p_cap <- 0.995
  a0 <- 2
  b0 <- 2
  ridge_lambda <- 1e-3

  gamma <- 2
  q0 <- 0.70
  r0_floor <- 0.55
  r0_cap <- 0.90

  delta_cell <- delta[meta_g$stratum]
  mu_cell <- (meta_g$psi + delta_cell) %% (2 * pi)
  mu_cell <- mu_cell + (2 * pi) * (mu_cell < 0)

  kappa_cell <- as.numeric(kappa[meta_g$stratum])
  p_cell <- as.numeric(p_mix_vec[meta_g$stratum])

  if (use_prior) {
    cos_mat <- cos(outer(mu_cell, theta_grid, FUN = function(a, b) b - a))
    w_vm <- exp(kappa_cell * cos_mat)
    w_vm <- w_vm / rowSums(w_vm)

    if (use_mix) {
      mix_prob <- p_cell * w_vm + (1 - p_cell) * (1 / M)
    } else {
      mix_prob <- w_vm
    }
  } else {
    mix_prob <- matrix(1 / M, nrow = N, ncol = M)
  }

  logQ <- log(mix_prob)
  logQ <- logQ - apply(logQ, 1, max)
  Q <- exp(logQ)
  Q <- Q / rowSums(Q)

  beta0 <- numeric(J)
  betac <- numeric(J)
  betas <- numeric(J)

  names(beta0) <- rownames(Yg)
  names(betac) <- rownames(Yg)
  names(betas) <- rownames(Yg)

  # Refit unstable NB genes at the Poisson limit, keeping fit and likelihood aligned.
  fit_nb_seed <- function(yj, E_cos, E_sin, X_design, log_s) {
    nb_warning <- FALSE
    fit <- tryCatch(
      withCallingHandlers(
        MASS::glm.nb(yj ~ E_cos + E_sin + offset(log_s)),
        warning = function(w) {
          nb_warning <<- TRUE
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )

    nb_ok <- !is.null(fit) && !nb_warning && isTRUE(fit$converged) &&
      is.null(fit$th.warn) && length(fit$theta) == 1L &&
      is.finite(fit$theta) && fit$theta > 0 &&
      all(is.finite(fit$coefficients))

    if (!nb_ok) {
      fit <- stats::glm.fit(
        x = X_design,
        y = yj,
        family = stats::poisson(),
        offset = log_s
      )
      fit$theta <- Inf
      if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients))) {
        stop("NB seed fitting and its Poisson-limit fallback could not produce a stable fit.")
      }
    }

    fit
  }

  size_seed <- rep(Inf, J)

  phase_loglik <- switch(
    family,
    Poisson = function(V) {
      term1_nm <- t(as.matrix(Matrix::t(V) %*% Yg))
      sumexpV <- colSums(exp(V))
      term2_nm <- outer(exp_s, sumexpV)
      term1_nm - term2_nm
    },
    NB = function(V) {
      loglik_nm <- matrix(0, nrow = N, ncol = M)
      for (j in seq_len(J)) {
        mu_j <- outer(exp_s, exp(V[j, ]))
        loglik_nm <- loglik_nm + matrix(
          stats::dnbinom(as.numeric(Yg[j, ]), mu = mu_j,
                         size = size_seed[j], log = TRUE),
          nrow = N, ncol = M
        )
      }
      loglik_nm
    }
  )

  r_cell <- rep(1, N)

  for (outer_it in seq_len(n_iter_outer)) {

    E_cos <- as.numeric(Q %*% cth)
    E_sin <- as.numeric(Q %*% sth)

    psi_now <- meta_g$psi
    c_shift <- E_cos * cos(psi_now) + E_sin * sin(psi_now)
    s_shift <- E_sin * cos(psi_now) - E_cos * sin(psi_now)

    for (it in seq_len(n_iter_seed)) {

      E_cos <- as.numeric(Q %*% cth)
      E_sin <- as.numeric(Q %*% sth)
      X_design <- cbind(Intercept = 1, cos = E_cos, sin = E_sin)

      for (j in seq_len(J)) {
        yj <- as.numeric(Yg[j, ])

        fit <- suppressWarnings(
          switch(
            family,
            NB = fit_nb_seed(yj, E_cos, E_sin, X_design, log_s),
            Poisson = stats::glm.fit(
              x = X_design,
              y = yj,
              family = stats::poisson(),
              offset = log_s
            )
          )
        )

        size_seed[j] <- switch(family, NB = fit$theta, Poisson = Inf)

        b <- fit$coefficients
        b[is.na(b)] <- 0
        b <- b / (1 + ridge_lambda)

        beta0[j] <- b[1]
        betac[j] <- b[2]
        betas[j] <- b[3]
      }

      V <- beta0 + outer(betac, cth) + outer(betas, sth)
      loglik_nm <- phase_loglik(V)

      if (use_prior) {
        cos_prior <- cos(outer(mu_cell, theta_grid, FUN = function(a, b) b - a))
        w_vm <- exp(kappa_cell * cos_prior)
        w_vm <- w_vm / rowSums(w_vm)

        if (use_mix) {
          mix_prob <- p_cell * w_vm + (1 - p_cell) * (1 / M)
        } else {
          mix_prob <- w_vm
        }

        log_prior <- log(mix_prob)
      } else {
        log_prior <- matrix(log(1 / M), nrow = N, ncol = M)
      }

      logQ_new <- log_prior + loglik_nm
      logQ_new <- logQ_new - apply(logQ_new, 1, max)
      Q <- exp(logQ_new)
      Q <- Q / rowSums(Q)
    }

    V <- beta0 + outer(betac, cth) + outer(betas, sth)
    loglik_nm <- phase_loglik(V)

    if (use_mix && use_prior) {
      cos_mat <- cos(outer(mu_cell, theta_grid, FUN = function(a, b) b - a))
      w_vm <- exp(kappa_cell * cos_mat)
      w_vm <- w_vm / rowSums(w_vm)

      tmp_vm <- loglik_nm + log(w_vm + 1e-300)
      m_vm <- apply(tmp_vm, 1, max)
      logZ_vm <- m_vm + log(rowSums(exp(tmp_vm - m_vm)))

      m_u <- apply(loglik_nm, 1, max)
      logZ_unif <- (m_u + log(rowSums(exp(loglik_nm - m_u)))) - log(M)

      logit_p <- log(p_cell + 1e-12) - log(1 - p_cell + 1e-12)
      z <- logit_p + (logZ_vm - logZ_unif)
      r_cell <- 1 / (1 + exp(-z))
      r_cell <- pmin(0.999, pmax(0.001, r_cell))
    } else {
      r_cell <- rep(1, N)
    }

    if (use_weight) {
      df_r <- tibble::tibble(
        stratum = meta_g$stratum,
        r = r_cell
      )

      r0_tab <- df_r %>%
        dplyr::group_by(stratum) %>%
        dplyr::summarise(
          r0 = as.numeric(stats::quantile(r, probs = q0, na.rm = TRUE, names = FALSE, type = 7)),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          r0 = pmin(r0_cap, pmax(r0_floor, r0))
        )

      r0_cell <- r0_tab$r0[match(meta_g$stratum, r0_tab$stratum)]
      r0_cell <- pmax(1e-6, r0_cell)

      w <- (pmin(1, pmax(0, r_cell / r0_cell)))^gamma

      Q <- Q * w + (1 - w) * (1 / M)
      Q <- Q / rowSums(Q)
    }

    df_m <- tibble::tibble(
      stratum = meta_g$stratum,
      c = c_shift,
      s = s_shift,
      r = r_cell
    )

    ct_n <- meta_g %>%
      dplyr::count(stratum, name = "n_ct")

    agg <- df_m %>%
      dplyr::group_by(stratum) %>%
      dplyr::summarise(
        r_sum = sum(r),
        mc = sum(r * c) / pmax(r_sum, 1e-12),
        ms = sum(r * s) / pmax(r_sum, 1e-12),
        n_all = dplyr::n(),
        r_tot = sum(r),
        .groups = "drop"
      ) %>%
      dplyr::left_join(ct_n, by = "stratum") %>%
      dplyr::mutate(
        Rbar = pmin(0.999999, sqrt(mc^2 + ms^2)),
        delta_new = (atan2(ms, mc) %% (2 * pi)),
        delta_new = delta_new + (2 * pi) * (delta_new < 0),
        kappa_new = dplyr::case_when(
          Rbar < 0.53 ~ 2 * Rbar + Rbar^3 + (5 * Rbar^5) / 6,
          Rbar >= 0.53 & Rbar < 0.85 ~ -0.4 + 1.39 * Rbar + 0.43 / (1 - Rbar),
          Rbar >= 0.85 ~ 1 / (Rbar^3 - 4 * Rbar^2 + 3 * Rbar)
        ),
        kappa_new = pmin(kappa_cap, pmax(kappa_floor, kappa_new)),
        p_new = (r_tot + a0 - 1) / (n_all + a0 + b0 - 2),
        p_new = pmin(p_cap, pmax(p_floor, p_new)),
        delta_new = ifelse(n_ct < min_n_ct, NA_real_, as.numeric(delta_new)),
        kappa_new = ifelse(n_ct < min_n_ct, NA_real_, as.numeric(kappa_new)),
        p_new = ifelse(n_ct < min_n_ct, NA_real_, as.numeric(p_new))
      )

    delta_old <- delta
    kappa_old <- kappa
    p_old <- p_mix_vec

    delta_new_vec <- stats::setNames(agg$delta_new, agg$stratum)
    kappa_new_vec <- stats::setNames(agg$kappa_new, agg$stratum)
    p_new_vec <- stats::setNames(agg$p_new, agg$stratum)

    delta_new_vec[is.na(delta_new_vec)] <- delta_old[names(delta_new_vec)][is.na(delta_new_vec)]
    kappa_new_vec[is.na(kappa_new_vec)] <- kappa_old[names(kappa_new_vec)][is.na(kappa_new_vec)]
    p_new_vec[is.na(p_new_vec)] <- p_old[names(p_new_vec)][is.na(p_new_vec)]

    if (use_prior) {
      delta <- delta_new_vec
      kappa <- kappa_new_vec
      p_mix_vec <- p_new_vec
    }

    delta_cell <- delta[meta_g$stratum]
    mu_cell <- (meta_g$psi + delta_cell) %% (2 * pi)
    mu_cell <- mu_cell + (2 * pi) * (mu_cell < 0)
    kappa_cell <- as.numeric(kappa[meta_g$stratum])
    p_cell <- as.numeric(p_mix_vec[meta_g$stratum])

    E_cos_now <- as.numeric(Q %*% cth)
    E_sin_now <- as.numeric(Q %*% sth)
    theta_hat <- (atan2(E_sin_now, E_cos_now) %% (2 * pi))
    theta_hat <- theta_hat + (2 * pi) * (theta_hat < 0)
    ZT_hat <- (theta_hat * 24 / (2 * pi)) %% 24
    Rbar_now <- sqrt(E_cos_now^2 + E_sin_now^2)

    common <- intersect(names(delta_old), names(delta))
    d_delta <- ((delta[common] - delta_old[common] + pi) %% (2 * pi)) - pi
    d_kappa <- kappa[common] - kappa_old[common]
    d_p <- p_mix_vec[common] - p_old[common]

    if (Verbose) {
      message(sprintf(
        "[OUTER %d/%d] mean(ZT)=%.2f | median(Rbar)=%.3f | mean(p)=%.3f | mean(r)=%.3f | max|d_delta|=%.4f | max|d_kappa|=%.4f | max|d_p|=%.4f",
        outer_it, n_iter_outer,
        mean(ZT_hat), median(Rbar_now), mean(p_mix_vec), mean(r_cell),
        max(abs(d_delta), na.rm = TRUE),
        max(abs(d_kappa), na.rm = TRUE),
        max(abs(d_p), na.rm = TRUE)
      ))
    }
  }

  E_cos_fin <- as.numeric(Q %*% cth)
  E_sin_fin <- as.numeric(Q %*% sth)
  theta_pred <- (atan2(E_sin_fin, E_cos_fin) %% (2 * pi))
  theta_pred <- theta_pred + (2 * pi) * (theta_pred < 0)
  ZT_pred <- (theta_pred * 24 / (2 * pi)) %% 24
  Rbar_pred <- sqrt(E_cos_fin^2 + E_sin_fin^2)

  out <- tibble::tibble(
    cell = meta_g$cell,
    ZT_input = meta_g$ZT_input,
    stratum = meta_g$stratum,
    theta_pred = theta_pred,
    ZT_pred = ZT_pred,
    Rbar = Rbar_pred,
    r_cell = r_cell
  )

  attr(out, "seed_genes") <- seed_in
  attr(out, "seed_diagnostics") <- data.frame(
    gene = names(seed_detect),
    detect_fraction = as.numeric(seed_detect),
    retained = names(seed_detect) %in% seed_in,
    row.names = NULL
  )
  attr(out, "n_seed_candidate") <- length(unique(seeds_use_input))
  attr(out, "n_seed_matched") <- length(seed_matched)
  attr(out, "n_seed_retained") <- length(seed_in)
  attr(out, "min_seed") <- min_seed
  attr(out, "family") <- family
  return(out)
}
