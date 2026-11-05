simulate_srb_example_data <- function(
    genes_per_module = 12L,
    cells_per_sample = 60L,
    seed = 42L
) {
  stopifnot(genes_per_module >= 8L, cells_per_sample >= 30L)
  set.seed(seed)

  zt_levels <- c(0, 4, 8, 12, 16, 20)
  design <- expand.grid(
    condition = c("A", "B"),
    ZT = zt_levels,
    replicate = seq_len(3L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design <- design[order(design$condition, design$ZT, design$replicate), ]
  design$sample <- sprintf(
    "%s_ZT%02d_R%d", design$condition, design$ZT, design$replicate
  )
  design$sample_phase_effect <- stats::rnorm(nrow(design), 0, 0.10)

  celltypes <- c("1", "2", "3")
  n_per_type <- rep(cells_per_sample %/% length(celltypes), length(celltypes))
  n_per_type[seq_len(cells_per_sample %% length(celltypes))] <-
    n_per_type[seq_len(cells_per_sample %% length(celltypes))] + 1L
  phase_sd <- c("1" = 0.25, "2" = 0.45, "3" = 0.70)
  noisy_sd <- c("1" = 0.12, "2" = 0.24, "3" = 0.40)

  meta <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) {
    ct <- rep(celltypes, n_per_type)
    theta_center <- 2 * pi * design$ZT[i] / 24 + design$sample_phase_effect[i]
    theta_true <- (theta_center + stats::rnorm(length(ct), 0, phase_sd[ct])) %% (2 * pi)
    theta_noisy <- (theta_true + stats::rnorm(length(ct), 0, noisy_sd[ct])) %% (2 * pi)
    phase_error <- abs(((theta_noisy - theta_true + pi) %% (2 * pi)) - pi)
    r_sim <- pmin(0.99, pmax(0.05, exp(-1.8 * phase_error) * stats::runif(length(ct), 0.85, 1)))

    data.frame(
      cell = sprintf("Cell_%04d_%03d", i, seq_along(ct)),
      sample = design$sample[i],
      replicate = design$replicate[i],
      ZT = design$ZT[i],
      celltype = ct,
      condition = design$condition[i],
      T_true = theta_true / (2 * pi) * 24,
      Theta_true = theta_true,
      Theta_noisy = theta_noisy,
      r_sim = r_sim,
      stringsAsFactors = FALSE
    )
  }))
  rownames(meta) <- meta$cell

  module_gene <- unlist(lapply(0:7, function(k) {
    sprintf("M%dG%02d", k, seq_len(genes_per_module))
  }), use.names = FALSE)
  module_truth <- rep(paste0("M", 0:7), each = genes_per_module)
  phase_base <- rep(seq(0, 2 * pi, length.out = genes_per_module + 1L)[-1L], 8L)

  amp_A <- amp_B <- rep(0, length(module_gene))
  phase_A <- phase_B <- phase_base
  module_id <- as.integer(sub("M([0-7]).*", "\\1", module_truth))
  amp_A[module_id == 1L] <- 0.85; amp_B[module_id == 1L] <- 0.85
  amp_A[module_id == 2L] <- 0;    amp_B[module_id == 2L] <- 1.05
  amp_A[module_id == 3L] <- 1.05; amp_B[module_id == 3L] <- 0
  amp_A[module_id == 4L] <- 0.38; amp_B[module_id == 4L] <- 1.15
  amp_A[module_id == 5L] <- 1.15; amp_B[module_id == 5L] <- 0.38
  amp_A[module_id == 6L] <- 0.90; amp_B[module_id == 6L] <- 0.90
  amp_A[module_id == 7L] <- 0.90; amp_B[module_id == 7L] <- 0.90
  phase_B[module_id == 6L] <- (phase_A[module_id == 6L] - pi / 2) %% (2 * pi)
  phase_B[module_id == 7L] <- (phase_A[module_id == 7L] + pi / 2) %% (2 * pi)

  seed_gene <- c(
    "Arntl", "Clock", "Cry1", "Cry2", "Npas2", "Per1", "Per2", "Per3",
    "Rgs16", "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Nfil3", "Rora",
    "Bhlhe40", "Bhlhe41", "Noct", "Nampt", "Wee1", "Cirbp", "Cpt1a", "Usp2"
  )
  seed_phase <- seq(0, 2 * pi, length.out = length(seed_gene) + 1L)[-1L]

  ovr_gene <- unlist(lapply(celltypes, function(ct) {
    sprintf("OVR%sG%02d", ct, seq_len(6L))
  }), use.names = FALSE)
  ovr_target <- rep(celltypes, each = 6L)
  ovr_phase <- rep(seq(0, 2 * pi, length.out = 7L)[-1L], length(celltypes))

  gene_truth <- rbind(
    data.frame(
      gene = module_gene,
      role = "module",
      module_truth = module_truth,
      target_celltype = NA_character_,
      mesor_log = stats::runif(length(module_gene), log(4), log(9)),
      amplitude_A = amp_A,
      amplitude_B = amp_B,
      phase_A = phase_A,
      phase_B = phase_B,
      stringsAsFactors = FALSE
    ),
    data.frame(
      gene = seed_gene,
      role = "phase_seed",
      module_truth = "seed",
      target_celltype = NA_character_,
      mesor_log = stats::runif(length(seed_gene), log(6), log(12)),
      amplitude_A = stats::runif(length(seed_gene), 0.95, 1.25),
      amplitude_B = stats::runif(length(seed_gene), 0.95, 1.25),
      phase_A = seed_phase,
      phase_B = seed_phase,
      stringsAsFactors = FALSE
    ),
    data.frame(
      gene = ovr_gene,
      role = "ovr_marker",
      module_truth = "ovr",
      target_celltype = ovr_target,
      mesor_log = stats::runif(length(ovr_gene), log(4), log(8)),
      amplitude_A = 1.25,
      amplitude_B = 1.25,
      phase_A = ovr_phase,
      phase_B = ovr_phase,
      stringsAsFactors = FALSE
    )
  )
  gene_truth$nb_size <- stats::runif(nrow(gene_truth), 8, 18)
  gene_truth$sample_mesor_sd <- 0.16
  gene_truth$sample_phase_sd <- 0.10

  n_gene <- nrow(gene_truth)
  n_sample <- nrow(design)
  sample_mesor <- matrix(stats::rnorm(n_gene * n_sample, 0, 0.16), n_gene, n_sample)
  sample_phase <- matrix(stats::rnorm(n_gene * n_sample, 0, 0.10), n_gene, n_sample)
  dimnames(sample_mesor) <- dimnames(sample_phase) <- list(gene_truth$gene, design$sample)

  sample_index <- match(meta$sample, design$sample)
  celltype_index <- match(meta$celltype, celltypes)
  library_factor <- stats::rlnorm(nrow(meta), meanlog = 0, sdlog = 0.22)
  celltype_mesor <- matrix(stats::rnorm(n_gene * length(celltypes), 0, 0.08), n_gene)

  count <- matrix(0L, nrow = n_gene, ncol = nrow(meta))
  for (g in seq_len(n_gene)) {
    condition_B <- meta$condition == "B"
    amp <- ifelse(condition_B, gene_truth$amplitude_B[g], gene_truth$amplitude_A[g])
    phase <- ifelse(condition_B, gene_truth$phase_B[g], gene_truth$phase_A[g])

    if (gene_truth$role[g] == "ovr_marker") {
      amp <- ifelse(meta$celltype == gene_truth$target_celltype[g], amp, 0)
    }

    phase <- phase + sample_phase[g, sample_index]
    log_mu <- gene_truth$mesor_log[g] + sample_mesor[g, sample_index] +
      celltype_mesor[g, celltype_index] + amp * cos(meta$Theta_true - phase) +
      log(library_factor)
    mu <- pmin(exp(log_mu), 500)
    y <- stats::rnbinom(nrow(meta), mu = mu, size = gene_truth$nb_size[g])
    zero_prob <- stats::plogis(0.8 - 1.25 * log1p(mu))
    y[stats::runif(nrow(meta)) < zero_prob] <- 0L
    count[g, ] <- y
  }

  rownames(count) <- gene_truth$gene
  colnames(count) <- meta$cell
  count <- Matrix::Matrix(count, sparse = TRUE)

  sample_effects <- data.frame(
    gene = rep(gene_truth$gene, times = n_sample),
    sample = rep(design$sample, each = n_gene),
    mesor_effect = as.numeric(sample_mesor),
    phase_effect = as.numeric(sample_phase),
    stringsAsFactors = FALSE
  )

  list(
    Count = count,
    meta = meta,
    truth = list(
      gene = gene_truth,
      sample = design,
      sample_effects = sample_effects,
      simulation = list(
        count_model = "negative binomial with mean-dependent zero inflation",
        cells_per_sample = cells_per_sample,
        celltype_phase_sd = phase_sd,
        noisy_phase_sd = noisy_sd,
        seed = seed
      )
    )
  )
}
