# Deterministic synthetic fixture for plotting regression only, not biological data.
heatmap_direction_fixture <- function(direction = c("PV", "CV")) {
  set.seed(20261018)
  theta <- rep(seq(0, 2 * pi, length.out = 49)[-49], 2)
  genes <- paste0("gene", seq_len(32))
  cells <- paste0("cell", seq_along(theta))
  count <- matrix(rpois(32 * length(theta), 20), nrow = 32,
                  dimnames = list(genes, cells))
  count <- Matrix::Matrix(count, sparse = TRUE)
  metadata <- data.frame(celltype = rep(c("CV", "PV"), each = 48), row.names = cells)
  res_phase <- data.frame(cell = cells, theta_pred = theta)
  gene_meta <- data.frame(gene = genes, module_map = rep(0:7, each = 4),
                          group1 = direction[1], group2 = direction[2])
  for (g in c("PV", "CV")) {
    gene_meta[[paste0("Phase_", g)]] <- rep(c(5, 1, 4, 2), 8)
    gene_meta[[paste0("Amplitude_", g)]] <- 1
    gene_meta[[paste0("FDR_", g)]] <- 0.01
  }
  gene_meta[gene_meta$module_map == 2, paste0("Phase_", direction[1])] <- NA_real_
  gene_meta[gene_meta$module_map == 3, paste0("Phase_", direction[2])] <- NA_real_
  gene_meta[gene_meta$module_map == 5, paste0("Phase_", direction[2])] <- NA_real_
  list(count = count, res_phase = res_phase, gene_meta = gene_meta, metadata = metadata)
}

build_direction_fixture <- function(x, display = c("CV", "PV"), ...) {
  builder <- getOption("scRhythmBayes.heatmap_test_builder", plot_module_heatmap)
  do.call(builder, c(x, list(group_by = "celltype", group1 = display[1],
                            group2 = display[2]), list(...)))
}

test_that("all module references and colors follow inference, not display order", {
  for (direction in list(c("PV", "CV"), c("CV", "PV"))) {
    x <- heatmap_direction_fixture(direction)
    a <- build_direction_fixture(x, module_ids = 0:7)
    b <- build_direction_fixture(x, display = c("PV", "CV"), module_ids = 0:7)
    expect_identical(a$summary$ref_group, direction[c(1, 1, 2, 1, 2, 1, 1, 1)])
    expect_identical(b$summary$ref_group, direction[c(1, 1, 2, 1, 2, 1, 1, 1)])
    expect_identical(a$palette$group, b$palette$group)
    expect_identical(a$group_levels, c("CV", "PV"))
    expect_identical(b$group_levels, c("PV", "CV"))
    for (id in paste0("module", 0:7)) {
      aa <- a$module_list[[id]]
      bb <- b$module_list[[id]]
      expect_equal(aa$mat_left_z, bb$mat_right_z, tolerance = 1e-10)
      expect_equal(aa$mat_right_z, bb$mat_left_z, tolerance = 1e-10)
      expect_equal(aa$phase_order$phase_wrap, c(1, 2, 4, 5))
      expect_length(aa$gene_keep, 4)
    }
    phase_anno <- as.character(a$module_list$module1$phase_annotation_left)
    expect_identical(phase_anno[c(5,15,25,35,45,55)], as.character(c(2,6,10,14,18,22)))
  }
})

test_that("legacy results require explicit comparison direction for M2-M7", {
  x <- heatmap_direction_fixture()
  x$gene_meta$group1 <- x$gene_meta$group2 <- NULL
  expect_error(build_direction_fixture(x, module_ids = 2:5), "need comparison direction")
  expect_error(build_direction_fixture(x, module_ids = 5), "need comparison direction")
  expect_error(build_direction_fixture(x, module_ids = 6:7), "need comparison direction")
  a <- build_direction_fixture(x, module_ids = 2:5, compare_groups = c("PV", "CV"))
  expect_identical(a$summary$ref_group, c("CV", "PV", "CV", "PV"))
})

test_that("unknown phases are retained last and complete failures are explicit", {
  x <- heatmap_direction_fixture()
  x$gene_meta$Phase_CV[9] <- NA_real_
  expect_warning(a <- build_direction_fixture(x, module_ids = 2), "retained last")
  expect_identical(a$summary$n_phase_ordered, 3L)
  expect_identical(a$summary$n_phase_missing, 1L)
  expect_identical(tail(a$module_list$module2$gene_keep, 1), "gene9")
  expect_length(a$module_list$module2$gene_keep, 4)
  x$gene_meta$Phase_CV[x$gene_meta$module_map == 2] <- NA_real_
  expect_error(build_direction_fixture(x, module_ids = 2), "no finite phases")
  x$gene_meta$Phase_CV[x$gene_meta$module_map == 0] <- NA_real_
  x$gene_meta$Phase_PV[x$gene_meta$module_map == 0] <- NA_real_
  x$gene_meta <- x$gene_meta[x$gene_meta$module_map != 7, , drop = FALSE]
  a <- build_direction_fixture(x, module_ids = c(0, 7))
  expect_identical(a$summary$status, c("ok", "no_genes"))
  expect_identical(a$summary$n_phase_missing[1], 4L)
})

test_that("non-liver labels and shuffled cell inputs are supported", {
  x <- heatmap_direction_fixture(c("CV", "PV"))
  names(x$gene_meta) <- gsub("_CV$", "_Control A", names(x$gene_meta))
  names(x$gene_meta) <- gsub("_PV$", "_Target B", names(x$gene_meta))
  x$gene_meta$group1 <- "Control A"
  x$gene_meta$group2 <- "Target B"
  x$metadata$celltype <- ifelse(x$metadata$celltype == "CV", "Control A", "Target B")
  a <- build_direction_fixture(x, display = c("Target B", "Control A"), module_ids = 2:3)
  x$metadata <- x$metadata[nrow(x$metadata):1, , drop = FALSE]
  x$res_phase <- x$res_phase[nrow(x$res_phase):1, , drop = FALSE]
  b <- build_direction_fixture(x, display = c("Target B", "Control A"), module_ids = 2:3)
  expect_identical(a$summary$ref_group, c("Target B", "Control A"))
  expect_equal(a$module_list$module2$mat_left_z, b$module_list$module2$mat_left_z)
  expect_equal(a$module_list$module3$mat_right_z, b$module_list$module3$mat_right_z)
})

test_that("mixed, mismatched and conflicting comparison directions are rejected", {
  x <- heatmap_direction_fixture()
  expect_error(build_direction_fixture(x, compare_groups = c("CV", "PV")), "conflicts")
  x$gene_meta$group1[1] <- "CV"
  expect_error(build_direction_fixture(x), "one complete comparison")
  x <- heatmap_direction_fixture(c("A", "B"))
  expect_error(build_direction_fixture(x), "matching the plotted groups")
  x <- heatmap_direction_fixture()
  x$gene_meta$group2 <- NULL
  expect_error(build_direction_fixture(x), "both group1 and group2")
})

test_that("duplicate genes and nonfinite cell phases do not silently misalign", {
  x <- heatmap_direction_fixture()
  x$gene_meta$gene[2] <- x$gene_meta$gene[1]
  expect_error(build_direction_fixture(x), "unique")
  x <- heatmap_direction_fixture()
  x$res_phase$theta_pred[1] <- NA_real_
  expect_error(build_direction_fixture(x), "finite numeric radians")
  x <- heatmap_direction_fixture()
  expect_error(build_direction_fixture(x, module_ids=9))
})
