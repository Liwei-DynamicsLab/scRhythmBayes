test_that("plot_ovr_curve uses fitted harmonic parameters", {
  result <- data.frame(
    gene = c("GeneA", "GeneB"),
    celltype = "TargetType",
    Mesor_target = log(c(5, 8)),
    Mesor_rest = log(c(4, 6)),
    Delta_Mesor = log(c(5, 8)) - log(c(4, 6)),
    Mean_target = c(5, 8),
    Mean_rest = c(4, 6),
    Delta_Mean = c(1, 2),
    Amplitude_target = c(1, 0.7),
    Amplitude_rest = c(0.1, 0.2),
    Delta_Amp = c(0.9, 0.5),
    Phase_target = c(pi / 2, pi),
    Phase_rest = c(0, pi / 3),
    Delta_Phase = c(pi / 2, 2 * pi / 3),
    FDR_target = c(0.001, 0.002),
    FDR_rest = c(0.4, 0.5),
    ovr_marker = TRUE
  )
  ovr <- list(ovr_result = result, group_by = "celltype")

  p <- plot_ovr_curve(
    ovr,
    genes = c("GeneA", "GeneB"),
    target_group = "TargetType",
    n_grid = 25
  )
  curve <- attr(p, "curve_data")

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(curve), 100L)
  expect_identical(levels(curve$group), c("Rest", "Target"))

  theta <- seq(0, 2 * pi, length.out = 25)
  target_absolute <- exp(log(5) + cos(theta - pi / 2))
  observed <- curve$fitted_mean[
    curve$gene == "GeneA" & curve$group == "Target" & curve$theta == 0
  ]
  expected <- target_absolute[1] / mean(target_absolute)

  expect_equal(observed, expected)
  expect_equal(
    mean(curve$fitted_mean[curve$gene == "GeneA" & curve$group == "Target"]),
    1
  )
  expect_equal(
    mean(curve$fitted_mean[curve$gene == "GeneA" & curve$group == "Rest"]),
    1
  )
  expect_equal(p$labels$y, "Model-fitted relative expression")
})

test_that("single-gene OVR curves reproduce the source figure labels", {
  result <- data.frame(
    gene = "GeneA",
    celltype = "TargetType",
    Mesor_target = log(5),
    Mesor_rest = log(4),
    Amplitude_target = 1,
    Amplitude_rest = 0.1,
    Phase_target = pi / 2,
    Phase_rest = 0,
    ovr_marker = TRUE
  )

  p <- plot_ovr_curve(
    list(
      ovr_result = result,
      group_by = "celltype",
      settings = list(groups = "TargetType")
    ),
    genes = "GeneA",
    target_group = "TargetType"
  )

  expect_equal(p$labels$title, "GeneA")
  expect_equal(p$labels$subtitle, "TargetType")
  expect_length(p$layers, 3L)
})
