test_that("srb_example_data has the sample-aware teaching design", {
  data("srb_example_data", package = "scRhythmBayes")

  count <- srb_example_data$Count
  meta <- srb_example_data$meta
  design <- unique(meta[c("sample", "condition", "ZT", "replicate")])

  expect_identical(colnames(count), meta$cell)
  expect_equal(dim(count), c(138L, 2160L))
  expect_equal(nrow(design), 36L)
  expect_true(all(with(design, table(condition, ZT)) == 3L))
  expect_true(all(table(meta$sample) == 60L))
  expect_identical(sort(unique(meta$celltype)), c("1", "2", "3"))
  expect_true(all(c("Arntl", "Clock", "Cry1", "Cry2", "Per1", "Per2", "Per3") %in%
                    rownames(count)))
  expect_identical(
    sort(unique(srb_example_data$truth$gene$module_truth[1:96])),
    paste0("M", 0:7)
  )
})
