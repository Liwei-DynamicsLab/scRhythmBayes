# Synthetic parameter fixture for API/definition regression, not research results.
module_parameter_fixture <- function(groups = c("Control A", "Target B")) {
  set.seed(20261018)
  module <- rep(0:7, each = 12)
  q1 <- ifelse(module %in% c(0,2), .8, .001)
  q2 <- ifelse(module %in% c(0,3), .8, .001)
  a1 <- c(.1,1,.1,1,1,3,1,1)[module+1L] * runif(length(module), .95, 1.05)
  a2 <- c(.1,1,1,.1,3,1,1,1)[module+1L] * runif(length(module), .95, 1.05)
  p1 <- runif(length(module),1,2)
  p2 <- (p1 + c(0,0,0,0,0,0,-pi/3,pi/3)[module+1L] + runif(length(module),-.02,.02)) %% (2*pi)
  w <- data.frame(gene = paste0("gene",seq_along(module)))
  for (i in 1:2) {
    w[[paste0("Mean_",groups[i])]] <- 1
    w[[paste0("FDR_",groups[i])]] <- if (i==1) q1 else q2
    w[[paste0("Amplitude_",groups[i])]] <- if (i==1) a1 else a2
    w[[paste0("Phase_",groups[i])]] <- if (i==1) p1 else p2
  }
  w
}

test_that("reference-to-target effects follow F-panel definitions", {
  w <- module_parameter_fixture()
  a <- infer_rhythm_module(w, group1="Control A", group2="Target B", Verbose=FALSE)
  expect_identical(a$module_result$seed_module, rep(0:7,each=12))
  expect_identical(a$comparison$reference_group, "Control A")
  expect_identical(a$comparison$target_group, "Target B")
  expect_true(all(a$input_diagnostics$complete_parameters))
  expect_setequal(a$module_result$module_map, 0:7)
  expect_true(all(a$module_result$rA[a$module_result$module_map==4] > 0))
  expect_true(all(a$module_result$rA[a$module_result$module_map==5] < 0))
  expect_true(all(a$module_result$delta_phi_rad[a$module_result$module_map==6] < 0))
  expect_true(all(a$module_result$delta_phi_rad[a$module_result$module_map==7] > 0))
  b <- infer_rhythm_module(w, group1="Target B", group2="Control A", Verbose=FALSE)
  invert <- c(0L,1L,3L,2L,5L,4L,7L,6L)
  expect_identical(b$module_result$seed_module, invert[a$module_result$seed_module+1L])
  expect_identical(b$module_result$module_map, invert[a$module_result$module_map+1L])
  # Input column order must not silently choose reference/target groups.
  c <- infer_rhythm_module(w[,rev(names(w))],group1="Control A",group2="Target B",Verbose=FALSE)
  expect_equal(a$module_result,c$module_result)
})

test_that("incomplete original parameters are flagged without dropping genes", {
  w <- module_parameter_fixture()
  w$`Phase_Control A`[1] <- NA_real_
  expect_warning(a <- infer_rhythm_module(w,group1="Control A",group2="Target B",Verbose=FALSE),
                 "incomplete original rhythm parameters")
  expect_equal(nrow(a$module_result),nrow(w))
  expect_true(a$input_diagnostics$missing_group1_parameters[1])
  expect_false(a$input_diagnostics$complete_parameters[1])
  w$`FDR_Control A` <- NA
  expect_warning(b <- infer_rhythm_module(w,group1="Control A",group2="Target B",Verbose=FALSE),
                 "incomplete original rhythm parameters")
  expect_equal(nrow(b$module_result),nrow(w))
  expect_true(all(b$input_diagnostics$missing_group1_parameters))
})

test_that("ambiguous genes, bad groups and nonnumeric parameters fail clearly", {
  w <- module_parameter_fixture()
  expect_error(infer_rhythm_module(w,group1=NA_character_,group2="Target B"),"single character")
  expect_error(infer_rhythm_module(w,group1="",group2="Target B"),"single character")
  expect_error(infer_rhythm_module(w,group1="Control A",group2="Control A"),"must be different")
  expect_error(infer_rhythm_module(w[1,,drop=FALSE],group1="Control A",group2="Target B"),"at least two")
  w$gene[2] <- w$gene[1]
  expect_error(infer_rhythm_module(w,group1="Control A",group2="Target B"),"unique")
  w <- module_parameter_fixture()
  w$`FDR_Control A` <- as.character(w$`FDR_Control A`)
  expect_error(infer_rhythm_module(w,group1="Control A",group2="Target B"),"must be numeric")
})
