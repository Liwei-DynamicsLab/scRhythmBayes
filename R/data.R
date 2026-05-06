#' Example single-cell circadian dataset
#'
#' A simulated example dataset used to demonstrate the scRhythmBayes workflow.
#' The object contains count data and cell-level metadata for testing cell-phase
#' inference, gene-level rhythmic parameter estimation, rhythm-state module
#' analysis, OVR analysis, DTW similarity analysis, and phase-rho summaries.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{Count}{A gene-by-cell count matrix. Rows are genes and columns are cells.}
#'   \item{meta}{A data.frame containing cell-level metadata, including
#'   \code{cell}, \code{ZT}, \code{celltype}, \code{condition},
#'   \code{T_true}, and \code{Theta_true}.}
#' }
#'
#' @source Simulated example data generated for scRhythmBayes package examples.
"srb_example_data"
