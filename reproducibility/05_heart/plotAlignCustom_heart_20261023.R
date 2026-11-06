# Adapted from the verified legacy plotting function; source file unchanged.
# Only the top padding is increased to separate the title from the matrix.
library(ComplexHeatmap)
library(grid)
plotAlignCustom <- function(alignment,
                            title = "Cardiomyocytes",
                            ref_name = "NCR",
                            query_name = "ACR") {

  costMat <- alignment$localCostMatrix
  costMat <- t(apply(costMat, 1, as.numeric))

  linearInd <- sub2ind(
    nrow(costMat),
    alignment$align[[1]]$index1,
    alignment$align[[1]]$index2
  )

  costMat[linearInd] <- NA
  costMat <- data.frame(costMat, row.names = 1:nrow(costMat))
  colnames(costMat) <- 1:ncol(costMat)

  ht <- ComplexHeatmap::Heatmap(
    as.matrix(costMat),
    name = "Cost",
    col = viridis::plasma(100),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = ref_name,
    row_title = query_name,
    top_annotation = NULL,
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    column_title_side = "bottom",
    row_title_side = "left",
    heatmap_legend_param = list(title = "Dissimilarity"),
    column_names_gp = grid::gpar(fontsize = 0),
    row_names_gp = grid::gpar(fontsize = 0)
  )

  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    padding = grid::unit(c(5, 2, 9, 2), "mm")
  )

  grid::grid.text(
    title,
    x = 0.45,
    y = 0.97,
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )

  invisible(NULL)
}
