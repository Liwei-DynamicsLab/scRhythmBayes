# Plot DTW-based phase-distribution similarity heatmap

This function plots a ComplexHeatmap heatmap from a
`scRhythmDTWSimilarity` object returned by
[`run_dtw_similarity()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/run_dtw_similarity.md).
Rows and columns represent celltype-by-group combinations. Heatmap
values represent DTW-based similarity or distance between circular
cell-phase distributions.

## Usage

``` r
plot_dtw_heatmap(
  dtw_res,
  matrix_type = c("similarity", "distance"),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  group_colors = NULL,
  celltype_colors = NULL,
  heatmap_colors = NULL,
  title = NULL,
  width = grid::unit(9, "cm"),
  height = grid::unit(9, "cm")
)
```

## Arguments

- dtw_res:

  A `scRhythmDTWSimilarity` object returned by
  [`run_dtw_similarity()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/run_dtw_similarity.md).

- matrix_type:

  A character value specifying which matrix to plot. Must be one of
  `"similarity"` or `"distance"`. Defaults to `"similarity"`.

- cluster_rows:

  A logical value indicating whether rows should be clustered. Defaults
  to `TRUE`.

- cluster_columns:

  A logical value indicating whether columns should be clustered.
  Defaults to `TRUE`.

- show_row_names:

  A logical value indicating whether row names should be shown. Defaults
  to `TRUE`.

- show_column_names:

  A logical value indicating whether column names should be shown.
  Defaults to `FALSE`.

- group_colors:

  A named character vector specifying colors for groups. If `NULL`, the
  default scRhythmBayes group palette is used. Defaults to `NULL`.

- celltype_colors:

  A named character vector specifying colors for cell types. If `NULL`,
  the default scRhythmBayes discrete palette is used. Defaults to
  `NULL`.

- heatmap_colors:

  A character vector specifying heatmap colors. If `NULL`, the default
  scRhythmBayes DTW heatmap palette is used. Defaults to `NULL`.

- title:

  A character vector specifying the heatmap title. If `NULL`, a default
  title is used. Defaults to `NULL`.

- width:

  A [`unit`](https://rdrr.io/r/grid/unit.html) object specifying heatmap
  body width. Defaults to `grid::unit(9, "cm")`.

- height:

  A [`unit`](https://rdrr.io/r/grid/unit.html) object specifying heatmap
  body height. Defaults to `grid::unit(9, "cm")`.

## Value

An object of class `scRhythmDTWHeatmap`, containing:

- `heatmap`:

  A ComplexHeatmap object.

- `matrix`:

  The plotted DTW similarity or distance matrix.

- `annotation`:

  Annotation data used for row and column annotations.

- `matrix_type`:

  The matrix type used for plotting.

- `palette`:

  Color palettes used for groups, cell types, and heatmap values.

- `params`:

  Plotting parameters.

## Examples

``` r
if (FALSE) { # \dontrun{
data("srb_example_data")

res_phase <- infer_cell_phase(
  count = srb_example_data$Count,
  metadata = srb_example_data$meta,
  ZT_by = "ZT",
  strata_by = "celltype",
  Species = "mouse",
  setseed = 1,
  Verbose = FALSE
)

dtw_res <- run_dtw_similarity(
  res_phase = res_phase,
  metadata = srb_example_data$meta,
  celltype_by = "celltype",
  group_by = "condition",
  min_cell = 10,
  Verbose = TRUE
)

ht_dtw <- plot_dtw_heatmap(dtw_res = dtw_res)

ht_dtw
plot(ht_dtw)
} # }
```
