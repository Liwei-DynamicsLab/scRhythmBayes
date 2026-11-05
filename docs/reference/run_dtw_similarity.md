# Run DTW-based phase-distribution similarity analysis

This function compares cell-phase distributions across celltype-by-group
combinations using dynamic time warping (DTW). For each combination,
inferred cell phases are converted into a circular histogram. Pairwise
DTW distances between histograms are then converted to a 0-to-1
similarity matrix.

## Usage

``` r
run_dtw_similarity(
  res_phase,
  metadata,
  theta_by = "theta_pred",
  celltype_by = "celltype",
  group_by = "condition",
  group_levels = NULL,
  celltype_levels = NULL,
  n_bin = 24,
  min_cell = 10,
  normalize_distance = TRUE,
  Verbose = FALSE
)
```

## Arguments

- res_phase:

  A data.frame returned by
  [`infer_cell_phase()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_cell_phase.md).
  It must contain `"cell"` and the column specified by `theta_by`.

- metadata:

  A data.frame containing cell-level metadata.

- theta_by:

  A character vector specifying the column name in `res_phase` that
  contains inferred cell phase in radians. Defaults to `"theta_pred"`.

- celltype_by:

  A character vector specifying the column name in `metadata` that
  defines cell types. Defaults to `"celltype"`.

- group_by:

  A character vector specifying the column name in `metadata` that
  defines groups or conditions. Defaults to `"condition"`.

- group_levels:

  A character vector specifying group order. If `NULL`, groups are
  sorted automatically. Defaults to `NULL`.

- celltype_levels:

  A character vector specifying cell-type order. If `NULL`, cell types
  are sorted automatically. Defaults to `NULL`.

- n_bin:

  A numeric value specifying the number of circular histogram bins.
  Defaults to `24`.

- min_cell:

  A numeric value specifying the minimum number of cells required for a
  celltype-by-group combination. Defaults to `10`.

- normalize_distance:

  A logical value indicating whether DTW distance should be normalized
  by warping path length. Defaults to `TRUE`.

- Verbose:

  A logical value indicating whether progress messages should be
  printed. Defaults to `FALSE`.

## Value

An object of class `scRhythmDTWSimilarity`, containing:

- `data`:

  Merged cell-level phase and metadata table.

- `histogram`:

  Circular histogram table for each celltype-by-group combination.

- `combo`:

  Cell-count summary for retained combinations.

- `annotation`:

  Annotation table matched to rows and columns of the matrices.

- `distance`:

  Pairwise DTW distance matrix.

- `similarity`:

  Pairwise 0-to-1 DTW similarity matrix.

- `params`:

  Parameters used for the analysis.

## Details

Cells are matched between `res_phase` and `metadata` by the `cell`
column.

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

dtw_res
dim(dtw_res$similarity)
dtw_res$similarity[1:3, 1:3]
head(dtw_res$histogram)
} # }
```
