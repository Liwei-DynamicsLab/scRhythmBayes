# Plot one-vs-rest rhythmic marker heatmap

This function builds a full-matrix heatmap for one-vs-rest target
rhythmic markers across multiple cell types. Genes are grouped by their
target cell type and ordered by target-group peak phase. Columns are
grouped by cell type and smoothed along inferred cell phase.

## Usage

``` r
plot_ovr_heatmap(
  count,
  res_phase,
  ovr_result,
  metadata,
  celltype_by = "celltype",
  theta_by = "theta_pred",
  scale_factor = 10000,
  palette = NULL,
  Verbose = FALSE
)
```

## Arguments

- count:

  A gene-by-cell count matrix. Row names must be gene names and column
  names must be cell names.

- res_phase:

  A data.frame returned by
  [`infer_cell_phase()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_cell_phase.md).
  It must contain `"cell"` and the column specified by `theta_by`.

- ovr_result:

  A data.frame returned by
  [`infer_ovr_rhythm()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_ovr_rhythm.md),
  usually `res_ovr$ovr_result`. Current results use the column named by
  `celltype_by`, `"Phase_target"`, and `"ovr_marker"`; legacy column
  names are also accepted.

- metadata:

  A data.frame containing cell-level metadata. Cells are matched to
  `count` by the `cell` column if available, otherwise by row names.

- celltype_by:

  A character vector specifying the column name in `metadata` used for
  cell-type groups. Defaults to `"celltype"`.

- theta_by:

  A character vector specifying the column name in `res_phase`
  containing inferred cell phase in radians. Defaults to `"theta_pred"`.

- scale_factor:

  A numeric value passed to
  [`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html).
  Defaults to `10000`.

- palette:

  A named list specifying optional palettes for `celltype`, `phase`, and
  `heatmap`. If `NULL`, default scRhythmBayes palettes are used.
  Defaults to `NULL`.

- Verbose:

  Logical; if `TRUE`, progress messages are printed. Defaults to
  `FALSE`.

## Value

An object of class `scRhythmOVRHeatmap`.

## Details

This function accepts a gene-by-cell count matrix and uses Seurat
internally for log-normalization and scaling. Users do not need to
provide a Seurat object.

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

res_ovr <- infer_ovr_rhythm(
  count = srb_example_data$Count,
  metadata = srb_example_data$meta,
  theta_pred = res_phase$theta_pred,
  r_cell = res_phase$r_cell,
  group_by = "celltype",
  ZT_by = "ZT",
  B = 5,
  n_worker = 1,
  setseed = 1,
  Verbose = FALSE
)

ht_ovr <- plot_ovr_heatmap(
  count = srb_example_data$Count,
  res_phase = res_phase,
  ovr_result = res_ovr$ovr_result,
  metadata = srb_example_data$meta,
  celltype_by = "celltype"
)

ht_ovr
plot(ht_ovr)
} # }
```
