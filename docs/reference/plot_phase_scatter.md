# Plot predicted cell phase against true time

This function plots predicted cell phase against true cell time for
simulation or example data with known ground-truth cell time. Cells are
matched between `res_phase` and `metadata` by the `cell` column.

## Usage

``` r
plot_phase_scatter(
  res_phase,
  metadata,
  truth_by = "T_true",
  pred_by = "ZT_pred",
  point_color = "#A8C7E6",
  point_size = 0.45,
  point_alpha = 0.6,
  title = "Cell phase inference"
)
```

## Arguments

- res_phase:

  A data.frame returned by
  [`infer_cell_phase()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_cell_phase.md).

- metadata:

  A data.frame containing cell-level metadata.

- truth_by:

  A character vector specifying the column name in `metadata` that
  contains true cell time in hours. Defaults to `"T_true"`.

- pred_by:

  A character vector specifying the column name in `res_phase` that
  contains predicted cell time in hours. Defaults to `"ZT_pred"`.

- point_color:

  A character vector specifying the point color. Defaults to
  `"#A8C7E6"`.

- point_size:

  A numeric value specifying the point size. Defaults to `0.45`.

- point_alpha:

  A numeric value specifying the point transparency. Defaults to `0.60`.

- title:

  A character vector specifying the plot title. Defaults to
  `"Cell phase inference"`.

## Value

Invisibly returns a list containing:

- `plot`:

  A ggplot object.

- `data`:

  The processed plotting data.

- `metrics`:

  A data.frame containing cell number, mean absolute circular error, and
  median absolute circular error.

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

p_scatter <- plot_phase_scatter(
  res_phase = res_phase,
  metadata = srb_example_data$meta,
  truth_by = "T_true"
)

p_scatter$plot
p_scatter$metrics
} # }
```
