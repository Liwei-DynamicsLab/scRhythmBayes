# Plot cell phase rose histogram

This function plots polar histograms of inferred cell phase. Cells are
matched between `res_phase` and `metadata` by the `cell` column.

## Usage

``` r
plot_phase_rose(
  res_phase,
  metadata,
  group_by = NULL,
  split_by = NULL,
  theta_by = "theta_pred",
  n_bin = 24,
  group_colors = NULL,
  show_count = TRUE,
  same_radius = TRUE,
  title = NULL
)
```

## Arguments

- res_phase:

  A data.frame returned by
  [`infer_cell_phase()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_cell_phase.md).

- metadata:

  A data.frame containing cell-level metadata.

- group_by:

  A character vector specifying the column name in `metadata` used for
  grouping. If `NULL`, all cells are plotted together. Defaults to
  `NULL`.

- split_by:

  A character vector specifying the column name in `metadata` used for
  faceting. If `NULL`, no facet is used. Defaults to `NULL`.

- theta_by:

  A character vector specifying the column name in `res_phase` that
  contains inferred cell phase in radians. Defaults to `"theta_pred"`.

- n_bin:

  A numeric value specifying the number of circular bins. Defaults to
  `24`.

- group_colors:

  A named character vector specifying colors for groups. If `NULL`, the
  default scRhythmBayes group palette is used. Defaults to `NULL`.

- show_count:

  A logical value indicating whether bin counts should be displayed.
  Defaults to `TRUE`.

- same_radius:

  A logical value indicating whether all panels should use the same
  radial scale. Defaults to `TRUE`.

- title:

  A character vector specifying the plot title. If `NULL`, no title is
  shown. Defaults to `NULL`.

## Value

Invisibly returns a list containing:

- `plot`:

  A ggplot object.

- `data`:

  The processed cell-level plotting data.

- `histogram`:

  The phase-bin histogram data.

- `cell_count`:

  Cell-count summary by group and facet.

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

p_rose <- plot_phase_rose(
  res_phase = res_phase,
  metadata = srb_example_data$meta,
  group_by = "condition",
  split_by = "celltype"
)

p_rose$plot
} # }
```
