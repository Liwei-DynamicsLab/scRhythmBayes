# Plot cell phases on a cosine-sine circle

This function plots inferred cell phases on the unit circle using cosine
and sine coordinates. Cells are matched between `res_phase` and
`metadata` by the `cell` column.

## Usage

``` r
plot_phase_circle(
  res_phase,
  metadata,
  zt_by = "ZT",
  group_by = "celltype",
  theta_by = "theta_pred",
  group_levels = NULL,
  zt_levels = NULL,
  zt_colors = NULL,
  point_size = 3.5,
  point_alpha = 0.4,
  title = NULL
)
```

## Arguments

- res_phase:

  A data.frame returned by
  [`infer_cell_phase()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_cell_phase.md).

- metadata:

  A data.frame containing cell-level metadata.

- zt_by:

  A character vector specifying the column name in `metadata` that
  contains sampling Zeitgeber time in hours. Defaults to `"ZT"`.

- group_by:

  A character vector specifying the column name in `metadata` used for
  faceting. If `NULL`, all cells are plotted in one panel. Defaults to
  `"celltype"`.

- theta_by:

  A character vector specifying the column name in `res_phase` that
  contains inferred cell phase in radians. Defaults to `"theta_pred"`.

- group_levels:

  A character vector specifying the order of facet groups. If `NULL`,
  groups are sorted automatically. Defaults to `NULL`.

- zt_levels:

  A numeric vector specifying the order of ZT groups. If `NULL`, ZT
  levels are sorted automatically. Defaults to `NULL`.

- zt_colors:

  A named character vector specifying colors for ZT groups. Names can be
  formatted ZT labels, such as `"ZT00"` and `"ZT04"`, or numeric hour
  labels, such as `"0"` and `"4"`. If `NULL`, the default scRhythmBayes
  phase palette is used. Defaults to `NULL`.

- point_size:

  A numeric value specifying the point size. Defaults to `3.5`.

- point_alpha:

  A numeric value specifying the point transparency. Defaults to `0.40`.

- title:

  A character vector specifying the plot title. If `NULL`, no title is
  shown. Defaults to `NULL`.

## Value

Invisibly returns a list containing:

- `plot`:

  A ggplot object.

- `data`:

  The processed plotting data.

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

p_circle <- plot_phase_circle(
  res_phase = res_phase,
  metadata = srb_example_data$meta,
  zt_by = "ZT",
  group_by = "celltype"
)

p_circle$plot
} # }
```
