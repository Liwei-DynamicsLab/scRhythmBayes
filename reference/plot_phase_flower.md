# Plot circular cell-phase density by sampling time

This function plots circular density summaries of inferred cell phases
across sampling Zeitgeber time groups. Cells are matched between
`res_phase` and `metadata` by the `cell` column.

## Usage

``` r
plot_phase_flower(
  res_phase,
  metadata,
  zt_by = "ZT",
  group_by = "celltype",
  theta_by = "theta_pred",
  group_levels = NULL,
  zt_levels = NULL,
  zt_colors = NULL,
  title = NULL,
  Verbose = FALSE
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

- title:

  A character vector specifying the plot title. If `NULL`, a default
  title is used. Defaults to `NULL`.

- Verbose:

  A logical value indicating whether progress messages should be
  printed. Defaults to `FALSE`.

## Value

Invisibly returns a list containing:

- `plot`:

  A ggplot object.

- `data`:

  The processed cell-level plotting data.

- `density`:

  The circular density data.

- `ring`:

  The ZT ring annotation data.

- `cell_count`:

  Cell-count summary by group.

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

p_flower <- plot_phase_flower(
  res_phase = res_phase,
  metadata = srb_example_data$meta,
  zt_by = "ZT",
  group_by = "celltype"
)

p_flower$plot
} # }
```
