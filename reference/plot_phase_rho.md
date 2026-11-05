# Plot cell-phase rho by cell type

This function plots circular concentration rho values from a
`scRhythmPhaseRho` object returned by
[`run_phase_rho()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/run_phase_rho.md).

## Usage

``` r
plot_phase_rho(
  rho_res,
  group_colors = NULL,
  celltype_order = NULL,
  show_sig = TRUE,
  facet_scales = "free_y",
  point_size = 2.2,
  point_alpha = 0.85,
  title = NULL
)
```

## Arguments

- rho_res:

  A `scRhythmPhaseRho` object returned by
  [`run_phase_rho()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/run_phase_rho.md).

- group_colors:

  A named character vector specifying colors for groups. If `NULL`,
  default scRhythmBayes group colors are used. Defaults to `NULL`.

- celltype_order:

  A character vector specifying the order of cell types. If `NULL`, cell
  types are ordered by mean rho. Defaults to `NULL`.

- show_sig:

  A logical value indicating whether significance labels should be
  displayed when sample-level significance results are available.
  Defaults to `TRUE`.

- facet_scales:

  A character vector specifying facet scales for sample-level faceted
  boxplots. Defaults to `"free_y"`.

- point_size:

  A numeric value specifying point size. Defaults to `2.2`.

- point_alpha:

  A numeric value specifying point transparency. Defaults to `0.85`.

- title:

  A character vector specifying the plot title. If `NULL`, no title is
  shown. Defaults to `NULL`.

## Value

Invisibly returns a list containing:

- `plot`:

  The ggplot object.

- `data`:

  The plotting data.

- `sig`:

  The significance-label data, if available.

- `mode`:

  The plotting mode, either `"sample_boxplot"` or `"group_dotplot"`.

- `params`:

  Plotting parameters used by the function.

## Details

If sample-level rho results are available, the function draws boxplots
across samples within each cell type. If no sample-level results are
available, the function draws a group-level dot plot with one rho value
per cell type and group. Group-level rho plots are descriptive only and
do not support statistical testing.

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

rho_res <- run_phase_rho(
  res_phase = res_phase,
  metadata = srb_example_data$meta,
  celltype_by = "celltype",
  group_by = "condition",
  theta_by = "theta_pred",
  group_levels = c("A", "B"),
  min_cell = 10,
  Verbose = FALSE
)

p_rho <- plot_phase_rho(rho_res)
p_rho$plot
} # }
```
