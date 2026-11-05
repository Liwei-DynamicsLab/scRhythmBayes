# Plot model-fitted Target-versus-Rest rhythmic curves

Reproduces the model-fitted Target-versus-Rest curve design used for the
original OVR figure. Target and pooled-Rest trajectories are calculated
from the fitted mesor, amplitude, and phase returned by
[`infer_ovr_rhythm()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_ovr_rhythm.md).
Each trajectory is divided by its mean so that the horizontal reference
value is one. No additional smoother is fitted.

## Usage

``` r
plot_ovr_curve(
  ovr,
  genes = NULL,
  target_group,
  group_by = NULL,
  n_grid = 400,
  target_color = NULL,
  rest_color = "#707070",
  ncol = NULL
)
```

## Arguments

- ovr:

  A result returned by
  [`infer_ovr_rhythm()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_ovr_rhythm.md),
  or its `ovr_result` table.

- genes:

  Optional character vector of genes to plot. If `NULL`, all OVR markers
  for `target_group` are used.

- target_group:

  Single target group to plot against pooled Rest.

- group_by:

  Group column in `ovr_result`. It is read from an OVR result list or
  inferred when the table contains exactly one non-result column.

- n_grid:

  Number of phase-grid points over 24 hours.

- target_color:

  Target-curve color. If `NULL`, the original OVR cell-type palette is
  matched to the target group order.

- rest_color:

  Color for the pooled-Rest curve.

- ncol:

  Optional number of facet columns for multiple genes.

## Value

A `ggplot` object. The model-derived values used for plotting are stored
in the `curve_data` attribute.

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
  Verbose = FALSE
)

res_ovr <- infer_ovr_rhythm(
  count = srb_example_data$Count,
  metadata = srb_example_data$meta,
  theta_pred = res_phase$theta_pred,
  r_cell = res_phase$r_cell,
  group_by = "celltype",
  ZT_by = "ZT",
  target_groups = "1",
  sample_by = "sample",
  B = 3,
  n_worker = 1,
  Verbose = FALSE
)

plot_ovr_curve(
  res_ovr,
  genes = res_ovr$ovr_markers$gene[1],
  target_group = "1"
)
} # }
```
