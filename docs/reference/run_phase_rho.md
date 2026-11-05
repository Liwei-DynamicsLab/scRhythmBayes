# Run cell-phase rho summary

This function summarizes the circular concentration of inferred cell
phases within each cell type and group. If a sample-level column is
provided, rho is also summarized at the cell type, group, and sample
level, and pairwise group comparisons are performed within each cell
type.

## Usage

``` r
run_phase_rho(
  res_phase,
  metadata,
  celltype_by = "celltype",
  group_by = "condition",
  sample_by = NULL,
  theta_by = "theta_pred",
  celltype_levels = NULL,
  group_levels = NULL,
  compare_groups = NULL,
  min_cell = 10,
  p_adjust_method = "BH",
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

- celltype_by:

  A character vector specifying the column name in `metadata` used as
  the cell-type annotation. Defaults to `"celltype"`.

- group_by:

  A character vector specifying the column name in `metadata` used as
  the group or condition annotation. Defaults to `"condition"`.

- sample_by:

  A character vector specifying the column name in `metadata` used as
  the sample annotation. If `NULL`, only group-level rho is computed.
  Defaults to `NULL`.

- theta_by:

  A character vector specifying the column name in `res_phase`
  containing inferred cell phase in radians. Defaults to `"theta_pred"`.

- celltype_levels:

  A character vector specifying the order of cell types. If `NULL`,
  levels are sorted automatically. Defaults to `NULL`.

- group_levels:

  A character vector specifying the order of groups. If `NULL`, levels
  are sorted automatically. Defaults to `NULL`.

- compare_groups:

  A character vector of length two specifying two groups to compare in
  sample-level rho tests. If `NULL`, the first two group levels are
  used. Defaults to `NULL`.

- min_cell:

  A numeric value specifying the minimum number of cells required for
  each rho estimate. Defaults to `10`.

- p_adjust_method:

  A character vector specifying the multiple-testing correction method
  passed to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).
  Defaults to `"BH"`.

- Verbose:

  A logical value indicating whether progress messages should be
  printed. Defaults to `FALSE`.

## Value

An object of class `scRhythmPhaseRho`, containing:

- `data`:

  Cleaned cell-level phase and metadata table.

- `rho_result`:

  Group-level rho estimates for each cell type and group.

- `rho_sample`:

  Sample-level rho estimates, returned only when `sample_by` is
  provided.

- `rho_sample_summary`:

  Sample-level rho summary by cell type and group, returned only when
  `sample_by` is provided.

- `sig_result`:

  Sample-level group-comparison results, returned only when `sample_by`
  is provided and comparison groups are available.

- `params`:

  Parameters used for rho summarization.

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

rho_res <- run_phase_rho(
  res_phase = res_phase,
  metadata = srb_example_data$meta,
  celltype_by = "celltype",
  group_by = "condition",
  group_levels = c("A", "B"),
  min_cell = 10,
  Verbose = TRUE
)

rho_res
head(rho_res$rho_result)
} # }
```
