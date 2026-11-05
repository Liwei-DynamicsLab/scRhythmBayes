# Plot paired module heatmaps across cell groups

This function builds paired module-level heatmaps showing smoothed
rhythmic expression patterns along inferred cell phase. For each
selected module, genes are ordered by module-specific reference phase
and expression is smoothed separately in two specified cell groups.

## Usage

``` r
plot_module_heatmap(
  count,
  res_phase,
  gene_meta,
  metadata,
  group_by,
  group1,
  group2,
  theta_by = "theta_pred",
  module_ids = 0:7,
  phase_levels = c(2, 6, 10, 14, 18, 22),
  scale_factor = 10000,
  palette = NULL,
  Verbose = FALSE,
  compare_groups = NULL
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

- gene_meta:

  A data.frame returned by
  [`infer_rhythm_module()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_rhythm_module.md).
  It must contain `"gene"` and `"module_map"` columns, as well as
  group-specific rhythm parameter columns such as `"Phase_A"`,
  `"Amplitude_A"`, and `"FDR_A"`. The alternative column style
  `"phi_A"`, `"A_A"`, and `"q_A"` is also supported.

- metadata:

  A data.frame containing cell-level metadata. Cells are matched to
  `count` by the `cell` column if available, otherwise by rownames.

- group_by:

  A character vector specifying the column name in `metadata` used to
  define the two paired cell groups.

- group1:

  A character vector specifying the first (left) group to plot. This
  controls display order, not the inference comparison direction.

- group2:

  A character vector specifying the second (right) group to plot. This
  controls display order, not the inference comparison direction.

- theta_by:

  A character vector specifying the column name in `res_phase`
  containing inferred cell phase in radians. Defaults to `"theta_pred"`.

- module_ids:

  An integer vector specifying modules to plot. Defaults to `0:7`.

- phase_levels:

  A numeric vector specifying displayed phase groups in hours. Defaults
  to `c(2, 6, 10, 14, 18, 22)`.

- scale_factor:

  A numeric value passed to
  [`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html).
  Defaults to `10000`.

- palette:

  A named list specifying optional palettes for `phase`, `group`, and
  `heatmap`. If `NULL`, black/white phase labels, default scRhythmBayes
  group colors, and a blue-red heatmap palette are used. Defaults to
  `NULL`.

- Verbose:

  Logical; if `TRUE`, progress messages are printed.

- compare_groups:

  Optional character vector `c(group1, group2)` specifying the original
  [`infer_rhythm_module()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_rhythm_module.md)
  comparison direction. Inferred from `gene_meta$group1/group2` when
  available; required for M2-M7 in legacy results lacking these columns.
  M2/M4 use comparison group2 and the other modules use comparison
  group1 for row ordering, regardless of display order. This selects the
  gained/stronger side for gain or amplitude changes. Default group
  colors follow comparison identities, not left/right position. Phase
  annotations label four-hour intervals by their left endpoints. Missing
  reference phases are retained last with a warning. Nonzero modules
  with no finite reference phase produce an error rather than an
  unsorted plot.

## Value

An object of class `scRhythmModuleHeatmap`. Each successful module
includes a `phase_order` table; `summary$n_phase_ordered` and
`summary$n_phase_missing` count reference phases before smoothing.

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

cell_use <- as.character(srb_example_data$meta$celltype) == "1"

res_param <- infer_gene_param(
  count = srb_example_data$Count[, cell_use, drop = FALSE],
  metadata = srb_example_data$meta[cell_use, , drop = FALSE],
  theta_pred = res_phase$theta_pred[cell_use],
  r_cell = res_phase$r_cell[cell_use],
  strata_by = "condition",
  ZT_by = "ZT",
  B = 5,
  n_worker = 1,
  setseed = 1,
  Verbose = FALSE
)

res_module <- infer_rhythm_module(
  param_wide = res_param$param_wide,
  compare_by = "condition",
  group1 = "A",
  group2 = "B",
  out_dir = NULL,
  Verbose = FALSE
)

hm <- plot_module_heatmap(
  count = srb_example_data$Count[, cell_use, drop = FALSE],
  res_phase = res_phase[cell_use, , drop = FALSE],
  gene_meta = res_module$module_result,
  metadata = srb_example_data$meta[cell_use, , drop = FALSE],
  group_by = "condition",
  group1 = "A",
  group2 = "B"
)

hm
plot(hm, module_id = 0)
} # }
```
