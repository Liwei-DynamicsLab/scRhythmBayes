# Plot module-level rhythm remodeling in polar coordinates

This function draws a polar scatter plot for one rhythm-state module.
Each gene is represented by two points, one for each compared group. The
angular coordinate represents peak phase, the radial coordinate
represents rhythmic amplitude, and line segments connect the same gene
across the two groups.

## Usage

``` r
plot_module_polar(
  gene_meta,
  module_id,
  group1,
  group2,
  q_cut = 0.05,
  amp_cut = 0.1,
  n_gene = 30,
  n_label = 8,
  pt_size = 3.2,
  label_size = 3.2,
  title = NULL,
  palette = NULL,
  Verbose = FALSE
)
```

## Arguments

- gene_meta:

  A data.frame returned by
  [`infer_rhythm_module()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_rhythm_module.md).
  It must contain `"gene"` and `"module_map"` columns, as well as
  group-specific rhythm parameter columns.

- module_id:

  An integer or character value specifying the module to plot.

- group1:

  A character vector specifying the first group. This should usually
  match `group1` used in
  [`infer_rhythm_module()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_rhythm_module.md).

- group2:

  A character vector specifying the second group. This should usually
  match `group2` used in
  [`infer_rhythm_module()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_rhythm_module.md).

- q_cut:

  A numeric value specifying the FDR cutoff used to define rhythmic
  genes. Defaults to `0.05`.

- amp_cut:

  A numeric value specifying the minimum amplitude required for selected
  genes. Defaults to `0.10`.

- n_gene:

  An integer specifying the number of top-ranked genes to plot. Defaults
  to `30`.

- n_label:

  An integer specifying the number of top-ranked genes to label.
  Defaults to `8`.

- pt_size:

  A numeric value specifying point size. Defaults to `3.2`.

- label_size:

  A numeric value specifying gene-label size. Defaults to `3.2`.

- title:

  A character value specifying the plot title. If `NULL`, the title is
  set to `"M"` followed by `module_id`. Defaults to `NULL`.

- palette:

  A named list specifying optional palettes for `module` and `group`. If
  `NULL`, default module colors and default scRhythmBayes group colors
  are used. Defaults to `NULL`.

- Verbose:

  Logical; if `TRUE`, progress messages are printed.

## Value

Invisibly returns a list containing the processed plotting data, paired
gene table, selected labels, top genes, parameter-column map, and a
ggplot object.

## Details

Genes are selected from a module result table returned by
[`infer_rhythm_module()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_rhythm_module.md).
The function automatically detects group-specific rhythm parameter
columns such as `"Phase_A"`, `"Amplitude_A"`, and `"FDR_A"`. The
alternative column style `"phi_A"`, `"A_A"`, and `"q_A"` is also
supported.

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

p_polar <- plot_module_polar(
  gene_meta = res_module$module_result,
  module_id = 7,
  group1 = "A",
  group2 = "B"
)

p_polar$plot
} # }
```
