# Infer rhythm-state modules

This function assigns genes to rhythm-state modules based on two-group
gene-level rhythmic parameter estimates.

## Usage

``` r
infer_rhythm_module(
  param_wide,
  compare_by = "condition",
  group1,
  group2,
  tau_q = 0.05,
  tau_A = 0.25,
  tau_dA = 1,
  tau_dphi = pi/6,
  out_dir = NULL,
  Verbose = TRUE
)
```

## Arguments

- param_wide:

  A wide-format data.frame returned by
  [`infer_gene_param()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_gene_param.md),
  with at least two unique genes.

- compare_by:

  Comparison variable name. Defaults to `"condition"`.

- group1:

  Reference group (baseline). This is not a plotting position.

- group2:

  Target group. Gain/loss, amplitude enhancement/attenuation and phase
  advance/delay describe group2 relative to group1. For CV-to-PV liver
  comparisons, use group1 = "Hep.CV", group2 = "Hep.PV".

- tau_q:

  FDR threshold for significant rhythmicity. Defaults to `0.05`.

- tau_A:

  Minimum rhythmic amplitude threshold. Defaults to `0.25`.

- tau_dA:

  Log2 amplitude-ratio threshold. Defaults to `1.00`.

- tau_dphi:

  Phase-shift threshold in radians. Defaults to `pi / 6`.

- out_dir:

  Output directory for result tables. If `NULL`, no files are written.
  Defaults to `NULL`.

- Verbose:

  Logical; if `TRUE`, progress messages are printed.

## Value

A list containing rhythm-state module assignment results and module
summaries. The main components are `module_result` and `module_summary`.
`comparison` records reference/target identities and effect conventions.
`input_diagnostics` flags incomplete original two-group parameter sets
before any non-rhythmic-side display masking. Incomplete inputs trigger
a warning but retain the existing classification algorithm; a posterior
assignment does not establish unestimated rhythm states.

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

cell_use <- srb_example_data$meta$celltype == "1"

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
  Verbose = FALSE
)

head(res_module$module_result)
} # }
```
