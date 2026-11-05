# Infer gene-level rhythmic parameters

Fits an unrestricted zero-inflated negative-binomial harmonic model
within each metadata stratum using fixed cell phase and rhythmic
probability. Time-resolved expression eligibility precedes testing, and
empirical-Bayes shrinkage provides stabilized harmonic effect-size
estimates.

## Usage

``` r
infer_gene_param(
  count,
  metadata,
  theta_pred,
  r_cell,
  strata_by,
  ZT_by,
  B = 100,
  n_worker = NULL,
  setseed = 1,
  Verbose = TRUE,
  sample_by = NULL,
  sample_boot_by = NULL
)
```

## Arguments

- count:

  Gene-by-cell raw count matrix with gene and cell names.

- metadata:

  Data frame containing cell-level metadata.

- theta_pred:

  Numeric vector of inferred cell phases in radians.

- r_cell:

  Numeric vector of cell-level posterior rhythmic probabilities.

- strata_by:

  Metadata column defining gene-level fitting strata.

- ZT_by:

  Metadata column defining zeitgeber-time or experimental-time bins.

- B:

  Number of bootstrap iterations for confidence intervals. Cells are
  resampled in cell-level mode; whole biological samples are resampled
  when `sample_by` is supplied.

- n_worker:

  Number of parallel workers; `NULL` uses available cores.

- setseed:

  Random seed.

- Verbose:

  Logical; whether to print progress messages.

- sample_by:

  Optional metadata column identifying independent biological samples.

- sample_boot_by:

  Optional metadata column stratifying sample-aware bootstrap
  resampling.

## Value

A list containing main result tables `param_long` and `param_wide`,
`bootstrap_res`, eligibility records for all input genes, fit
diagnostics, and stratum-level shrinkage metadata.

## Details

Within each fitting stratum, genes are eligible when their maximum
time-bin pseudobulk relative UMI exceeds `1e-5`; only eligible genes are
tested and bootstrapped. P values and BH FDR are calculated from the
unrestricted fit. Empirical-Bayes shrinkage stabilizes amplitude but
does not alter significance testing.

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

head(res_param$param_long)
head(res_param$param_wide)
} # }
```
