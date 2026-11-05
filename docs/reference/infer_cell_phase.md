# Infer cell-level circadian phase

This function infers cell-level circadian phase from a gene-by-cell
count matrix and a cell-level metadata table.

## Usage

``` r
infer_cell_phase(
  count,
  metadata,
  ZT_by,
  strata_by,
  seedgene = NULL,
  Species = c("auto", "human", "mouse"),
  n_grid = 60,
  n_iter_outer = 5,
  n_iter_seed = 3,
  Kappa = 5,
  p_mix = 0.85,
  Min_n = 200,
  setseed = 1,
  Verbose = TRUE,
  min_seed = 0.01,
  family = c("NB", "Poisson")
)
```

## Arguments

- count:

  A gene-by-cell raw count matrix. `count` must be either a base R
  matrix or a sparse `Matrix` object, with genes in rows and cells in
  columns.

- metadata:

  A data.frame containing cell-level metadata.

- ZT_by:

  A character string specifying the column in `metadata` that contains
  sample-level Zeitgeber time in hours.

- strata_by:

  A character string specifying the column in `metadata` used as the
  grouping variable for cell-level phase inference.

- seedgene:

  An optional character vector of seed genes used for phase inference.
  If `NULL`, built-in default seed genes are used.

- Species:

  Species used to select the built-in default seed-gene set when
  `seedgene = NULL`. One of `"auto"`, `"human"`, or `"mouse"`.

- n_grid:

  Number of discrete phase grid points used to approximate the circular
  latent phase space.

- n_iter_outer:

  Number of outer iterations for phase inference.

- n_iter_seed:

  Number of seed-gene fitting iterations within each outer iteration.

- Kappa:

  Initial concentration parameter for the von Mises prior.

- p_mix:

  Initial mixture weight between the structured von Mises component and
  the uniform background component.

- Min_n:

  Minimum number of cells required in a stratum for stratum-level
  parameter updating.

- setseed:

  Integer random seed used for reproducibility.

- Verbose:

  Logical; if `TRUE`, progress messages are printed.

- min_seed:

  Minimum fraction of supplied cells with nonzero counts for a matched
  seed gene. Defaults to 0.01.

- family:

  Count distribution used for gene fitting and phase likelihood. One of
  `"NB"` (default) or `"Poisson"`.

## Value

A tibble with one row per cell. The returned columns include `cell`,
`ZT_input`, `stratum`, `theta_pred`, `ZT_pred`, `Rbar`, and `r_cell`.
Attributes seed_genes, seed_diagnostics (all matched genes before
filtering), n_seed_candidate, n_seed_matched, n_seed_retained, min_seed
and family record the exact seed selection and count distribution.

## Details

The inferred cell phase may depend on the input analysis scope. Running
`infer_cell_phase()` on a subset of cells is not guaranteed to produce
identical results to running it on the full dataset and then subsetting,
because seed-gene fitting and stratum-level prior updating are performed
on the supplied data matrix.

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
  Verbose = TRUE
)

head(res_phase)
} # }
```
