# Infer one-vs-rest rhythmic markers

For each target group, all remaining cells are pooled as Rest. Target
and Rest are fitted with infer_gene_param() using the same inputs.

## Usage

``` r
infer_ovr_rhythm(
  count,
  metadata,
  theta_pred,
  r_cell,
  group_by,
  ZT_by,
  target_groups = NULL,
  q_cut = 0.05,
  amp_cut = 0.25,
  n_worker = NULL,
  setseed = 1,
  Verbose = TRUE,
  sample_by = NULL,
  ...
)
```

## Arguments

- count:

  Gene-by-cell raw count matrix.

- metadata:

  Cell-level metadata aligned by the `cell` column.

- theta_pred:

  Inferred cell phases in radians.

- r_cell:

  Cell-level rhythmic reliability values between 0 and 1.

- group_by:

  Metadata column containing the groups to compare.

- ZT_by:

  Metadata column containing zeitgeber-time bins.

- target_groups:

  Optional groups to test. By default, all groups are used.

- q_cut:

  FDR threshold for target-specific rhythmic markers.

- amp_cut:

  Minimum fitted target amplitude.

- n_worker:

  Number of parallel workers passed to
  [`infer_gene_param()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_gene_param.md).

- setseed:

  Random seed.

- Verbose:

  Logical; whether to print progress messages.

- sample_by:

  Optional metadata column identifying independent biological samples
  for sample-aware inference.

- ...:

  Additional arguments passed to
  [`infer_gene_param()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_gene_param.md),
  such as `B` and `sample_boot_by`.

## Value

A list containing `ovr_result`, `ovr_markers`, `ovr_summary`, the name
of `group_by`, and analysis settings. The result table retains fitted
Target and Rest Mesor, mean, amplitude, and phase values so that
model-based curves can be drawn with
[`plot_ovr_curve()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/plot_ovr_curve.md).
