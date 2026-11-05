# scRhythmBayes example-data validation

Date: 2026-11-05

## Data design

- Count: 138 genes by 2,160 cells; sparse `.rda` size 342 KB.
- Samples: 36 independent biological samples.
- Design: 2 conditions by 6 ZT values by 3 replicates.
- Cells: 60 per sample; 20 from each of three cell types.
- Truth: M0-M7 modules, 24 mouse phase-seed genes, 18 OVR marker genes,
  cell phases, and native sample-specific mesor and phase effects.
- Count model: negative binomial with mean-dependent zero inflation.

## End-to-end results

The candidate data were tested in a clean R session before replacing
`data/srb_example_data.rda`.

| Check | Result |
|---|---:|
| `infer_cell_phase()` runtime | 4.16 s |
| Matched/retained default mouse seeds | 23 / 23 |
| Median absolute phase error | 0.37 h |
| 90th percentile absolute phase error | 0.96 h |
| Sample-aware gene test runtime | 34.89 s |
| Independent samples per condition | 18 / 18 |
| Finite sample-aware P values and FDR | 276 / 276 |
| Sample test status | 276 `ok_cluster_multiplier_score` |
| M0-M7 seed/module accuracy | 97.9% / 97.9% |
| Full three-target OVR runtime | 201.11 s |
| OVR markers by cell type | 5 / 6 / 6 |
| True target-matched markers among calls | 17 / 17 |

`plot_ovr_curve()` was tested with four genes and produced finite Target and
Rest trajectories directly from fitted Mesor, amplitude, and phase values.

## Automated checks

- `devtools::load_all()` completes without missing OVR export or S3 warnings.
- `devtools::test()` passes 141 tests with 0 failures, 0 warnings, and 0 skips.
- `data-raw/validate_srb_example_data.R` contains the reproducible full
  end-to-end validation used for future data regeneration.

The complete OVR workflow is intentionally kept out of routine unit tests
because it takes about three minutes on the tested Mac. Structural and curve
formula checks remain in `tests/testthat/`.

## Tutorial validation

- All executable README analysis chunks were run sequentially in a clean
  loaded package session, including the model-fitted OVR curve.
- The four focused vignettes were rendered with the sample-aware parameters,
  real biological sample IDs, and generated tutorial figures.
- `R CMD check` completed with 0 errors, 0 warnings, and 0 notes while
  rebuilding all vignette outputs with `SRB_RUN_OVR_TUTORIAL=true`.
