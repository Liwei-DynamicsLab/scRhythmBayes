# Example single-cell circadian dataset

A simulated example dataset used to demonstrate the scRhythmBayes
workflow. The object contains count data and cell-level metadata for
testing cell-phase inference, gene-level rhythmic parameter estimation,
rhythm-state module analysis, OVR analysis, DTW similarity analysis, and
phase-rho summaries.

## Usage

``` r
srb_example_data
```

## Format

A list with three elements:

- Count:

  A sparse 138 gene by 2,160 cell count matrix.

- meta:

  Cell metadata containing `cell`, biological `sample`, `replicate`,
  `ZT`, `celltype`, `condition`, true phase (`T_true`, `Theta_true`),
  and simulated noisy phase inputs (`Theta_noisy`, `r_sim`). The data
  contain 36 samples: two conditions by six ZT values by three
  replicates.

- truth:

  A list containing gene-level M0–M7 and OVR truth, sample design,
  native sample-specific mesor and phase effects, and the simulation
  settings.

## Source

Simulated for the scRhythmBayes teaching workflow with
`data-raw/generate_srb_example_data.R`.

## Details

Counts were generated from a sample-aware negative-binomial model with
mean-dependent zero inflation. Biological samples have native
sample-specific mesor and phase perturbations. The 24 phase-seed genes
use mouse gene symbols, including `Arntl`, `Clock`, `Cry1`, `Cry2`, and
`Per1--Per3`. Simulated phase columns are retained as ground truth and
are not outputs from
[`infer_cell_phase()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_cell_phase.md).
