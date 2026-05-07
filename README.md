# scRhythmBayes <img src="man/figures/logo.png" align="right" width="160">

[![R-CMD-check](https://github.com/Liwei-DynamicsLab/scRhythmBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Liwei-DynamicsLab/scRhythmBayes/actions/workflows/R-CMD-check.yaml)

**Phase-resolved circadian analysis for single-cell transcriptomes**

`scRhythmBayes` is an R package for decoding cell-level circadian phase and characterizing rhythm-state remodeling in heterogeneous single-cell transcriptomic data.

The package provides tools for cell-level circadian phase inference, gene-level rhythmic modelling, rhythm module analysis, one-vs-rest cell-type-specific rhythm-state analysis, DTW-based phase-distribution similarity analysis, phase-concentration rho analysis, and visualization.

<p align="center">
  <img src="man/figures/readme.png" width="760">
</p>
## Installation

The package can be installed from GitHub:

```r
install.packages("remotes")
remotes::install_github("Liwei-Zhang/scRhythmBayes")
```

Load the package:

```r
library(scRhythmBayes)
```

Check the installed version:

```r
packageVersion("scRhythmBayes")
```

Note: if the GitHub repository is private, installation requires GitHub authentication and repository access.

## Quick start

Load the example dataset:

```r
data("srb_example_data")
names(srb_example_data)
```

The example dataset contains a gene-by-cell expression matrix and cell-level metadata:

```r
count <- srb_example_data$Count
metadata <- srb_example_data$meta
```

## Infer cell-level circadian phase

```r
res_phase <- infer_cell_phase(
  count = srb_example_data$Count,
  metadata = srb_example_data$meta,
  ZT_by = "ZT",
  strata_by = "celltype",
  Species = "mouse",
  setseed = 1,
  Verbose = FALSE
)

head(res_phase)
```

## Estimate gene-level rhythmic parameters

```r
res_param <- infer_gene_param(
  count = srb_example_data$Count,
  res_phase = res_phase,
  metadata = srb_example_data$meta,
  celltype_by = "celltype",
  group_by = "condition",
  theta_by = "theta_pred",
  Verbose = FALSE
)

head(res_param)
```

## Infer rhythm modules

```r
res_module <- infer_rhythm_module(
  res_param = res_param,
  Verbose = FALSE
)

head(res_module$module_result)
```

## One-vs-rest rhythm-state analysis

```r
ovr_result <- infer_ovr_rhythm(
  res_param = res_param,
  target_by = "celltype",
  Verbose = FALSE
)

head(ovr_result$result)
```

## DTW-based phase-distribution similarity

```r
dtw_res <- run_dtw_similarity(
  res_phase = res_phase,
  metadata = srb_example_data$meta,
  celltype_by = "celltype",
  group_by = "condition",
  theta_by = "theta_pred",
  min_cell = 10,
  Verbose = FALSE
)

dtw_res$similarity
```

```r
ht_dtw <- plot_dtw_heatmap(dtw_res)
plot(ht_dtw)
```

## Phase-concentration rho analysis

```r
rho_res <- run_phase_rho(
  res_phase = res_phase,
  metadata = srb_example_data$meta,
  celltype_by = "celltype",
  group_by = "condition",
  theta_by = "theta_pred",
  min_cell = 10,
  Verbose = FALSE
)

rho_res
```

```r
p_rho <- plot_phase_rho(rho_res)
p_rho$plot
```

## Main functions

| Function | Description |
|---|---|
| `infer_cell_phase()` | Infer cell-level circadian phase |
| `infer_gene_param()` | Estimate gene-level rhythmic parameters |
| `infer_rhythm_module()` | Infer rhythm modules |
| `infer_ovr_rhythm()` | Identify one-vs-rest cell-type-specific rhythmic programs |
| `plot_phase_scatter()` | Compare inferred phase with sampling time |
| `plot_phase_circle()` | Visualize circular cell-phase distributions |
| `plot_phase_flower()` | Visualize phase distributions by sampling time |
| `plot_phase_rose()` | Plot rose-style phase distributions |
| `plot_module_heatmap()` | Plot rhythm module heatmap |
| `plot_module_polar()` | Plot rhythm modules on polar coordinates |
| `plot_ovr_heatmap()` | Plot one-vs-rest rhythm-state heatmap |
| `run_dtw_similarity()` | Run DTW-based phase-distribution similarity analysis |
| `plot_dtw_heatmap()` | Plot DTW similarity or distance heatmap |
| `run_phase_rho()` | Summarize circular phase concentration rho |
| `plot_phase_rho()` | Plot rho by cell type and group |

## Citation

If you use `scRhythmBayes`, please cite the associated manuscript:

Zhang Q, et al. Decoding cell-level circadian phase and rhythm-state remodeling in single-cell transcriptomes with scRhythmBayes.

## License

This package is released under the MIT License.

## Authors

Qian Zhang and Liwei Zhang.
