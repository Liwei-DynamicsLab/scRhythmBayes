---
title: "scRhythmBayes"
output: github_document
---

# scRhythmBayes <img src="man/figures/logo.png" align="right" width="160">

[![R-CMD-check](https://github.com/Liwei-DynamicsLab/scRhythmBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Liwei-DynamicsLab/scRhythmBayes/actions/workflows/R-CMD-check.yaml)

**Cell-level circadian phase inference and rhythm-state remodeling for single-cell transcriptomes**

`scRhythmBayes` is an R package for phase-resolved circadian analysis in heterogeneous single-cell transcriptomic data.

The package is designed to infer latent cell-level circadian phase, estimate gene-level rhythmic parameters, and characterize condition- and cell-type-resolved rhythm-state remodeling. It provides a core workflow for cell phase inference, rhythmic parameter estimation, and differential rhythm-state module analysis, together with extension analyses for one-vs-rest rhythmic program discovery, DTW-based phase-distribution similarity, and circular phase-concentration rho analysis.

<p align="center">
  <img src="man/figures/readme.png" width="760">
</p>

## Overview

`scRhythmBayes` organizes single-cell circadian analysis into one core workflow and several extension workflows.

### Core workflow

The core workflow follows the main analysis logic of the method:

1. **Cell-level circadian phase inference**  
   Infer latent circadian phase for each cell from seed rhythmic genes.

2. **Gene-level rhythmic parameter inference**  
   Estimate rhythmic parameters, including mesor, mean expression, amplitude, phase, dropout rate, likelihood-ratio statistics, and FDR.

3. **Differential rhythm-state module analysis**  
   Classify genes into interpretable rhythm-state remodeling modules between two biological conditions.

### Extension analyses

Additional analysis branches are provided for heterogeneous or multi-cell-type datasets:

1. **One-vs-rest rhythmic program discovery**  
   Identify cell-type-specific rhythmic programs by comparing one target cell type against all remaining cell types.

2. **DTW-based phase-distribution similarity analysis**  
   Compare circular cell-phase distributions across celltype-by-condition combinations.

3. **Circular phase-concentration rho analysis**  
   Summarize the concentration of inferred cell phases within cell types, groups, and optional sample-level units.

## Installation

The package can be installed from GitHub:

```{r install, eval=FALSE}
install.packages("remotes")
remotes::install_github("Liwei-DynamicsLab/scRhythmBayes")
```

Load the package:

```{r load-package, eval=FALSE}
library(scRhythmBayes)
```

Check the installed version:

```{r version, eval=FALSE}
packageVersion("scRhythmBayes")
```

Note: if the GitHub repository is private, installation requires GitHub authentication and repository access.

## Example data

`scRhythmBayes` includes an example dataset for testing the complete workflow.

```{r example-data, eval=FALSE}
data("srb_example_data")

names(srb_example_data)
dim(srb_example_data$Count)
dim(srb_example_data$meta)
head(srb_example_data$meta)
```

The example object contains:

| Object | Description |
|---|---|
| `srb_example_data$Count` | Gene-by-cell sparse count matrix (1,000 genes x 4,000 cells) |
| `srb_example_data$meta` | Cell-level metadata (cell ID, sampling time, cell type, condition, simulated true phase) |

Basic input objects:

```{r basic-inputs, eval=FALSE}
count <- srb_example_data$Count
metadata <- srb_example_data$meta
```

## Core workflow

### 1. Infer cell-level circadian phase

```{r infer-phase, eval=FALSE}
res_phase <- infer_cell_phase(
  count = count,
  metadata = metadata,
  ZT_by = "ZT",
  strata_by = "celltype",
  Species = "mouse",
  setseed = 1,
  Verbose = FALSE
)
```

Output columns:

| Column | Description |
|---|---|
| `cell` | Cell barcode or ID |
| `ZT_input` | Observed sampling Zeitgeber time |
| `stratum` | Stratum used during phase inference |
| `theta_pred` | Inferred cell-level circadian phase in radians |
| `ZT_pred` | Converted 24-hour circadian time |
| `Rbar` | Aggregate phase concentration |
| `r_cell` | Cell-level phase reliability weight |

### 2. Estimate gene-level rhythmic parameters

```{r infer-param, eval=FALSE}
cell_use <- metadata$celltype == "1"

res_param <- infer_gene_param(
  count = count[, cell_use, drop = FALSE],
  metadata = metadata[cell_use, , drop = FALSE],
  theta_pred = res_phase$theta_pred[cell_use],
  r_cell = res_phase$r_cell[cell_use],
  strata_by = "condition",
  B = 10,
  n_worker = 1,
  setseed = 1,
  Verbose = FALSE
)
```

Outputs:

| Output | Description |
|---|---|
| `param_long` | Long-format gene-level parameters |
| `param_wide` | Wide-format paired table for module inference |
| `bootstrap_res` | Optional bootstrap confidence intervals |

### 3. Infer differential rhythm-state modules

```{r infer-module, eval=FALSE}
res_module <- infer_rhythm_module(
  param_wide = res_param$param_wide,
  compare_by = "condition",
  group1 = "A",
  group2 = "B",
  out_dir = NULL,
  Verbose = FALSE
)
```

### 4. Visualize rhythm-state modules

```{r plot-module-heatmap, eval=FALSE}
hm_module <- plot_module_heatmap(
  count = count[, cell_use, drop = FALSE],
  res_phase = res_phase[cell_use, , drop = FALSE],
  gene_meta = res_module$module_result,
  metadata = metadata[cell_use, , drop = FALSE],
  group_by = "condition",
  group1 = "A",
  group2 = "B"
)

hm_module
plot(hm_module, module_id = 0)
```

```{r plot-module-polar, eval=FALSE}
p_polar <- plot_module_polar(
  gene_meta = res_module$module_result,
  module_id = 7,
  group1 = "A",
  group2 = "B"
)

p_polar$plot
```

## Extension workflows

### One-vs-rest rhythmic program discovery

The one-vs-rest branch is designed for heterogeneous tissues or multi-cell-type datasets. It compares each target cell type against all remaining cell types to identify cell-type-specific rhythmic programs.

```{r ovr-workflow, eval=FALSE}
res_ovr <- infer_ovr_rhythm(
  param_long = res_param$param_long,
  group_by = "celltype",
  q_cut = 0.05,
  amp_cut = 0.25,
  Verbose = FALSE
)
```

### DTW-based phase-distribution similarity

The DTW branch compares circular phase distributions between celltype-by-condition combinations.

```{r dtw-workflow, eval=FALSE}
dtw_res <- run_dtw_similarity(
  res_phase = res_phase,
  metadata = metadata,
  celltype_by = "celltype",
  group_by = "condition",
  min_cell = 10,
  Verbose = FALSE
)
```

### Circular phase-concentration rho analysis

The rho branch summarizes how concentrated inferred cell phases are within each cell type and condition.

```{r rho-workflow, eval=FALSE}
rho_res <- run_phase_rho(
  res_phase = res_phase,
  metadata = metadata,
  celltype_by = "celltype",
  group_by = "condition",
  theta_by = "theta_pred",
  min_cell = 10,
  Verbose = FALSE
)
```

## Tutorials

Detailed tutorials are available through the pkgdown documentation website. These tutorials mirror the main analysis branches of `scRhythmBayes` and include step-by-step examples, output interpretation, and figure-generation workflows.

### Core workflow

- [Basic workflow: from cell phase inference to rhythm-state modules](https://liwei-dynamicslab.github.io/scRhythmBayes/articles/basic-workflow.html)

### Extension analyses

- [Rhythm-state module analysis](https://liwei-dynamicslab.github.io/scRhythmBayes/articles/module-analysis.html)
- [One-vs-rest rhythmic program discovery](https://liwei-dynamicslab.github.io/scRhythmBayes/articles/ovr-analysis.html)
- [DTW and rho phase-prioritization analysis](https://liwei-dynamicslab.github.io/scRhythmBayes/articles/phase-prioritization-dtw-rho.html)

The full documentation website is available at:

- [scRhythmBayes documentation](https://liwei-dynamicslab.github.io/scRhythmBayes/)

## Benchmark and analysis scripts

Analysis and benchmark scripts will be added to this repository prior to publication. These scripts will include:

- simulation scripts used to generate benchmark datasets;
- phase-inference benchmark scripts;
- gene-level rhythmic parameter benchmark scripts;
- rhythm-state module classification benchmark scripts;
- robustness analyses for seed-gene number, zeitgeber-prior sensitivity, temporal sampling density and cell number;
- runtime benchmark scripts;
- figure-generation scripts for the main and supplementary figures.

Large intermediate files and processed outputs will not be stored directly in this repository when file size is prohibitive. These files will be deposited separately or linked through a public repository upon publication.

## Installation and dependencies

All analyses were performed in R version 4.2.0 or higher. Core package dependencies include Matrix, dplyr, tibble, tidyr, future, furrr, purrr, readr, mvtnorm, ggplot2, ggrepel, scales, mgcv, ComplexHeatmap, circlize, RColorBrewer, dtw, circular and rlang.

Downstream analyses additionally used Seurat, SeuratObject and pathway-analysis packages including ReactomePA, clusterProfiler, org.Mm.eg.db and org.Hs.eg.db. Full package versions and reproducibility information will be provided in the repository.

## Main functions

### Core workflow

| Function | Description |
|---|---|
| `infer_cell_phase()` | Infer cell-level circadian phase |
| `infer_gene_param()` | Estimate gene-level rhythmic parameters |
| `infer_rhythm_module()` | Infer differential rhythm-state modules |
| `plot_module_heatmap()` | Plot rhythm-state module heatmap |
| `plot_module_polar()` | Plot rhythm modules on polar coordinates |

### Phase visualization

| Function | Description |
|---|---|
| `plot_phase_scatter()` | Compare inferred phase with sampling time |
| `plot_phase_circle()` | Visualize circular cell-phase distributions |
| `plot_phase_flower()` | Visualize phase distributions by sampling time |
| `plot_phase_rose()` | Plot rose-style phase distributions |

### Extension analyses

| Function | Description |
|---|---|
| `infer_ovr_rhythm()` | Identify one-vs-rest cell-type-specific rhythmic programs |
| `plot_ovr_heatmap()` | Plot one-vs-rest rhythm-state heatmap |
| `run_dtw_similarity()` | Run DTW-based phase-distribution similarity analysis |
| `plot_dtw_heatmap()` | Plot DTW similarity or distance heatmap |
| `run_phase_rho()` | Summarize circular phase concentration rho |
| `plot_phase_rho()` | Plot rho by cell type and group |

## Notes

`scRhythmBayes` is intended for phase-resolved analysis of single-cell transcriptomic datasets where cells are collected across circadian or time-of-day sampling points. The core workflow focuses on latent cell-phase inference and rhythm-state remodeling analysis. The extension workflows support additional biological questions, including cell-type-specific rhythmic program discovery and comparison of inferred phase distributions across cell types and conditions.

The package is under active development. Function interfaces and tutorials may be updated as the manuscript and reproducibility workflows are finalized.

## Citation

If you use `scRhythmBayes`, please cite the associated manuscript:

Zhang Q, Li C, Zhang L. Decoding cell-level circadian phase and rhythm-state remodeling in single-cell transcriptomes with scRhythmBayes.

## License

This package is released under the MIT License.

## Authors

Liwei Zhang, Qian Zhang, and Chenyu Li.


