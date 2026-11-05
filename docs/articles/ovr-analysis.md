# One-vs-rest rhythmic program discovery

## Overview

This tutorial describes the one-vs-rest rhythmic program discovery
workflow in `scRhythmBayes`.

The one-vs-rest branch is designed for heterogeneous tissues and
multi-cell-type datasets. It compares one target cell type against all
remaining cell types and identifies cell-type-specific rhythmic
programs.

This branch starts from the gene-level rhythmic parameter layer. It is
different from pairwise rhythm-state module analysis:

| Analysis | Main question | Typical use case |
|----|----|----|
| Pairwise module analysis | How do rhythmic states differ between two predefined groups? | Condition `A` versus condition `B` within one cell type |
| One-vs-rest analysis | Which rhythmic programs are specific to one target cell type relative to all other cell types? | Cell-type-specific rhythmic program discovery in heterogeneous tissue |

In this tutorial, the workflow is:

1.  infer cell-level circadian phase;
2.  fit each target cell type and the pooled rest with sample-aware
    inference;
3.  inspect the OVR output table;
4.  visualize target-associated programs with a heatmap or fitted
    curves;
5.  extract target-specific candidate genes.

## Load package and example data

``` r
library(scRhythmBayes)

data("srb_example_data")

count <- srb_example_data$Count
metadata <- srb_example_data$meta
```

``` r
dim(count)
#> NULL
head(metadata)
#>                        cell    sample replicate ZT celltype condition
#> Cell_0001_001 Cell_0001_001 A_ZT00_R1         1  0        1         A
#> Cell_0001_002 Cell_0001_002 A_ZT00_R1         1  0        1         A
#> Cell_0001_003 Cell_0001_003 A_ZT00_R1         1  0        1         A
#> Cell_0001_004 Cell_0001_004 A_ZT00_R1         1  0        1         A
#> Cell_0001_005 Cell_0001_005 A_ZT00_R1         1  0        1         A
#> Cell_0001_006 Cell_0001_006 A_ZT00_R1         1  0        1         A
#>                   T_true Theta_true Theta_noisy     r_sim
#> Cell_0001_001 23.7745644 6.22416640   6.0883578 0.7240322
#> Cell_0001_002 23.7111107 6.20755425   6.0324486 0.6303893
#> Cell_0001_003 22.2182691 5.81672924   5.8263271 0.8666633
#> Cell_0001_004  0.5581621 0.14612650   0.2245110 0.8593562
#> Cell_0001_005  0.7203817 0.18859549   0.3327113 0.6901345
#> Cell_0001_006  0.1788832 0.04683152   0.1722017 0.7564754
table(metadata$celltype)
#> 
#>   1   2   3 
#> 720 720 720
table(metadata$celltype, metadata$condition)
#>    
#>       A   B
#>   1 360 360
#>   2 360 360
#>   3 360 360
```

The two required input objects are:

| Object | Description |
|----|----|
| `count` | Gene-by-cell count matrix |
| `metadata` | Cell metadata containing cell IDs, sample IDs, cell types, conditions, and sampling time |

For the example dataset, `metadata$celltype` defines the cell-type
groups used for OVR testing.

## Step 1: infer cell-level circadian phase

The OVR workflow first requires a cell-level circadian phase estimate.
Here, phase is inferred using `celltype` as the stratum.

``` r
res_phase <- infer_cell_phase(
  count = count,
  metadata = metadata,
  ZT_by = "ZT",
  strata_by = "celltype",
  Species = "mouse",
  setseed = 1,
  Verbose = FALSE
)

head(res_phase)
#> # A tibble: 6 × 7
#>   cell          ZT_input stratum theta_pred ZT_pred  Rbar r_cell
#>   <chr>            <dbl> <chr>        <dbl>   <dbl> <dbl>  <dbl>
#> 1 Cell_0001_001        0 1            6.19   23.6   0.992  0.999
#> 2 Cell_0001_002        0 1            6.23   23.8   0.991  0.999
#> 3 Cell_0001_003        0 1            6.00   22.9   0.991  0.999
#> 4 Cell_0001_004        0 1            6.14   23.4   0.989  0.999
#> 5 Cell_0001_005        0 1            0.142   0.543 0.992  0.999
#> 6 Cell_0001_006        0 1            0.137   0.524 0.991  0.999
```

The key output columns used downstream are:

| Column | Description |
|----|----|
| `cell` | Cell identifier |
| `theta_pred` | Inferred cell-level circadian phase in radians |
| `ZT_pred` | Inferred cell-level circadian phase converted to 24-hour time |
| `r_cell` | Cell-level phase reliability weight |

## Step 2: run one-vs-rest analysis

[`infer_ovr_rhythm()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_ovr_rhythm.md)
fits each target cell type and the pooled rest directly. Supplying
`sample_by = "sample"` uses independent biological samples for
rhythmicity testing. The complete three-target example takes a few
minutes; use `target_groups = "1"` for a shorter first run.

``` r
res_ovr <- infer_ovr_rhythm(
  count = count,
  metadata = metadata,
  theta_pred = res_phase$theta_pred,
  r_cell = res_phase$r_cell,
  group_by = "celltype",
  ZT_by = "ZT",
  target_groups = c("1", "2", "3"),
  q_cut = 0.05,
  amp_cut = 0.25,
  sample_by = "sample",
  sample_boot_by = "ZT",
  B = 3,
  n_worker = 1,
  setseed = 1,
  Verbose = FALSE
)

head(res_ovr$ovr_result)
res_ovr$ovr_summary
```

The main thresholds are:

| Argument | Description |
|----|----|
| `q_cut` | FDR cutoff for rhythmicity |
| `amp_cut` | Minimum target-group amplitude cutoff used for target-specific calls |

``` r
table(res_ovr$ovr_result$ovr_marker, useNA = "ifany")
```

## Interpret the OVR result

The OVR result is a gene-by-target-cell-type table.

| Column | Meaning |
|----|----|
| `gene` | Gene name |
| `celltype` | Target cell type being tested |
| `FDR_target`, `FDR_rest` | Sample-aware rhythmicity FDR in Target and Rest |
| `Mesor_target`, `Mesor_rest` | Fitted log-scale mesor values |
| `Mean_target`, `Mean_rest` | Fitted baseline means |
| `Amplitude_target`, `Amplitude_rest` | Fitted rhythmic amplitudes |
| `Phase_target`, `Phase_rest` | Fitted phases in radians |
| `Delta_Mean` | Target-minus-rest mean difference |
| `Delta_Amp` | Target-minus-rest amplitude difference |
| `ovr_marker` | Whether the gene is classified as target-specific |

A typical target-specific rhythmic program is expected to have strong
rhythmic evidence and sufficient amplitude in the target cell type, with
weaker or absent rhythmic evidence in the pooled rest.

## Step 4: plot the OVR heatmap

The OVR heatmap visualizes target-associated rhythmic programs across
the heterogeneous dataset.

``` r
ht_ovr <- plot_ovr_heatmap(
  count = count,
  res_phase = res_phase,
  ovr_result = res_ovr$ovr_result,
  metadata = metadata,
  celltype_by = "celltype",
  Verbose = FALSE
)

ht_ovr
plot(ht_ovr)
```

The heatmap is intended to help inspect whether OVR-selected genes show
cell-type-associated phase-resolved expression patterns.

## Step 5: plot fitted Target-versus-Rest curves

The curve function reproduces the original Figure 5F design from the
fitted harmonic parameters retained by
[`infer_ovr_rhythm()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_ovr_rhythm.md).
Each Target and Rest trajectory is normalized to mean one; no second
smoother is fitted.

``` r
gene_1 <- res_ovr$ovr_markers$gene[res_ovr$ovr_markers$celltype == "1"][1]

plot_ovr_curve(
  res_ovr,
  genes = gene_1,
  target_group = "1"
)
```

## Step 6: extract target-specific candidates

The `ovr_marker` column can be used to extract candidate
cell-type-specific rhythmic genes.

``` r
target_specific <- res_ovr$ovr_markers

head(target_specific)
```

The number of target-specific genes depends on signal strength, the
selected thresholds, and the structure of the dataset.

``` r
if (nrow(target_specific) > 0) {
  table(target_specific$celltype)
}
```

If no target-specific genes are detected in a small example dataset,
this does not indicate function failure. In real datasets, users may
adjust `q_cut` and `amp_cut` based on the expected signal and desired
stringency.

## Recommended minimal workflow

``` r
res_phase <- infer_cell_phase(
  count = count,
  metadata = metadata,
  ZT_by = "ZT",
  strata_by = "celltype",
  Species = "mouse",
  setseed = 1,
  Verbose = FALSE
)

res_ovr <- infer_ovr_rhythm(
  count = count,
  metadata = metadata,
  theta_pred = res_phase$theta_pred,
  r_cell = res_phase$r_cell,
  group_by = "celltype",
  ZT_by = "ZT",
  target_groups = "1",
  sample_by = "sample",
  sample_boot_by = "ZT",
  B = 3,
  n_worker = 1,
  setseed = 1,
  Verbose = FALSE
)

ht_ovr <- plot_ovr_heatmap(
  count = count,
  res_phase = res_phase,
  ovr_result = res_ovr$ovr_result,
  metadata = metadata,
  celltype_by = "celltype",
  Verbose = FALSE
)

plot(ht_ovr)
```

## Notes for real datasets

For real datasets, the OVR workflow is most useful when:

1.  the dataset contains multiple biologically meaningful cell types;
2.  cell types have sufficient cell numbers for rhythmic parameter
    estimation;
3.  the user wants cell-type-specific rhythmic programs rather than a
    predefined pairwise condition comparison;
4.  the target cell type may have rhythmic genes that are absent,
    weaker, or phase-shifted in the rest of the tissue.

For direct comparison of two selected groups, use
[`infer_rhythm_module()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_rhythm_module.md)
instead. For target-versus-all discovery across many cell types, use
[`infer_ovr_rhythm()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/infer_ovr_rhythm.md).

## Session information

``` r
sessionInfo()
#> R version 4.5.1 (2025-06-13)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Sequoia 15.5
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Asia/Shanghai
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] scRhythmBayes_2.0.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-6       gtable_0.3.6       jsonlite_2.0.0     dplyr_1.2.1       
#>  [5] compiler_4.5.1     Rcpp_1.1.1         tidyselect_1.2.1   jquerylib_0.1.4   
#>  [9] systemfonts_1.3.2  scales_1.4.0       textshaping_1.0.5  yaml_2.3.12       
#> [13] fastmap_1.2.0      lattice_0.22-9     ggplot2_4.0.3      R6_2.6.1          
#> [17] generics_0.1.4     knitr_1.51         MASS_7.3-65        htmlwidgets_1.6.4 
#> [21] ggrepel_0.9.8      tibble_3.3.1       desc_1.4.3         bslib_0.10.0      
#> [25] pillar_1.11.1      RColorBrewer_1.1-3 rlang_1.1.7        utf8_1.2.6        
#> [29] cachem_1.1.0       xfun_0.57          fs_1.6.7           sass_0.4.10       
#> [33] S7_0.2.1           otel_0.2.0         cli_3.6.5          pkgdown_2.2.0     
#> [37] magrittr_2.0.4     digest_0.6.39      grid_4.5.1         rstudioapi_0.18.0 
#> [41] lifecycle_1.0.5    vctrs_0.7.2        evaluate_1.0.5     glue_1.8.0        
#> [45] farver_2.1.2       ragg_1.5.2         rmarkdown_2.30     tools_4.5.1       
#> [49] pkgconfig_2.0.3    htmltools_0.5.9
```
