# scRhythmBayes <img src="man/figures/logo.png" align="right" width="160">

[![R-CMD-check](https://github.com/Liwei-DynamicsLab/scRhythmBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Liwei-DynamicsLab/scRhythmBayes/actions/workflows/R-CMD-check.yaml)

**Cell-level circadian phase inference and rhythm-state remodeling for single-cell transcriptomes**

`scRhythmBayes` is an R package for phase-resolved analysis of heterogeneous
single-cell transcriptomes. It infers latent cell-level circadian phase,
estimates gene-level rhythmic parameters, performs biological-sample-aware
rhythmicity testing, and classifies differential rhythm states between groups.

<p align="center">
  <img src="man/figures/readme.png" width="760">
</p>

## Analysis map

| Layer | Main question | Function |
|---|---|---|
| Cell phase | What is the latent circadian phase of each cell? | `infer_cell_phase()` |
| Gene rhythm | Is a gene rhythmic, and what are its mesor, amplitude, and phase? | `infer_gene_param()` |
| Two-group remodeling | Does rhythmic state change between two predefined groups? | `infer_rhythm_module()` |
| One-vs-rest | Is a rhythmic program specific to one target cell type? | `infer_ovr_rhythm()` |
| Phase distribution | Are cell-phase distributions similar across groups? | `run_dtw_similarity()` |
| Phase concentration | Are phases concentrated or dispersed within samples/groups? | `run_phase_rho()` |

The detailed, unified tutorial is available at
[`vignettes/complete-workflow.Rmd`](vignettes/complete-workflow.Rmd). It
covers input validation, sample-aware inference, M0--M7 interpretation, OVR
curves, DTW, rho, simulation-truth checks, and troubleshooting.

## Installation

```r
install.packages("remotes")
remotes::install_github("Liwei-DynamicsLab/scRhythmBayes")
```

```r
library(scRhythmBayes)
packageVersion("scRhythmBayes")
```

If the repository is private, GitHub authentication and repository access are
required.

## Built-in teaching data

```r
data("srb_example_data")

names(srb_example_data)
dim(srb_example_data$Count)
dim(srb_example_data$meta)
head(srb_example_data$meta)
```

The current `srb_example_data` contains:

- 138 genes and 2,160 cells;
- 36 independent biological samples;
- two conditions, six ZT values, and three replicates per condition-by-ZT cell;
- 60 cells per sample and three cell types;
- negative-binomial counts with mean-dependent zero inflation;
- native sample-specific mesor and phase effects;
- 24 mouse phase-seed genes;
- M0--M7, OVR, cell-phase, and sample-effect truth.

The object has three components:

| Component | Description |
|---|---|
| `Count` | Sparse gene-by-cell raw count matrix |
| `meta` | Cell metadata and simulation phase truth |
| `truth` | Gene truth, sample design, sample effects, and simulation settings |

Important metadata columns:

| Column | Description |
|---|---|
| `cell` | Cell ID matched to `colnames(Count)` |
| `sample` | Independent biological sample |
| `replicate` | Replicate within condition and ZT |
| `ZT` | Sampling Zeitgeber time |
| `celltype` | Cell-type stratum |
| `condition` | Condition A or B |
| `T_true`, `Theta_true` | Simulated true phase in hours and radians |
| `Theta_noisy`, `r_sim` | Simulated noisy phase inputs retained as truth |

`Theta_noisy` and `r_sim` are not outputs from phase inference. The
end-to-end workflow below uses `theta_pred` and `r_cell` inferred from
the count matrix.

Verify the sample design:

```r
count <- srb_example_data$Count
metadata <- srb_example_data$meta

stopifnot(identical(colnames(count), metadata$cell))

sample_design <- unique(
  metadata[c("sample", "condition", "ZT", "replicate")]
)

nrow(sample_design)
with(sample_design, table(condition, ZT))
table(metadata$sample)
```

## Core workflow

### 1. Infer cell-level circadian phase

```r
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
```

The main output columns are:

| Column | Meaning |
|---|---|
| `theta_pred` | Inferred phase in radians |
| `ZT_pred` | Inferred phase in hours |
| `Rbar` | Posterior circular concentration |
| `r_cell` | Cell-level rhythmic reliability used downstream |

Inspect the exact seed-gene selection:

```r
attr(res_phase, "n_seed_matched")
attr(res_phase, "n_seed_retained")
attr(res_phase, "seed_genes")
head(attr(res_phase, "seed_diagnostics"))
```

Basic diagnostic plots:

```r
plot_phase_scatter(
  res_phase,
  metadata,
  truth_by = "T_true",
  pred_by = "ZT_pred"
)$plot

plot_phase_circle(
  res_phase,
  metadata,
  zt_by = "ZT",
  group_by = "celltype"
)$plot

plot_phase_rose(
  res_phase,
  metadata,
  group_by = "condition",
  split_by = "celltype"
)$plot
```

The scatter plot uses simulated truth. Circle, flower, and rose plots can be
used with real datasets that do not contain true phase.

### 2. Run sample-aware gene-level inference

For a condition comparison, select one cell type and keep its biological
sample IDs.

```r
cell_use <- metadata$celltype == "1"

count_use <- count[, cell_use, drop = FALSE]
metadata_use <- metadata[cell_use, , drop = FALSE]

length(unique(metadata_use$sample))
table(metadata_use$condition)
```

```r
res_param <- infer_gene_param(
  count = count_use,
  metadata = metadata_use,
  theta_pred = res_phase$theta_pred[cell_use],
  r_cell = res_phase$r_cell[cell_use],
  strata_by = "condition",
  ZT_by = "ZT",
  sample_by = "sample",
  sample_boot_by = "ZT",
  B = 3,
  n_worker = 1,
  setseed = 1,
  Verbose = FALSE
)

head(res_param$param_long)
```

With `sample_by`, P values and FDR are based on independent biological
samples through the restricted-null cluster multiplier score test. Multiple
cells from the same sample are not treated as independent replicates.

Check that sample-aware inference was actually used:

```r
res_param$inference

table(
  res_param$fit_diagnostics$Sample_status,
  useNA = "ifany"
)

stopifnot(
  identical(
    res_param$inference$mode,
    "restricted_null_cluster_multiplier_score"
  )
)
```

`B = 3` keeps the teaching example fast. It controls bootstrap confidence
intervals, not the sample-aware score-test P values. Increase `B` for final
confidence intervals.

Main outputs:

| Object | Description |
|---|---|
| `param_long` | Gene parameters by condition |
| `param_wide` | Paired wide table for module analysis |
| `fit_diagnostics` | Convergence and sample-score diagnostics |
| `bootstrap_res` | Confidence intervals for top-ranked genes |
| `gene_eligibility` | Time-resolved expression eligibility |
| `inference` | Statistical mode, sample columns, and sample counts |

Important parameters:

| Column | Description |
|---|---|
| `Mesor` | Log-scale harmonic intercept |
| `Mean` | Exponentiated baseline mean |
| `Amplitude` | Fitted harmonic amplitude |
| `Amplitude_shrunk` | Stabilized amplitude effect size |
| `Phase` | Peak phase in radians |
| `P_value` | P value from the active inference mode |
| `FDR` | BH-adjusted rhythmicity significance |

### 3. Infer M0--M7 rhythm-state modules

```r
res_module <- infer_rhythm_module(
  param_wide = res_param$param_wide,
  compare_by = "condition",
  group1 = "A",
  group2 = "B",
  tau_q = 0.05,
  tau_A = 0.25,
  tau_dA = 1,
  tau_dphi = pi / 6,
  Verbose = FALSE
)

res_module$comparison
res_module$module_summary
```

All labels describe `group2` relative to `group1`.

| ID | Module | Interpretation of B relative to A |
|---:|---|---|
| M0 | Non-rhythmic | Neither group has supported rhythmicity |
| M1 | Stable Rhythm | Rhythmic in both groups without a large change |
| M2 | Rhythm Gain | Rhythmicity gained in B |
| M3 | Rhythm Loss | Rhythmicity lost in B |
| M4 | Amplitude-enhanced rhythmic | Amplitude increased in B |
| M5 | Amplitude-attenuated rhythmic | Amplitude decreased in B |
| M6 | Phase-advanced rhythmic | Phase is earlier in B |
| M7 | Phase-delayed rhythmic | Phase is later in B |

`seed_module` is the rule-based initial class; `module_map` is the final
posterior assignment. `post_max` and `margin` summarize assignment
certainty.

### 4. Visualize rhythm-state modules

```r
hm_module <- plot_module_heatmap(
  count = count_use,
  res_phase = res_phase[cell_use, , drop = FALSE],
  gene_meta = res_module$module_result,
  metadata = metadata_use,
  group_by = "condition",
  group1 = "A",
  group2 = "B",
  Verbose = FALSE
)

plot(hm_module, module_id = 2)
```

```r
plot_module_polar(
  gene_meta = res_module$module_result,
  module_id = 7,
  group1 = "A",
  group2 = "B"
)$plot
```

## One-vs-rest rhythmic programs

`infer_ovr_rhythm()` directly refits Target and pooled Rest groups. It does
not take a separately fitted `param_long` table.

```r
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

res_ovr$ovr_summary
head(res_ovr$ovr_markers)
```

A gene is called an OVR marker when the target is rhythmic and sufficiently
strong, while the pooled Rest is not significant:

- `FDR_target < q_cut`;
- `Amplitude_target >= amp_cut`;
- `FDR_rest >= q_cut`.

For a faster first check, use `target_groups = "1"`.

### Model-fitted Target-versus-Rest curves

```r
genes_1 <- head(
  res_ovr$ovr_markers$gene[
    res_ovr$ovr_markers$celltype == "1"
  ],
  4
)

plot_ovr_curve(
  res_ovr,
  genes = genes_1,
  target_group = "1",
  ncol = 2
)
```

The curves use the fitted model
`exp(Mesor + Amplitude * cos(theta - Phase))`; no second smoother is fitted.

### OVR heatmap

```r
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

## Phase-distribution extensions

### DTW similarity

```r
dtw_res <- run_dtw_similarity(
  res_phase = res_phase,
  metadata = metadata,
  celltype_by = "celltype",
  group_by = "condition",
  n_bin = 24,
  min_cell = 10,
  normalize_distance = TRUE,
  Verbose = FALSE
)

dtw_res$similarity

ht_dtw <- plot_dtw_heatmap(
  dtw_res,
  matrix_type = "similarity"
)

plot(ht_dtw)
```

DTW compares complete circular phase-distribution shapes. It is descriptive;
it is not a biological-replicate hypothesis test.

### Sample-level circular rho

```r
rho_res <- run_phase_rho(
  res_phase = res_phase,
  metadata = metadata,
  celltype_by = "celltype",
  group_by = "condition",
  sample_by = "sample",
  theta_by = "theta_pred",
  group_levels = c("A", "B"),
  compare_groups = c("A", "B"),
  min_cell = 10,
  Verbose = FALSE
)

head(rho_res$rho_sample)
rho_res$rho_sample_summary
rho_res$sig_result

plot_phase_rho(rho_res)$plot
```

With `sample_by = "sample"`, between-condition comparisons use sample-level
rho values instead of individual cells as replicates.

## Using a real dataset

Minimum metadata:

| Required field | Example |
|---|---|
| Cell ID | `cell` |
| Collection time | `ZT` |
| Phase stratum | `celltype` |
| Comparison group | `condition` |
| Independent sample | `sample` |

Recommended checks:

```r
stopifnot(
  !anyDuplicated(rownames(count)),
  !anyDuplicated(colnames(count)),
  !anyDuplicated(metadata$cell),
  identical(colnames(count), metadata$cell),
  all(c("ZT", "celltype", "condition", "sample") %in% names(metadata))
)
```

Do not randomly assign sample IDs to cells. Sample IDs must represent the
original animal, patient, culture, or other independent experimental unit.

## Runtime and validation

On the bundled 138-by-2,160 example data, the tested teaching settings took
approximately:

| Step | Runtime |
|---|---:|
| Cell-phase inference | 4 seconds |
| Two-condition sample-aware gene inference | 35 seconds |
| Module inference | less than 1 second |
| Full three-target sample-aware OVR | 3.4 minutes |

The bundled simulation recovered all eight M0--M7 categories. In the recorded
validation run, all sample-aware tests reported
`ok_cluster_multiplier_score`, and OVR returned 17 target-matched markers.

See
[`data-raw/validation_report.md`](data-raw/validation_report.md) and
[`data-raw/validate_srb_example_data.R`](data-raw/validate_srb_example_data.R)
for the exact validation checks.

## Troubleshooting

| Problem | What to inspect |
|---|---|
| Count/metadata mismatch | `identical(colnames(count), metadata$cell)` |
| No phase seeds | Species, gene symbols, aliases, and row names |
| Too few retained seeds | `attr(res_phase, "seed_diagnostics")` |
| Too few independent samples | Experimental design; cells cannot replace samples |
| Rank-deficient phase model | Circadian coverage within each fitting stratum |
| Sample-aware fallback suspected | `res_param$inference` and `Sample_status` |
| No OVR markers | Target/rest FDR and amplitude before changing thresholds |
| Reversed module interpretation | `group1` and `group2` direction |
| Slow run | Test one cell type or OVR target with `B = 3` first |

## Tutorials

- [Complete workflow](vignettes/complete-workflow.Rmd)
- [Basic workflow](vignettes/basic-workflow.Rmd)
- [Rhythm-state modules](vignettes/module-analysis.Rmd)
- [One-vs-rest analysis](vignettes/ovr-analysis.Rmd)
- [DTW and rho phase prioritization](vignettes/phase-prioritization-dtw-rho.Rmd)

## Main functions

| Function | Description |
|---|---|
| `infer_cell_phase()` | Infer cell-level circadian phase |
| `infer_gene_param()` | Estimate gene parameters and run sample-aware tests |
| `infer_rhythm_module()` | Infer M0--M7 differential rhythm states |
| `plot_module_heatmap()` | Plot phase-resolved module heatmaps |
| `plot_module_polar()` | Plot module genes in amplitude-phase space |
| `infer_ovr_rhythm()` | Identify one-vs-rest rhythmic markers |
| `plot_ovr_curve()` | Plot fitted Target-versus-Rest trajectories |
| `plot_ovr_heatmap()` | Plot OVR rhythmic programs |
| `run_dtw_similarity()` | Compare phase-distribution shape |
| `plot_dtw_heatmap()` | Plot DTW distance or similarity |
| `run_phase_rho()` | Compute group- or sample-level circular concentration |
| `plot_phase_rho()` | Plot circular rho summaries |

## Citation

If you use `scRhythmBayes`, please cite the associated manuscript:

Zhang Q, Li C, Zhang L. *Decoding cell-level circadian phase and rhythm-state
remodeling in single-cell transcriptomes with scRhythmBayes.*

## License

MIT License.

## Authors

Liwei Zhang, Qian Zhang, and Chenyu Li.
