# scRhythmBayes manuscript analysis code

Analysis notebooks, benchmark runners and helper functions accompanying the
scRhythmBayes manuscript. Install the package from the parent repository; its
source and dependency declarations are in `../R/` and `../DESCRIPTION`.

## Analyses

| Directory | Main entry points | Analysis |
|---|---|---|
| `01_cell_phase` | `scRITMO_simulation_6CT.ipynb`, method notebooks, `6CT_all_benchmark_cell_phase.Rmd` | Cell-phase simulation and method comparison |
| `02_gene_benchmark` | `simulated_data_sample_specific.Rmd`, `all_gene_level_benchmark.Rmd` | Sample-aware gene benchmark |
| `03_liver` | `scCircadianLiver_new.Rmd` | Main liver analysis and appended GSE145197 analysis |
| `04_ovr` | `scliver_ovr_run_example_updated.Rmd` | One-versus-rest, leave-one-cell-type-out and pairwise comparisons |
| `05_heart` | `scHeart_complete_reanalysis_20261022.Rmd`, `scHeart_plot_20261023.Rmd` | Cardiac phase organization and plots |
| `06_drosophila` | `GSE157504_LD_scRhythmBayes_validation.Rmd` | Drosophila analysis |
| `07_runtime` | See its `README.md` | Runtime and memory profiling |
| `support/legacy` | Functions explicitly used by historical analysis routes | Analysis-specific helper implementations |

## Running an analysis

1. Read `DATA_AND_ENVIRONMENT.md` for input data and environment requirements.
2. Configure file paths in the selected notebook. `/path/to/...` strings are
   explicit placeholders: replace them with your data, repository and output
   locations. They are not directories supplied by this repository. Home-relative
   paths (`~/...`) also need review. Use a separate output directory for a new run.
3. Install the packages imported by that notebook. Python comparators and Seurat
   v4/v5 routes require their respective environments.
4. Run the relevant setup, inference and plotting chunks in order. Several
   notebooks start from saved results, contain `eval=FALSE` chunks, or preserve
   alternative plotting sections. Select the intended analysis route before
   running; these are analysis notebooks, not a single-command pipeline.
5. Run `capture_environment.R` in the R session used for the analysis to record
   its package versions. Save comparator versions and input provenance alongside
   new outputs.

Notebook panel names include historical layouts; use the manuscript legends for
final panel numbering. Enrichment sections use Reactome.

## Included material and validation

Large data objects, fitted-result caches and historical environments are not
included. Some plotting routes and historical inference variants require
additional external files identified in `DATA_AND_ENVIRONMENT.md`.
R and Python sources have been checked for syntax after path anonymization;
full analyses and clean installation have not been rerun for this collection.
`source_manifest.tsv` lists distributed files and SHA-256 checksums.
