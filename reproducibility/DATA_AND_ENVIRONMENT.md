# Inputs and software

Large `.qs`, `.rds`, count matrices and raw patient-level data are not bundled. GEO accessions identify source studies; they do not imply that an author's processed object can be downloaded directly from GEO with identical contents.

| Analysis | Source accession / required object | Preparation notes |
|---|---|---|
| Empirical phase simulation | GSE145197; empirical parameters/library sizes and simulation export | Run the simulation notebook and method notebooks with consistent cells and phase references. |
| Gene benchmark | `simulation_data_replicates.qs`; benchmark caches | The included simulator generates the 36-sample design. Its save paths and the benchmark's load paths differ in the original sources and must be aligned. |
| Main liver / OVR | GSE288579; multiome files / processed RNA Seurat object | Preserve sample/ZT labels and Hep.CV → Hep.PV comparison order. The liver notebook contains both preprocessing and saved-object routes. |
| Independent liver | GSE145197; processed liver object | The appended notebook section imports an existing object. A complete provenance chain for that processed object is required for raw-data reproduction. |
| Heart | GSE270258; processed cardiac object, saved phases, sample labels | NCR/ACR × AM/PM, 16 samples. Representative 8/16 h values encode collection windows rather than exact measured ZT. Cached phases/plots are not independent validation. |
| Drosophila | GSE157504; raw counts, original cluster annotations and mapping | Follow the notebook's LD/DD handling and 17 retained neuronal clusters. |

Official data entry points:

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288579
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145197
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270258
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157504

These links are accession entry points, not verified direct-download manifests. File checksums and exact processed-object preparation steps should accompany an archived data release.

## Environment boundaries

- Package source and dependencies: see the parent repository (`../R/`, `../DESCRIPTION`). Historical cached results may have been produced with different revisions.
- R analyses: install the libraries imported in the selected notebook. Configure `.libPaths()` and environment activation statements for your own installation.
- Fly preprocessing and downstream sections require separate Seurat versions as documented in the notebook. Do not mix their assumptions in one session.
- Python comparators use method-specific environments. Jupyter source is preserved, while execution outputs were removed to avoid bundling stale logs or embedded data.
- The original notebook seeds and thresholds are retained. Historical full-session package versions were not reconstructed or fabricated.
- `capture_environment.R` records the environment of a future run; its output must not be represented as the historical analysis environment.

## Path configuration and additional files

Paths have been anonymized without changing model parameters, seeds or
statistical thresholds. Replace these roots consistently within the selected
analysis and its helper scripts:

| Placeholder | Required location |
|---|---|
| `/path/to/scRhythmBayes` | This repository checkout |
| `/path/to/analysis_home` | Your analysis workspace containing processed data, simulation exports, comparator checkouts and environments; preserve or edit the suffix for each input |
| `/path/to/data_volume` | Data directory containing the study input objects |
| `/path/to/source_archive` | Historical analysis inputs/source archive, if that route is used |
| `/path/to/analysis_workspace` | Additional historical analysis files referenced by a selected route |
| `/path/to/legacy_study` | Inputs for legacy Windows analysis sections |
| `/path/to/external_plot_helpers` | External `panel_alignment.R` / `audit_panel_alignment.py` layout utilities, not bundled |

The included `support/legacy` functions are distinct historical implementations;
do not replace them with current package functions without checking the calling
interface. Historical `infer_cell_phase.R` files under `outputs/seed_revision`
or `outputs/` are not supplied merely because another function with the same
name is included. Those routes require the matching source version.

Some figure sections additionally reference `sc_gene_curve.R`, an existing fly
figure-assembly script, layout images, and precomputed plot inputs. These are
not bundled. Review `source()`, input-loading calls and file-existence checks
before selecting those sections. The heart notebook currently requires the
external panel-layout utilities even for its setup chunk. This source collection
therefore does not provide complete raw-data-to-figure execution for every route.

The simulator and benchmark may use different cache filenames/directories;
point the benchmark to the simulator output explicitly. Runtime summary scripts
also require measured metrics and provenance files; see `07_runtime/README.md`.
