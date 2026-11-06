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
- R analyses: install the libraries imported in the selected notebook into the environment used to run it.
- Fly preprocessing and downstream sections require separate Seurat versions as documented in the notebook. Do not mix their assumptions in one session.
- Python comparators use method-specific environments. Jupyter source is preserved, while execution outputs were removed to avoid bundling stale logs or embedded data.
- The original notebook seeds and thresholds are retained. Historical full-session package versions were not reconstructed or fabricated.
- `capture_environment.R` records the environment of a future run; its output must not be represented as the historical analysis environment.

## Path configuration and additional files

See `PATHS.md` for the relative directory layout. Use the R/Python environment
appropriate for the selected analysis; scripts no longer activate a personal
renv project, append a machine-specific R library, or require a fixed Python
executable. The Seurat version requirements in the analysis remain applicable.

The `support/legacy` functions are distinct historical implementations. The
files under `external/inference_original` and `external/seed_revision` are
additional source versions that are not bundled; a same-named function in the
package is not necessarily interchangeable.

Some plotting sections require external layout utilities (`panel_alignment.R`,
`audit_panel_alignment.py`), `sc_gene_curve.R`, figure-assembly scripts, layout
images or cached plot inputs. These are not bundled. The heart setup currently
requires the layout utilities. Review file-existence checks and `source()` calls
before choosing a route. Missing dependencies are not removed or silently
replaced by unrelated functions.

This source collection does not supply complete raw-data-to-figure execution
for every historical route. Point benchmark input paths to the corresponding
simulator outputs. Runtime summaries additionally require measured metrics and
provenance files described in `07_runtime/README.md`.
