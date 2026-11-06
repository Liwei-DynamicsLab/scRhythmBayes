# Paths and working directory

Run all analysis code with `reproducibility/` as the working directory. The
scripts use relative paths and do not activate the author's local environment.

| Location | Contents |
|---|---|
| `../R/` | Package source in this repository |
| `01_cell_phase/` through `07_runtime/` | Included notebooks and scripts |
| `support/legacy/` | Included historical helper functions |
| `data/` | Input counts, metadata, seed sets and processed objects supplied by the reader |
| `data/gene_benchmark/` | Simulated gene-benchmark input objects |
| `outputs/simulation/` | Cell-phase simulation exports; historical subdirectory names retained |
| `outputs/gene_benchmark/` | Gene-benchmark cached results and figures |
| `outputs/liver/`, `outputs/liver_GSE145197/` | Liver cached results and figures |
| `outputs/heart/`, `outputs/drosophila/` | Biological analysis caches and outputs |
| `external/` | Comparator checkouts and additional helper scripts, where required |

Input files, caches, external repositories and generated results are ignored by
Git. Create the required directories and provide the files needed by the selected
route. A path in a script does not imply that its target is included.

## R scripts and R Markdown

Open R with this directory as the working directory, then source an R script.
For rendering a notebook, explicitly keep the same root for its chunks:

```r
analysis_root <- normalizePath(".")
stopifnot(file.exists("PATHS.md"), dir.exists("01_cell_phase"))
rmarkdown::render("05_heart/scHeart_complete_reanalysis_20261022.Rmd",
                  knit_root_dir = analysis_root)
```

Rendering honors each notebook's existing `eval` flags. This does not execute
disabled analysis chunks. For interactive chunks, configure the IDE's chunk
working directory to this folder.

## Python scripts and Jupyter

Launch Python scripts from this directory. For a notebook opened in
`01_cell_phase/`, run this setup cell before the supplied analysis cells:

```python
from pathlib import Path
import os

here = Path.cwd().resolve()
root = next((p for p in (here, *here.parents)
             if (p / "PATHS.md").is_file() and (p / "01_cell_phase").is_dir()), None)
if root is None:
    raise RuntimeError("Open the notebook inside the reproducibility directory")
os.chdir(root)
```

Install comparator dependencies in the chosen kernel. The code no longer tests
for a particular personal conda executable. Historical source-version checks
remain in place where they are part of a comparator runner.
