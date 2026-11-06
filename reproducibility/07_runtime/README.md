# Computational-resource benchmark code

`runners/` contains method-specific CellPhase runners. `scripts/` contains the
process profiler, gene-level runners, summary-table builder and plotting code.
Configure their input, output and environment paths before execution.

1. Prepare the shared simulation input and method-specific environments.
2. Use `python scripts/profile_process.py --help` to inspect profiler arguments,
   then profile each runner in its required environment.
3. Retain the resulting metrics JSON, sampled resource CSV and prediction outputs.
4. The summary-table builder expects `metrics/`, `runs/` and per-runner
   `*.provenance.json` records under this directory. These historical measurement
   files are not included. Supply the corresponding records from a measured run
   before running the builder. The plotting script requires the summary tables
   and external layout utilities.

Wall time covers runner launch through result export, excluding installation.
Peak aggregate RSS sums resident memory over the monitored processes; shared
pages may be counted more than once. Effective CPU cores are summed process CPU
percentages divided by 100. Results depend on hardware, environments and native
method workloads. The code alone does not establish a runtime ranking.
