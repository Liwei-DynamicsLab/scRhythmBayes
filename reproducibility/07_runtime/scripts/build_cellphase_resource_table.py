#!/usr/bin/env python3
"""Build the audited CellPhase resource table from measured JSON and CSV files."""

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

METHODS = [
    {
        "method": "scRhythmBayes",
        "fit_scope": "by cell type",
        "metric": "cellphase_scRhythmBayes.json",
        "output": "runs/scRhythmBayes/scRhythmBayes_cell_phase.csv",
        "prediction": "ZT_pred",
        "provenance": "run_scRhythmBayes.provenance.json",
        "model_input_note": "1000-gene input; current infer_cell_phase implementation",
    },
    {
        "method": "scRITMO",
        "fit_scope": "global",
        "metric": "cellphase_scRITMO.json",
        "output": "runs/scRITMO/scRITMO_cell_phase.csv",
        "prediction": "pred_hour",
        "provenance": "run_scRITMO.provenance.json",
        "model_input_note": "1000-gene input; CPU torch path",
    },
    {
        "method": "Tempo",
        "fit_scope": "global",
        "metric": "cellphase_Tempo.json",
        "output": "runs/Tempo/Tempo_cell_phase.csv",
        "prediction": "pred_hour",
        "provenance": "run_Tempo.provenance.json",
        "model_input_note": "1000-gene input; notebook final solution",
    },
    {
        "method": "CoPhaser",
        "fit_scope": "global",
        "metric": "cellphase_CoPhaser.json",
        "output": "runs/CoPhaser/CoPhaser_cell_phase.csv",
        "prediction": "pred_hour",
        "provenance": "run_CoPhaser.provenance.json",
        "model_input_note": "1000-gene input; 200 epochs with original restart logic",
    },
    {
        "method": "CYCLOPS",
        "fit_scope": "global",
        "metric": "cellphase_CYCLOPS.json",
        "output": "runs/CYCLOPS/CYCLOPS_cell_phase.csv",
        "prediction": "pred_hour",
        "provenance": "run_CYCLOPS.provenance.json",
        "model_input_note": "1000-gene input; 17 retained seeds, 5 eigengenes, 80 models",
    },
    {
        "method": "reCAT",
        "fit_scope": "by cell type",
        "metric": "cellphase_reCAT.json",
        "output": "runs/reCAT/reCAT_6CT_cell_phase.csv",
        "prediction": "reCAT_phase_hour",
        "provenance": "run_reCAT.provenance.json",
        "model_input_note": "1000-gene input; 15 predefined clock genes passed to reCAT",
    },
    {
        "method": "Cyclum",
        "fit_scope": "by cell type",
        "metric": "cellphase_Cyclum.json",
        "output": "runs/Cyclum/Cyclum_cell_phase.csv",
        "prediction": "pred_hour",
        "provenance": "run_Cyclum.provenance.json",
        "model_input_note": "1000-gene input; six cell-type fits",
    },
    {
        "method": "scPrisma",
        "fit_scope": "by cell type",
        "metric": "cellphase_scPrisma.json",
        "output": "runs/scPrisma/scPrisma_cell_phase.csv",
        "prediction": "pred_hour",
        "provenance": "run_scPrisma.provenance.json",
        "model_input_note": "1000-gene input; CPU/Numba path, six cell-type fits",
    },
]


def finite_number(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def source_from_provenance(record: dict) -> str:
    path = ROOT / "runners" / record["provenance"]
    provenance = json.loads(path.read_text(encoding="utf-8"))
    return provenance.get("executable_source") or provenance.get("source")


def validate_output(record: dict) -> tuple[int, int, int]:
    output = ROOT / record["output"]
    with output.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "cell" not in rows[0] or record["prediction"] not in rows[0]:
        raise RuntimeError(f"Missing validation columns in {output}")
    unique_cells = len({row["cell"] for row in rows})
    predicted = sum(finite_number(row[record["prediction"]]) for row in rows)
    return len(rows), unique_cells, predicted


def main() -> None:
    output_rows = []
    for record in METHODS:
        metric_path = ROOT / "metrics" / record["metric"]
        metric = json.loads(metric_path.read_text(encoding="utf-8"))
        n_rows, n_unique, n_predicted = validate_output(record)
        validation = (
            "pass"
            if metric["status"] == "complete" and n_rows == 9000 and n_unique == 9000 and n_predicted > 0
            else "fail"
        )
        output_rows.append(
            {
                "method": record["method"],
                "status": metric["status"],
                "validation_status": validation,
                "fit_scope": record["fit_scope"],
                "n_cells_input": 9000,
                "n_genes_input": 1000,
                "n_output_rows": n_rows,
                "n_unique_cells": n_unique,
                "n_predicted_phase": n_predicted,
                "prediction_coverage": n_predicted / n_rows,
                "runtime_minutes": metric["wall_minutes"],
                "peak_rss_gb": metric["peak_process_tree_rss_gb"],
                "configured_cores": metric.get("configured_cores"),
                "mean_effective_cpu_cores": metric["mean_effective_cpu_cores"],
                "peak_effective_cpu_cores": metric["peak_effective_cpu_cores"],
                "sample_interval_seconds": metric["sample_interval_seconds"],
                "resource_scope": "complete isolated runner; setup/package installation excluded",
                "rss_definition": metric["rss_definition"],
                "cpu_definition": metric["cpu_definition"],
                "model_input_note": record["model_input_note"],
                "source_file": source_from_provenance(record),
                "runner_file": str(ROOT / "runners" / f"run_{record['method']}.{'R' if record['method'] in {'scRhythmBayes', 'reCAT', 'CYCLOPS'} else 'py'}"),
                "metric_file": str(metric_path),
                "result_file": str(ROOT / record["output"]),
            }
        )

    if any(row["validation_status"] != "pass" for row in output_rows):
        raise RuntimeError("At least one CellPhase output failed validation")

    out = ROOT / "data" / "cellphase_resource_current_evidence.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    print(out)


if __name__ == "__main__":
    main()
