# Generated from the named notebook; original file is unchanged.
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# ---- notebook code cell 0 ----
# 01 · Imports
from pathlib import Path
import inspect
import json
import hashlib
import random
import shutil
import numpy as np
import pandas as pd
import anndata as ad
import torch
from scipy import sparse
import tempo
from tempo import unsupervised_alg, cell_posterior

print("Tempo:", tempo.__file__)

# ---- notebook code cell 1 ----
# 02 · Paths
RUN_DIR = Path("outputs/simulation/output/scRITMO_extended_simulation/sigma_gradient_1p5_2p5_3p5_4p5_5p5_6p5_empirical_fixed_v1")
H5AD_PATH = RUN_DIR / "scRITMO_extended_simulation.h5ad"
METADATA_PATH = RUN_DIR / "R_input/metadata.csv"
GENE_PARAMETERS_PATH = RUN_DIR / "baseline_gene_parameters.csv"
GENE_PRIOR_PATH = Path("data/gene_acrophase_prior.csv")
CORE_CLOCK_PATH = Path("data/Tempo_core_clock_genes.txt")
OUT_DIR = Path("07_runtime/runs/Tempo")
reference_gene = "Bmal1"
SEED = 1

for path in [H5AD_PATH, METADATA_PATH, GENE_PARAMETERS_PATH, GENE_PRIOR_PATH, CORE_CLOCK_PATH]:
    assert path.is_file(), str(path)

# ---- notebook code cell 2 ----
# 03 · Raw counts
source = ad.read_h5ad(H5AD_PATH)
assert source.shape == (9000, 1000), source.shape
assert source.obs_names.is_unique and source.var_names.is_unique

cell_ids = source.obs_names.copy()
gene_ids = source.var_names.copy()
raw_counts = source.X.copy()
values = raw_counts.data if sparse.issparse(raw_counts) else np.asarray(raw_counts)

assert np.isfinite(values).all()
assert (values >= 0).all()
assert np.equal(values, np.floor(values)).all()
assert (np.asarray(raw_counts.sum(axis=1)).ravel() > 0).all()

adata = ad.AnnData(
    X=raw_counts,
    obs=pd.DataFrame(index=cell_ids),
    var=pd.DataFrame(index=gene_ids),
)

assert adata.obs_names.equals(cell_ids)
assert adata.var_names.equals(gene_ids)

print("Raw counts:", adata.shape)

# ---- notebook code cell 3 ----
# 04 · Clock priors
prior = pd.read_csv(GENE_PRIOR_PATH, index_col="gene")
core_genes = [x for x in CORE_CLOCK_PATH.read_text().splitlines() if x]

expected_clock = {
    "Bhlhe40", "Bhlhe41", "Bmal1", "Ciart", "Clock", "Cry1", "Dbp",
    "Npas2", "Nr1d1", "Nr1d2", "Per1", "Per2", "Per3", "Rorc", "Tef",
}

assert len(core_genes) == 15
assert len(set(core_genes)) == 15
assert set(core_genes) == expected_clock
assert prior.index.is_unique
assert len(prior) == 15
assert set(prior.index) == expected_clock
assert list(prior.columns) == ["prior_acrophase_loc", "prior_acrophase_95_interval"]
assert expected_clock.issubset(adata.var_names)
assert np.isfinite(prior.to_numpy()).all()
assert prior["prior_acrophase_loc"].between(0, 2 * np.pi, inclusive="left").all()
assert prior["prior_acrophase_95_interval"].between(1e-3, np.pi).all()
assert np.isclose(prior.loc[reference_gene, "prior_acrophase_loc"], 0.0)
assert (np.asarray(adata[:, core_genes].X.sum(axis=0)).ravel() > 0).all()

print("Clock genes:", len(core_genes), "/ 15")

# ---- notebook code cell 4 ----
# 05 · Fixed reference rotation
parameters = pd.read_csv(GENE_PARAMETERS_PATH, index_col="gene_id")
assert parameters.index.equals(gene_ids)

clock = parameters.loc[core_genes]
gene_phase = np.arctan2(clock["b_1"], clock["a_1"]) % (2 * np.pi)

BMAL1_SIMULATION_PHASE_RAD = float(gene_phase.loc[reference_gene])
TEMPO_TO_SIMULATION_ROTATION_RAD = BMAL1_SIMULATION_PHASE_RAD

expected_prior = (gene_phase - BMAL1_SIMULATION_PHASE_RAD) % (2 * np.pi)
prior_difference = (
    prior.loc[core_genes, "prior_acrophase_loc"]
    - expected_prior
    + np.pi
) % (2 * np.pi) - np.pi

assert np.max(np.abs(prior_difference)) < 1e-5

rotation_record = {
    "reference_gene": reference_gene,
    "simulation_reference_rad": BMAL1_SIMULATION_PHASE_RAD,
    "tempo_reference_rad": 0.0,
    "rotation_rad": TEMPO_TO_SIMULATION_ROTATION_RAD,
    "formula": "pred_rad = (tempo_map_rad + rotation_rad) % (2*pi)",
    "alignment_uses_cell_truth": False,
    "gene_prior_source": "simulation generator clock phases",
}

print(json.dumps(rotation_record, indent=2))

# ---- notebook code cell 5 ----
# 06 · Fresh official Tempo run
signature = inspect.signature(unsupervised_alg.run)
assert signature.parameters["use_noninformative_phase_prior"].default is True
assert signature.parameters["cell_phase_prior_path"].default is None

random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

if OUT_DIR.exists():
    shutil.rmtree(OUT_DIR)

OUT_DIR.mkdir(parents=True, exist_ok=False)

(OUT_DIR / "phase_rotation.json").write_text(
    json.dumps(rotation_record, indent=2)
)

provenance = {
    "h5ad": str(H5AD_PATH),
    "tempo_module": tempo.__file__,
    "official_repository": "https://github.com/bauerbach95/tempo",
    "tutorials": [
        "tutorial_running_tempo.ipynb",
        "tutorial_tempo_inputs.ipynb",
        "tutorial_tempo_outputs.ipynb",
    ],
    "input_shape": list(adata.shape),
    "cell_phase_prior": "noninformative",
    "seed": SEED,
}

for label, path in [
    ("gene_prior", GENE_PRIOR_PATH),
    ("core_clock", CORE_CLOCK_PATH),
    ("unsupervised_alg", Path(inspect.getfile(unsupervised_alg))),
]:
    provenance[label + "_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()

(OUT_DIR / "input_provenance.json").write_text(
    json.dumps(provenance, indent=2)
)

unsupervised_alg.run(
    adata=adata,
    folder_out=str(OUT_DIR),
    gene_acrophase_prior_path=str(GENE_PRIOR_PATH),
    core_clock_gene_path=str(CORE_CLOCK_PATH),
    reference_gene=reference_gene,
    use_noninformative_phase_prior=True,
    cell_phase_prior_path=None,
)

# ---- notebook code cell 6 ----
# 07 · Tempo run status
opt_posterior_path = OUT_DIR / "tempo_results/opt/cell_posterior.tsv"
step0_posterior_path = OUT_DIR / "tempo_results/0/cell_phase_estimation/cell_posterior.tsv"
clock_evidence_path = OUT_DIR / "evidence/step_0_clock_evidence.txt"
null_evidence_path = OUT_DIR / "evidence/null_log_evidence_vec.txt"

assert step0_posterior_path.is_file(), f"Missing Tempo step-0 posterior: {step0_posterior_path}"
assert clock_evidence_path.is_file(), str(clock_evidence_path)
assert null_evidence_path.is_file(), str(null_evidence_path)

clock_evidence = float(clock_evidence_path.read_text())
null_evidence = np.atleast_1d(np.loadtxt(null_evidence_path))
log10_bf_tempo_vs_null = float(
    (clock_evidence - np.max(null_evidence)) / np.log(10)
)

final_solution_available = opt_posterior_path.is_file()

if final_solution_available:
    posterior_path = opt_posterior_path
    tempo_status = "final_opt_available"
else:
    posterior_path = step0_posterior_path
    tempo_status = "stopped_by_default_clock_evidence_gate"

run_status = {
    "final_solution_available": bool(final_solution_available),
    "tempo_status": tempo_status,
    "log10_bf_tempo_vs_null": log10_bf_tempo_vs_null,
    "clock_evidence": clock_evidence,
    "null_evidence_max": float(np.max(null_evidence)),
    "posterior_source": str(posterior_path),
    "official_threshold_modified": False,
    "evidence_gate_bypassed": False,
}

(OUT_DIR / "Tempo_run_status.json").write_text(
    json.dumps(run_status, indent=2)
)
pd.DataFrame([run_status]).to_csv(
    OUT_DIR / "Tempo_run_status.csv",
    index=False,
)

print(json.dumps(run_status, indent=2))

# ---- notebook code cell 7 ----
# 08 · Official MAP phase
posterior = pd.read_table(
    posterior_path,
    sep="\t",
    index_col="barcode",
)

assert posterior.index.is_unique
assert len(posterior) == 9000
assert set(posterior.index) == set(cell_ids)

posterior = posterior.loc[cell_ids]

assert posterior.index.equals(cell_ids)
assert list(posterior.columns) == [
    f"bin_{i}" for i in range(posterior.shape[1])
]

probabilities = posterior.to_numpy(dtype=float)

assert probabilities.shape[1] >= 2
assert np.isfinite(probabilities).all()
assert (probabilities >= 0).all()

row_sums = probabilities.sum(axis=1)
assert np.allclose(row_sums, 1.0, rtol=1e-5, atol=1e-6)

posterior_obj = cell_posterior.ThetaPosteriorDist(
    torch.tensor(probabilities, dtype=torch.float32)
)

tempo_map_rad = (
    posterior_obj.map_phase.detach().cpu().numpy().astype(float)
)

pred_rad = (
    tempo_map_rad + TEMPO_TO_SIMULATION_ROTATION_RAD
) % (2 * np.pi)

pred_hour = pred_rad * 24 / (2 * np.pi)

assert np.isfinite(pred_rad).all()
assert ((pred_rad >= 0) & (pred_rad < 2 * np.pi)).all()
assert ((pred_hour >= 0) & (pred_hour < 24)).all()

# ---- notebook code cell 8 ----
# 09 · Evaluation
metadata_columns = [
    "cell", "sample_name", "ZT", "celltype",
    "sigma_bio_hour", "T_true", "Theta_true",
]

metadata = pd.read_csv(
    METADATA_PATH,
    dtype={"cell": str},
)

assert set(metadata_columns).issubset(metadata.columns)
assert len(metadata) == 9000
assert metadata["cell"].is_unique
assert set(metadata["cell"]) == set(cell_ids)

metadata = metadata.set_index("cell").loc[cell_ids]

assert metadata.index.equals(cell_ids)
assert not metadata[metadata_columns[1:]].isna().any().any()
assert set(metadata["celltype"]) == {f"CT{i}" for i in range(1, 7)}

assert np.allclose(
    metadata["T_true"],
    source.obs.loc[cell_ids, "true_phase_hour"],
)
assert np.allclose(
    metadata["Theta_true"],
    source.obs.loc[cell_ids, "true_phase_rad"],
)

result = metadata[metadata_columns[1:]].copy()
result.index.name = "cell"
result["pred_hour"] = pred_hour
result["pred_rad"] = pred_rad
result["wrapped_error_hour"] = (
    (result["pred_hour"] - result["T_true"] + 12) % 24
) - 12
result["abs_error_hour"] = result["wrapped_error_hour"].abs()
result["tempo_map_rad"] = tempo_map_rad
result["fixed_rotation_rad"] = TEMPO_TO_SIMULATION_ROTATION_RAD
result["tempo_status"] = tempo_status
result["posterior_source"] = str(posterior_path)
result["log10_bf_tempo_vs_null"] = log10_bf_tempo_vs_null
result = result.reset_index()

assert len(result) == 9000
assert result["cell"].tolist() == cell_ids.tolist()
assert not result.isna().any().any()

# ---- notebook code cell 9 ----
# 10 · Exports
def make_summary(frame):
    rows = []
    for group in ["overall", *[f"CT{i}" for i in range(1, 7)]]:
        x = frame if group == "overall" else frame.loc[frame["celltype"] == group]
        error = x["wrapped_error_hour"].to_numpy()
        rows.append({
            "group": group,
            "n_cells": len(x),
            "circular_MAE_hour": float(np.mean(np.abs(error))),
            "circular_MedAE_hour": float(np.median(np.abs(error))),
            "circular_RMSE_hour": float(np.sqrt(np.mean(error ** 2))),
        })
    return pd.DataFrame(rows)

summary = make_summary(result)
assert summary.loc[1:, "n_cells"].sum() == 9000

formal_csv = OUT_DIR / "Tempo_cell_phase.csv"
formal_summary_csv = OUT_DIR / "Tempo_circular_error_summary.csv"
diagnostic_csv = OUT_DIR / "Tempo_step0_cell_phase_diagnostic.csv"
diagnostic_summary_csv = OUT_DIR / "Tempo_step0_circular_error_summary.csv"

if final_solution_available:
    result.to_csv(formal_csv, index=False)
    summary.to_csv(formal_summary_csv, index=False)

    if diagnostic_csv.exists():
        diagnostic_csv.unlink()
    if diagnostic_summary_csv.exists():
        diagnostic_summary_csv.unlink()

    output_csv = formal_csv
    output_summary_csv = formal_summary_csv
else:
    if formal_csv.exists():
        formal_csv.unlink()
    if formal_summary_csv.exists():
        formal_summary_csv.unlink()

    result.to_csv(diagnostic_csv, index=False)
    summary.to_csv(diagnostic_summary_csv, index=False)

    output_csv = diagnostic_csv
    output_summary_csv = diagnostic_summary_csv

qc = {
    "n_cells": 9000,
    "n_genes": 1000,
    "n_clock_genes": len(core_genes),
    "cell_alignment": True,
    "raw_counts": True,
    "n_output_cells": len(result),
    "posterior_n_bins": probabilities.shape[1],
    "posterior_row_sum_min": float(row_sums.min()),
    "posterior_row_sum_max": float(row_sums.max()),
    "pred_hour_min": float(pred_hour.min()),
    "pred_hour_max": float(pred_hour.max()),
    "n_missing_values": int(result.isna().sum().sum()),
    "tempo_status": tempo_status,
    "final_solution_available": bool(final_solution_available),
    "posterior_source": str(posterior_path),
    "log10_bf_tempo_vs_null": log10_bf_tempo_vs_null,
}

pd.DataFrame([qc]).to_csv(
    OUT_DIR / "Tempo_QC.csv",
    index=False,
)

print(summary.to_string(index=False))

# ---- notebook code cell 10 ----
# 11 · Result paths
print("Tempo status:", tempo_status)
print("Final solution available:", final_solution_available)
print("log10 BF Tempo vs null:", log10_bf_tempo_vs_null)
print("Result CSV:", output_csv.resolve())
print("Summary CSV:", output_summary_csv.resolve())

if final_solution_available:
    print(
        f'tempo <- read.csv("{formal_csv.resolve()}", check.names = FALSE)'
    )
else:
    print("Tempo stopped under its default evidence criterion.")
    print("No formal Tempo_cell_phase.csv was generated.")
    print(
        f'tempo_step0 <- read.csv("{diagnostic_csv.resolve()}", check.names = FALSE)'
    )

# ---- notebook code cell 11 ----
# 12 · Posterior circular mean

phase_grid = np.arange(probabilities.shape[1]) * 2 * np.pi / probabilities.shape[1]

posterior_cos = probabilities @ np.cos(phase_grid)
posterior_sin = probabilities @ np.sin(phase_grid)

tempo_mean_rad = np.arctan2(posterior_sin, posterior_cos) % (2 * np.pi)
tempo_mean_R = np.sqrt(posterior_cos**2 + posterior_sin**2)

pred_rad_mean = (tempo_mean_rad + TEMPO_TO_SIMULATION_ROTATION_RAD) % (2 * np.pi)
pred_hour_mean = pred_rad_mean * 24 / (2 * np.pi)

Tempo_mean = metadata[["sample_name", "ZT", "celltype", "sigma_bio_hour", "T_true", "Theta_true"]].copy()
Tempo_mean.index.name = "cell"
Tempo_mean["pred_hour"] = pred_hour_mean
Tempo_mean["pred_rad"] = pred_rad_mean
Tempo_mean["wrapped_error_hour"] = ((Tempo_mean["pred_hour"] - Tempo_mean["T_true"] + 12) % 24) - 12
Tempo_mean["abs_error_hour"] = Tempo_mean["wrapped_error_hour"].abs()
Tempo_mean["tempo_circular_mean_rad"] = tempo_mean_rad
Tempo_mean["posterior_R"] = tempo_mean_R
Tempo_mean["fixed_rotation_rad"] = TEMPO_TO_SIMULATION_ROTATION_RAD
Tempo_mean = Tempo_mean.reset_index()

Tempo_mean.to_csv(OUT_DIR / "Tempo_cell_phase_circular_mean.csv", index=False)

Tempo_mean.head()

# ---- notebook code cell 12 ----
# 13 · Circular mean summary

summary_mean = []

for group in ["overall", "CT1", "CT2", "CT3", "CT4", "CT5", "CT6"]:
    x = Tempo_mean if group == "overall" else Tempo_mean.loc[Tempo_mean["celltype"] == group]
    e = x["wrapped_error_hour"].to_numpy()
    summary_mean.append({"group": group, "n_cells": len(x), "circular_MAE_hour": np.mean(np.abs(e)), "circular_MedAE_hour": np.median(np.abs(e)), "circular_RMSE_hour": np.sqrt(np.mean(e**2))})

Tempo_mean_summary = pd.DataFrame(summary_mean)
Tempo_mean_summary.to_csv(OUT_DIR / "Tempo_circular_mean_error_summary.csv", index=False)

Tempo_mean_summary
