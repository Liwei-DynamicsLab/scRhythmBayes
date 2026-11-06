# Generated from the named notebook; original file is unchanged.
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# ---- notebook code cell 0 ----
from pathlib import Path
import os, sys, json, random, hashlib, inspect, datetime, platform, gc
os.environ["TF_DETERMINISTIC_OPS"] = "1"
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.io import mmread
from sklearn.preprocessing import scale
import sklearn
import matplotlib.pyplot as plt

INPUT_DIR = Path("outputs/simulation/output/scRITMO_extended_simulation/sigma_gradient_1p5_2p5_3p5_4p5_5p5_6p5_empirical_fixed_v1/R_input")
COUNT_PATH = INPUT_DIR / "count_gene_by_cell.mtx"
GENE_PATH = INPUT_DIR / "genes.csv"
METADATA_PATH = INPUT_DIR / "metadata.csv"
OUT_DIR = Path("07_runtime/runs/Cyclum")
OUT_DIR.mkdir(parents=True, exist_ok=True)
SEED = 1
np.random.seed(SEED)
random.seed(SEED)

def require(condition, message):
    if not condition:
        raise ValueError(message)

def sha256(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def save_json(name, obj):
    (OUT_DIR / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False))

# A rerun invalidates any prior completion marker. R checks this status before importing.
run_info = {"status": "started", "started_utc": datetime.datetime.utcnow().isoformat() + "Z",
            "fit_scope": "per_celltype", "seed": SEED, "python": sys.version, "executable": sys.executable}
save_json("run_manifest.json", run_info)
print("Output:", OUT_DIR)

CT_ORDER = [f"CT{i}" for i in range(1, 7)]


# ---- notebook code cell 1 ----
import tensorflow as tf
import keras
sys.path.insert(0, "external/Cyclum")
import cyclum
from cyclum.tuning import CyclumAutoTune

require(tf.__version__ == "2.15.1", "Expected TensorFlow 2.15.1, got " + tf.__version__)
require(keras.__version__.startswith("2.15."), "Expected Keras 2.15.x")
require(platform.machine() == "arm64", "Select the Cyclum ARM kernel")
tf.config.experimental.enable_op_determinism()
keras.backend.clear_session()
tf.random.set_seed(SEED)
np.random.seed(SEED)
random.seed(SEED)
TUNE_KWARGS = dict(max_linear_dims=5, epochs=500, verbose=100, rate=5e-4, early_stop=False)
inspect.signature(CyclumAutoTune).bind(np.zeros((5, 5)), **TUNE_KWARGS)
TRAIN_KWARGS = dict(batch_size=None, epochs=100, verbose=10, rate=1e-4)
inspect.signature(CyclumAutoTune.train).bind(None, np.zeros((5, 5)), **TRAIN_KWARGS)
inspect.signature(CyclumAutoTune.predict_pseudotime).bind(None, np.zeros((5, 5)))
api_names = ["__init__", "train", "predict_pseudotime", "show_elbow", "show_bar", "save"]
api = {name: str(inspect.signature(getattr(CyclumAutoTune, name))) for name in api_names}
source_text = "\n\n".join(inspect.getsource(getattr(CyclumAutoTune, name)) for name in api_names)
(OUT_DIR / "Cyclum_installed_API_source.txt").write_text(source_text)
run_info.update({"versions": {"tensorflow": tf.__version__, "keras": keras.__version__, "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "sklearn": sklearn.__version__}, "cyclum_path": cyclum.__file__, "api_signatures": api, "autotune_kwargs": TUNE_KWARGS, "final_train_kwargs": TRAIN_KWARGS, "final_train_signature": api["train"], "api_source_sha256": sha256(OUT_DIR / "Cyclum_installed_API_source.txt")})
save_json("run_manifest.json", run_info)
print("TensorFlow:", tf.__version__)
print("Keras:", keras.__version__)
print("Cyclum:", cyclum.__file__)
print(json.dumps(api, indent=2))

run_info["determinism"] = "Fixed per-CT seeds and TensorFlow deterministic ops; bitwise equality across hardware/library versions is not guaranteed."
run_info["ct_seeds"] = {ct: SEED + i for i, ct in enumerate(CT_ORDER)}
run_info["cyclum_source_sha256"] = {str(p.relative_to(Path(cyclum.__file__).parent)): sha256(p) for p in sorted(Path(cyclum.__file__).parent.rglob("*.py"))}
save_json("run_manifest.json", run_info)


# ---- notebook code cell 2 ----
count_gene_by_cell = mmread(str(COUNT_PATH)).tocsr()
genes = pd.read_csv(GENE_PATH, dtype={"gene": str})
require("gene" in genes, "genes.csv must contain gene")
require(len(genes) == count_gene_by_cell.shape[0], "Gene dimension mismatch")
require(genes["gene"].notna().all() and not genes["gene"].duplicated().any(), "Missing/duplicate genes")
require(np.isfinite(count_gene_by_cell.data).all(), "Non-finite counts")
require((count_gene_by_cell.data >= 0).all(), "Negative counts")
require((count_gene_by_cell.data == np.floor(count_gene_by_cell.data)).all(), "Noninteger counts")
X_counts = count_gene_by_cell.T.toarray().astype(np.float64)
n_cells, n_genes = X_counts.shape
require((n_cells, n_genes) == (9000, 1000), "This is the fixed 6CT dataset: expected 9000 cells × 1000 genes")
row_index = np.arange(n_cells)
run_info.update({"shape_cells_genes": [n_cells, n_genes],
    "input_sha256": {"count_gene_by_cell.mtx": sha256(COUNT_PATH), "genes.csv": sha256(GENE_PATH)},
    "row_order": "Matrix column order; checked against original h5ad before per-CT inference"})
save_json("run_manifest.json", run_info)
print("Cyclum input:", X_counts.shape)


# ---- notebook code cell 3 ----
DESIGN_COLS = ["cell", "sample_name", "ZT", "celltype"]
metadata = pd.read_csv(METADATA_PATH, usecols=DESIGN_COLS, dtype={"cell": str, "sample_name": str, "celltype": str})
require(set(DESIGN_COLS).issubset(metadata.columns), "Missing metadata columns")
metadata = metadata[DESIGN_COLS].copy()
require(len(metadata) == n_cells and not metadata.isna().any().any(), "Metadata size/missing-value error")
require(not metadata["cell"].duplicated().any(), "Duplicate cell IDs")
require(set(metadata["celltype"]) == {"CT" + str(i) for i in range(1, 7)}, "Expected CT1–CT6")
require((metadata.groupby("sample_name")["ZT"].nunique() == 1).all(), "A sample has inconsistent ZT")
# Matrix Market has no IDs. Verify the export against its original count matrix and IDs.
import h5py
H5AD_PATH = INPUT_DIR.parent / "scRITMO_extended_simulation.h5ad"
def decode_index(group):
    key = group.attrs.get("_index", "_index")
    if isinstance(key, bytes): key = key.decode()
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in group[key][:]])
with h5py.File(str(H5AD_PATH), "r") as h:
    source_cells = decode_index(h["obs"])
    source_genes = decode_index(h["var"])
    require(np.array_equal(source_cells, metadata["cell"].values), "Metadata order differs from original h5ad")
    require(np.array_equal(source_genes, genes["gene"].values), "Gene order differs from original h5ad")
    matrix_key = "layers/counts" if "layers/counts" in h else "X"
    x = h[matrix_key]
    if isinstance(x, h5py.Dataset):
        original_counts = sparse.csr_matrix(x[:])
    else:
        encoding = x.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes): encoding = encoding.decode()
        require(encoding in ["csr_matrix", "csc_matrix"], "Unsupported h5ad sparse encoding")
        constructor = sparse.csr_matrix if encoding == "csr_matrix" else sparse.csc_matrix
        original_counts = constructor((x["data"][:], x["indices"][:], x["indptr"][:]), shape=tuple(x.attrs["shape"]))
    require(original_counts.shape == (n_cells, n_genes), "Original matrix shape mismatch")
    require((original_counts.tocsr() != count_gene_by_cell.T.tocsr()).nnz == 0,
            "Counts do not match original h5ad; cell order cannot be verified")

require(np.isfinite(metadata["ZT"].to_numpy(dtype=float)).all(), "Invalid ZT")
require((metadata.groupby("celltype").size().reindex(CT_ORDER) >= TUNE_KWARGS["max_linear_dims"] + 2).all(), "Too few cells for PCA within a CT")
run_info["input_sha256"]["metadata.csv"] = sha256(METADATA_PATH)
run_info["input_sha256"]["source_h5ad"] = sha256(H5AD_PATH)
run_info["cell_order_verified_against"] = str(H5AD_PATH) + ":" + matrix_key
save_json("run_manifest.json", run_info)
del original_counts


# ---- notebook code cell 4 ----
def align_to_sample_zt(p, sample_zt_rad, n_s=200):
    p = np.asarray(p, dtype=float)
    p0 = np.asarray(sample_zt_rad, dtype=float)
    require(p.ndim == 1 and p.shape == p0.shape and len(p) > 0, "Invalid alignment shapes")
    require(np.isfinite(p).all() and np.isfinite(p0).all(), "Non-finite alignment inputs")
    def circular_deviation2(x, y, period=2 * np.pi):
        x, y = x % period, y % period
        v1 = np.abs(x - y)
        v2 = (period - v1) % period
        return np.minimum(v1, v2)
    shifts = np.linspace(0, 2 * np.pi, n_s)
    theta_cs = (p[:, None] - shifts[None, :]) % (2 * np.pi)
    theta_cs_neg = (-p[:, None] - shifts[None, :]) % (2 * np.pi)
    v = np.median(circular_deviation2(p0[:, None], theta_cs), axis=0)
    v_neg = np.median(circular_deviation2(p0[:, None], theta_cs_neg), axis=0)
    i, j = np.argmin(v), np.argmin(v_neg)
    if v[i] < v_neg[j]:
        aligned, mad, orientation, shift = theta_cs[:, i], v[i], 1, shifts[i]
    else:
        aligned, mad, orientation, shift = theta_cs_neg[:, j], v_neg[j], -1, shifts[j]
    candidates = pd.DataFrame({"shift_rad": np.tile(shifts, 2),
        "orientation": np.repeat([1, -1], n_s),
        "ZT_circular_MedAE_rad": np.concatenate([v, v_neg])})
    return aligned, float(mad), orientation, float(shift), candidates


# ---- notebook code cell 5 ----
raw_parts, aligned_parts, alignment_rows, ct_qc_rows = [], [], [], []
run_info["ct_status"] = {}
for ct in CT_ORDER:
    ct_dir = OUT_DIR / ct
    ct_dir.mkdir(parents=True, exist_ok=True)
    run_info["ct_status"][ct] = "started"
    save_json("run_manifest.json", run_info)
    try:
        keras.backend.clear_session()
        gc.collect()
        ct_seed = run_info["ct_seeds"][ct]
        tf.keras.utils.set_random_seed(ct_seed)
        idx = np.flatnonzero(metadata["celltype"].to_numpy() == ct)
        X_log = np.log2(X_counts[idx, :] + 1.0)
        log_mean, log_sd = X_log.mean(axis=0), X_log.std(axis=0, ddof=0)
        X_ct = np.asarray(scale(X_log, axis=0), dtype=np.float32)
        variable = log_sd > 0
        require(variable.any() and np.isfinite(X_ct).all(), ct + ": invalid scaled expression")
        require(np.max(np.abs(X_ct.mean(axis=0, dtype=np.float64))) < 1e-5, "Centering failed")
        require(np.allclose(X_ct[:, variable].std(axis=0, dtype=np.float64), 1, atol=1e-5), "Scaling failed")
        pd.DataFrame({"gene": genes["gene"], "log2_mean": log_mean, "log2_sd": log_sd,
                      "zero_variance": ~variable, "used": True}).to_csv(ct_dir / "Cyclum_gene_preprocessing_QC.csv", index=False)
        np.save(ct_dir / "Cyclum_input_scaled.npy", X_ct)
        np.save(ct_dir / "Cyclum_count_column_index.npy", idx)
        model = CyclumAutoTune(X_ct, **TUNE_KWARGS)
        require(np.isfinite(model.pca_loss).all() and np.isfinite(model.ae_loss).all(), "Invalid AutoTune losses")
        bypass = [layer for layer in model.model.layers if layer.name == "bypass_out"]
        selected = int(bypass[0].units) if bypass else 0
        boundary = selected == TUNE_KWARGS["max_linear_dims"]
        losses = ([{"model": "PCA", "dimension": i+1, "loss": float(v)} for i,v in enumerate(model.pca_loss)] +
                  [{"model": "circular_AE", "dimension": i, "loss": float(v)} for i,v in enumerate(model.ae_loss)])
        pd.DataFrame(losses).assign(celltype=ct).to_csv(ct_dir / "Cyclum_autotune_losses.csv", index=False)
        for name, fig in [("elbow", model.show_elbow()), ("bar", model.show_bar())]:
            fig.savefig(ct_dir / ("Cyclum_autotune_" + name + ".png"), dpi=180, bbox_inches="tight")
            plt.close(fig)
        history = model.train(X_ct, **TRAIN_KWARGS)
        require(np.isfinite(history.history["loss"]).all(), "Invalid training loss")
        pd.DataFrame(history.history).rename_axis("epoch_zero_based").to_csv(ct_dir / "Cyclum_training_history.csv")
        pseudotime = np.asarray(model.predict_pseudotime(X_ct)).reshape(-1)
        require(pseudotime.shape == (len(idx),) and np.isfinite(pseudotime).all(), "Invalid pseudotime")
        raw_rad = np.mod(pseudotime.astype(np.float64), 2*np.pi)
        raw_ct = pd.DataFrame({"cell": metadata.iloc[idx]["cell"].to_numpy(), "count_column_index": idx,
            "cyclum_pseudotime_raw": pseudotime, "cyclum_theta_rad_raw": raw_rad,
            "cyclum_theta_hour_raw": raw_rad * 24/(2*np.pi), "ct_seed": ct_seed,
            "selected_linear_dims": selected, "selected_at_search_boundary": boundary})
        raw_ct.to_csv(ct_dir / "Cyclum_cell_phase_raw.csv", index=False)
        model.save(str(ct_dir / "Cyclum_final_model.h5"))
        np.save(ct_dir / "Cyclum_circular_loadings.npy", model.get_circular_loadings())
        zt_rad = (metadata.iloc[idx]["ZT"].to_numpy(dtype=float) % 24) * (2*np.pi/24)
        pred_rad, mad, orientation, shift, candidates = align_to_sample_zt(raw_rad, zt_rad)
        candidates.to_csv(ct_dir / "Cyclum_alignment_candidates.csv", index=False)
        alignment = {"celltype": ct, "n_cells": len(idx), "reference": "within-CT sample ZT only",
            "objective": "cell-weighted circular median absolute deviation", "n_shifts": 200,
            "grid_includes_endpoint": True, "orientation": orientation, "shift_rad": shift,
            "formula": "pred_rad = (orientation * raw_rad - shift_rad) % (2*pi)",
            "ZT_alignment_MAD_hour": mad * 24/(2*np.pi)}
        pd.DataFrame([alignment]).to_csv(ct_dir / "Cyclum_alignment_summary.csv", index=False)
        save_json(ct + "/Cyclum_alignment.json", alignment)
        aligned_ct = pd.DataFrame({"cell": raw_ct["cell"], "pred_hour": pred_rad*24/(2*np.pi),
            "pred_rad": pred_rad, "alignment_orientation": orientation, "alignment_shift_rad": shift})
        aligned_ct.to_csv(ct_dir / "Cyclum_aligned_predictions.csv", index=False)
        raw_parts.append(raw_ct); aligned_parts.append(aligned_ct); alignment_rows.append(alignment)
        qc = {"celltype": ct, "n_cells": len(idx), "seed": ct_seed, "n_genes": n_genes,
            "zero_variance_genes": int((~variable).sum()), "selected_linear_dims": selected,
            "selected_at_search_boundary": boundary, "autotune_internal_solution_found": not boundary,
            "final_loss": float(history.history["loss"][-1]), "train_epochs": len(history.history["loss"]),
            "raw_circular_resultant_length": float(abs(np.mean(np.exp(1j*raw_rad))))}
        ct_qc_rows.append(qc)
        save_json(ct + "/Cyclum_QC.json", qc)
        run_info["ct_status"][ct] = "complete"
        save_json("run_manifest.json", run_info)
        print(ct, qc)
        if boundary: print("QC WARNING:", ct, "AutoTune reached max_linear_dims=5; no internal solution found.")
        del model, history, X_ct, X_log
    except Exception as exc:
        run_info.update(status="failed", failed_celltype=ct, error=repr(exc))
        run_info["ct_status"][ct] = "failed"
        save_json("run_manifest.json", run_info)
        raise
raw = pd.concat(raw_parts, ignore_index=True).sort_values("count_column_index").reset_index(drop=True)
aligned = pd.concat(aligned_parts, ignore_index=True).set_index("cell").loc[metadata["cell"]].reset_index()
require(len(raw) == 9000 and np.array_equal(raw["cell"], metadata["cell"]), "Incomplete or misordered predictions")
require(not aligned["cell"].duplicated().any() and len(aligned) == 9000, "Invalid prediction coverage")
raw.to_csv(OUT_DIR / "Cyclum_cell_phase_raw.csv", index=False)
aligned.to_csv(OUT_DIR / "Cyclum_aligned_predictions.csv", index=False)
pd.DataFrame(alignment_rows).to_csv(OUT_DIR / "Cyclum_alignment_summary.csv", index=False)
pd.DataFrame(ct_qc_rows).to_csv(OUT_DIR / "Cyclum_QC_by_celltype.csv", index=False)
run_info.update(status="predictions_frozen", alignment=alignment_rows, celltype_qc=ct_qc_rows,
                frozen_prediction_sha256=sha256(OUT_DIR / "Cyclum_aligned_predictions.csv"))
save_json("run_manifest.json", run_info)


# ---- notebook code cell 6 ----
require(run_info["status"] == "predictions_frozen", "Freeze all CT predictions first")
for name, path in [("metadata.csv", METADATA_PATH), ("genes.csv", GENE_PATH), ("count_gene_by_cell.mtx", COUNT_PATH), ("source_h5ad", H5AD_PATH)]:
    require(sha256(path) == run_info["input_sha256"][name], "Input changed during run: " + name)
META_COLS = ["cell", "sample_name", "ZT", "celltype", "sigma_bio_hour", "T_true", "Theta_true"]
evaluation_metadata = pd.read_csv(METADATA_PATH, usecols=META_COLS, dtype={"cell": str, "sample_name": str, "celltype": str})
require(evaluation_metadata[DESIGN_COLS].equals(metadata[DESIGN_COLS]), "Design metadata changed")
require(not evaluation_metadata[META_COLS].isna().any().any(), "Missing evaluation metadata")
require(np.isfinite(evaluation_metadata[["ZT", "sigma_bio_hour", "T_true", "Theta_true"]].to_numpy(dtype=float)).all(), "Nonfinite evaluation metadata")
metadata = evaluation_metadata[META_COLS].copy()
def wrap_hour(x):
    return (np.asarray(x) + 12) % 24 - 12

Cyclum_phase = metadata.merge(raw, on="cell", how="left", validate="one_to_one").merge(
    aligned, on="cell", how="left", validate="one_to_one")
require(np.isfinite(Cyclum_phase[["T_true", "pred_hour", "pred_rad"]].to_numpy(dtype=float)).all(), "Invalid evaluation values")
Cyclum_phase["wrapped_error_hour"] = wrap_hour(Cyclum_phase["pred_hour"] - Cyclum_phase["T_true"])
Cyclum_phase["abs_error_hour"] = np.abs(Cyclum_phase["wrapped_error_hour"])
require(((Cyclum_phase["pred_hour"] >= 0) & (Cyclum_phase["pred_hour"] < 24)).all(), "Prediction outside [0,24)")
rows = []
for group in ["overall"] + ["CT" + str(i) for i in range(1, 7)]:
    subset = Cyclum_phase if group == "overall" else Cyclum_phase[Cyclum_phase["celltype"] == group]
    e = subset["wrapped_error_hour"].to_numpy()
    require(len(e) > 0, "Empty evaluation group " + group)
    rows.append({"group": group, "n_cells": len(e), "circular_MAE_hour": float(np.mean(np.abs(e))),
        "circular_MedAE_hour": float(np.median(np.abs(e))), "circular_RMSE_hour": float(np.sqrt(np.mean(e ** 2)))})
Cyclum_summary = pd.DataFrame(rows)
STANDARD_COLS = META_COLS + ["pred_hour", "pred_rad", "wrapped_error_hour", "abs_error_hour"]
Cyclum_phase = Cyclum_phase[STANDARD_COLS + [c for c in Cyclum_phase if c not in STANDARD_COLS]]
Cyclum_phase.to_csv(OUT_DIR / "Cyclum_cell_phase.csv", index=False)
Cyclum_summary.to_csv(OUT_DIR / "Cyclum_circular_error_summary.csv", index=False)
run_info.update({"status": "complete", "finished_utc": datetime.datetime.utcnow().isoformat() + "Z",
    "alignment": alignment_rows, "any_boundary_selection": any(q["selected_at_search_boundary"] for q in ct_qc_rows), "prediction_sha256": sha256(OUT_DIR / "Cyclum_cell_phase.csv")})
save_json("run_manifest.json", run_info)
print(Cyclum_summary.to_string(index=False))
print("R benchmark input:", OUT_DIR / "Cyclum_cell_phase.csv")

