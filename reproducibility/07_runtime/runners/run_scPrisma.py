# Generated from the benchmark notebook; original file is unchanged.
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# ---- notebook code cell 0 ----
from pathlib import Path
import os, sys, json, random, hashlib, inspect, datetime, platform, gc, time, contextlib
INPUT_DIR = Path("/path/to/analysis_home/Projects/scRhythmBayes_simulation/output/scRITMO_extended_simulation/sigma_gradient_1p5_2p5_3p5_4p5_5p5_6p5_empirical_fixed_v1/R_input")
OUT_DIR = Path("/path/to/scRhythmBayes/reproducibility/07_runtime/runs/scPrisma")
OUT_DIR.mkdir(parents=True, exist_ok=True)
COUNT_PATH, GENE_PATH, METADATA_PATH = [INPUT_DIR / n for n in ("count_gene_by_cell.mtx", "genes.csv", "metadata.csv")]
SEED = 1
CT_ORDER = [f"CT{i}" for i in range(1, 7)]
ZT_ORDER = [0, 4, 8, 12, 16, 20]
N_TOP_GENES = 7000
REORDER_KWARGS = dict(iterNum=300, batch_size=20, gamma=0, lr=0.1, optimize_alpha=True)

def require(condition, message):
    if not condition: raise ValueError(message)
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""): h.update(block)
    return h.hexdigest()
def save_json(name, obj):
    p = OUT_DIR / name
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False))
    tmp.replace(p)
run_info = dict(status="started", started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    method="scPrisma", variant="per_CT_hard_sample_ZT_prior", fit_scope="per_celltype", seed=SEED,
    python=sys.version, executable=sys.executable, input_dir=str(INPUT_DIR), output_dir=str(OUT_DIR),
    reorder_kwargs=REORDER_KWARGS, n_top_genes=N_TOP_GENES, ct_status={},
    preprocessing="per CT: min_cells=3; total=1e4; log1p; remove constants; <=7000 HVG; no z-score; internal row L2",
    truth_policy="T_true/Theta_true/sigma_bio_hour loaded only after all predictions frozen",
    calibration="existing Cyclum 200-shift two-orientation cell-weighted median circular deviation to sample ZT",
    rank_conversion="raw_rad[cell] = 2*pi*inverse_permutation[cell]/n_CT",
    determinism="Fixed numpy seed per CT and shuffled starting cell order; cross-platform bitwise identity is not guaranteed")
save_json("run_manifest.json", run_info)
print("Output:", OUT_DIR)


# ---- notebook code cell 1 ----
try:
    import numpy as np
    import pandas as pd
    import scipy
    from scipy import sparse
    from scipy.io import mmread
    import scanpy as sc
    import anndata as ad
    import scPrisma
    import matplotlib.pyplot as plt
    import importlib.metadata as im
    from threadpoolctl import threadpool_limits
    require(sys.version_info[:2] == (3, 10), "Select the scPrisma Python 3.10 kernel")
    np.random.seed(SEED); random.seed(SEED)
    inspect.signature(scPrisma.algorithms.reorder_indicator).bind(np.zeros((6, 10)), np.eye(6), **REORDER_KWARGS)
    inspect.signature(scPrisma.algorithms.e_to_range).bind(np.eye(6))
    api_names = ["reorder_indicator", "e_to_range", "reconstruct_e", "sga_matrix_momentum_indicator", "cell_normalization", "ge_to_spectral_matrix"]
    api_source = "\n\n".join(inspect.getsource(getattr(scPrisma.algorithms, k)) for k in api_names)
    (OUT_DIR / "scPrisma_installed_API_source.txt").write_text(api_source)
    run_info.update(scprisma_path=scPrisma.__file__,
        api_signatures={k:str(inspect.signature(getattr(scPrisma.algorithms,k))) for k in api_names},
        api_source_sha256=sha256(OUT_DIR / "scPrisma_installed_API_source.txt"),
        versions={k:im.version(k) for k in ["scPrisma", "numpy", "pandas", "scipy", "scanpy", "anndata", "numba", "scikit-learn", "threadpoolctl"]},
        source_sha256={str(p.relative_to(Path(scPrisma.__file__).parent)):sha256(p) for p in sorted(Path(scPrisma.__file__).parent.rglob("*.py"))},
        ct_seeds={ct: SEED+i for i,ct in enumerate(CT_ORDER)})
    print("scPrisma:", scPrisma.__file__)
    print("API:", run_info["api_signatures"]["reorder_indicator"])

    run_info["tutorial_sources_at_authoring"] = [{'path': '/path/to/analysis_home/Projects/scPrisma/tests/cpu/test_cpu.py', 'sha256': '9e78ffeef7e066454a1a80cb1031c8c53fb68537d6becdf2a5af020b220c47bb'}, {'path': '/path/to/analysis_home/Projects/scPrisma/tutorials/cpu/hematopoietic_progenitors_reconstruction_with_partial_prior_knowledge.ipynb', 'sha256': '53644323903bc46bc27c02bdcb8fba237fb615e762a510ba6d452ef1dc65a52f'}, {'path': '/path/to/analysis_home/Projects/scPrisma/tutorials/cpu/tutorial_de_novo_reconstruction.ipynb', 'sha256': '5b145efa9ba7382b078cb85a6d3b7aaa7e416bdf06b8c0976b3a2c22015c5ed1'}, {'path': '/path/to/analysis_home/Projects/scPrisma/tutorials/cpu/tutorial_prior_knowledge_linear_and_cyclic_cpu.ipynb', 'sha256': '5a99eb822494a3fc1f3379dd2df964510e4037adfdbd5644c91a952d34b62c48'}]
    run_info["tutorial_upload_identity_verified"] = False
    save_json("run_manifest.json", run_info)
except BaseException as exc:
    run_info.update(status="failed", error=repr(exc))
    save_json("run_manifest.json", run_info)
    raise


# ---- notebook code cell 2 ----
try:
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
    print("scPrisma input:", X_counts.shape)
except BaseException as exc:
    run_info.update(status="failed", error=repr(exc))
    save_json("run_manifest.json", run_info)
    raise


# ---- notebook code cell 3 ----
try:
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
    require((metadata.groupby("celltype").size().reindex(CT_ORDER) >= 12).all(), "Too few cells within a CT")
    for ct in CT_ORDER:
        require(set(metadata.loc[metadata.celltype == ct, "ZT"]) == set(ZT_ORDER), ct + ": unexpected ZT labels")
    run_info["input_sha256"]["metadata.csv"] = sha256(METADATA_PATH)
    run_info["input_sha256"]["source_h5ad"] = sha256(H5AD_PATH)
    run_info["cell_order_verified_against"] = str(H5AD_PATH) + ":" + matrix_key
    save_json("run_manifest.json", run_info)
    del original_counts
except BaseException as exc:
    run_info.update(status="failed", error=repr(exc))
    save_json("run_manifest.json", run_info)
    raise


# ---- notebook code cell 4 ----
try:
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
except BaseException as exc:
    run_info.update(status="failed", error=repr(exc))
    save_json("run_manifest.json", run_info)
    raise


# ---- notebook code cell 5 ----
try:
    def prepare_ct(counts, design, gene_names, seed):
        require(set(design.columns) == set(DESIGN_COLS), "Inference must receive only design metadata")
        n = len(design)
        idx = np.random.default_rng(seed).permutation(n)
        idx = idx[np.argsort(design.iloc[idx]["ZT"].to_numpy(), kind="stable")]
        d = design.iloc[idx].reset_index(drop=True)
        C = counts[idx, :].copy()
        lib = C.sum(axis=1)
        detected_cell = (C > 0).sum(axis=1)
        detected_gene = (C > 0).sum(axis=0)
        require(np.all(lib > 0), "Zero-library cell: cannot preserve all cells with valid inference")
        keep = detected_gene >= 3
        require(keep.sum() >= 2, "Too few detected genes")
        a = ad.AnnData(C[:, keep].copy())
        a.var_names = np.asarray(gene_names)[keep]
        require(np.all(np.asarray(a.X.sum(axis=1)).ravel() > 0), "A cell has zero counts after gene filtering")
        sc.pp.normalize_total(a, target_sum=1e4)
        require(np.allclose(np.asarray(a.X.sum(axis=1)).ravel(), 1e4), "Library normalization failed")
        sc.pp.log1p(a)
        variable = np.var(a.X, axis=0) > 1e-12
        a = a[:, variable].copy()
        require(a.n_vars >= 2, "Too few nonconstant genes")
        if a.n_vars > N_TOP_GENES:
            sc.pp.highly_variable_genes(a, n_top_genes=N_TOP_GENES, flavor="seurat")
            a = a[:, a.var.highly_variable].copy()
        A = np.ascontiguousarray(a.X, dtype=np.float64)
        require(np.isfinite(A).all() and (A >= 0).all() and np.all(np.linalg.norm(A, axis=1) > 0), "Invalid scPrisma input")
        zt = d.ZT.to_numpy(dtype=float)
        indicator = (zt[:, None] == zt[None, :]).astype(np.float64)
        gqc = pd.DataFrame({"gene":gene_names, "detected_cells":detected_gene,
            "passed_min_cells":keep, "selected":np.isin(gene_names, a.var_names)})
        cqc = d.copy()
        cqc["library_size"] = lib
        cqc["detected_genes"] = detected_cell
        cqc["below_tutorial_200_genes"] = detected_cell < 200
        cqc["below_liver_500_genes"] = detected_cell < 500
        cqc["retained"] = True
        return A, indicator, d, idx, gqc, cqc

    def fit_one_ct(counts, design, gene_names, seed, ct_dir, kwargs):
        start = time.monotonic()
        A, indicator, d, idx, gqc, cqc = prepare_ct(counts, design, gene_names, seed)
        n = len(d)
        gqc.to_csv(ct_dir / "gene_preprocessing_QC.csv", index=False)
        cqc.to_csv(ct_dir / "cell_QC.csv", index=False)
        np.save(ct_dir / "input_log_normalized.npy", A)
        np.save(ct_dir / "indicator.npy", indicator.astype(np.uint8))
        np.random.seed(seed)
        # Limit numerical threads to avoid oversubscription on CPU. No API monkey-patching.
        with threadpool_limits(limits=4), open(ct_dir / "optimization.log", "w") as log:
            with contextlib.redirect_stdout(log):
                E, E_rec = scPrisma.algorithms.reorder_indicator(A, indicator, **kwargs)
        require(E.shape == E_rec.shape == (n,n), "Unexpected return shapes")
        require(np.isfinite(E).all() and np.isfinite(E_rec).all(), "Nonfinite reconstruction")
        np.savez_compressed(ct_dir / "reconstruction_matrices.npz", E=E, E_rec=E_rec.astype(np.uint8))
        require(np.isin(E_rec, [0,1]).all(), "E_rec is not binary")
        require(np.all(E_rec.sum(axis=0) == 1) and np.all(E_rec.sum(axis=1) == 1), "E_rec is not a permutation")
        require(np.max(np.abs(E[indicator == 0])) <= 1e-8, "Soft E violates indicator")
        require(np.sum(E_rec[indicator == 0]) == 0, "Rounded permutation crossed ZT blocks; stop, do not silently repair")
        order = np.asarray(scPrisma.algorithms.e_to_range(E_rec))
        require(order.shape == (n,) and np.array_equal(np.sort(order), np.arange(n)), "Invalid e_to_range output")
        require(np.array_equal(E_rec.argmax(axis=1), order), "Permutation direction mismatch")
        require(np.array_equal(d.ZT.to_numpy()[order], d.ZT.to_numpy()), "ZT block preservation failed")
        rank = np.empty(n, dtype=int)
        rank[order] = np.arange(n)
        raw_rad = 2*np.pi*rank/n
        raw = d.copy()
        raw["ct_input_index"] = np.arange(n)
        raw["ct_original_index"] = idx
        raw["rank_zero_based"] = rank
        raw["raw_rad"] = raw_rad
        raw["raw_hour"] = raw_rad*24/(2*np.pi)
        raw["seed"] = seed
        raw.to_csv(ct_dir / "raw_cell_ordering.csv", index=False)
        raw.iloc[order].assign(output_position=np.arange(n)).to_csv(ct_dir / "raw_ordering_by_position.csv", index=False)
        aligned, mad, orientation, shift, candidates = align_to_sample_zt(raw_rad, d.ZT.to_numpy()*2*np.pi/24)
        candidates.to_csv(ct_dir / "alignment_candidates.csv", index=False)
        alignment = dict(celltype=str(d.celltype.iloc[0]), n_cells=n, orientation=orientation, shift_rad=shift,
            ZT_alignment_MAD_hour=mad*24/(2*np.pi), reference="sample ZT only", n_shifts=200,
            objective="cell-weighted median absolute circular deviation", formula="(orientation*raw_rad-shift_rad) mod 2pi")
        pd.DataFrame([alignment]).to_csv(ct_dir / "alignment_summary.csv", index=False)
        pred = d.copy()
        pred["pred_rad"] = aligned
        pred["pred_hour"] = aligned*24/(2*np.pi)
        pred.to_csv(ct_dir / "aligned_predictions.csv", index=False)
        qc = dict(celltype=str(d.celltype.iloc[0]), n_cells=n, n_genes_input=int(counts.shape[1]), n_genes_selected=A.shape[1],
            zero_library_cells=int((cqc.library_size==0).sum()), low_complexity_200=int(cqc.below_tutorial_200_genes.sum()),
            E_row_sum_max_abs_deviation=float(np.max(np.abs(E.sum(axis=1)-1))),
            E_col_sum_max_abs_deviation=float(np.max(np.abs(E.sum(axis=0)-1))), E_min=float(E.min()),
            E_max=float(E.max()), prior_violations=0, permutation_valid=True, seed=seed,
            elapsed_seconds=time.monotonic()-start, iterNum=kwargs["iterNum"],
            convergence_assessed=False, identity_fraction=float(np.mean(order==np.arange(n))))
        blocks = d.groupby("ZT", sort=True).size().rename("n_cells").reset_index()
        blocks["rank_span_hour"] = blocks.n_cells/n*24
        blocks.to_csv(ct_dir / "ZT_blocks.csv", index=False)
        fig, ax = plt.subplots(1,3, figsize=(12,3.5))
        ax[0].imshow(indicator, interpolation="nearest", aspect="auto", cmap="Greys")
        ax[0].set(xlabel="Input cell", ylabel="Output position", title="Hard ZT indicator")
        ax[1].scatter(np.arange(n), order, c=d.ZT.to_numpy()[order], s=2, cmap="twilight", vmin=0, vmax=24)
        ax[1].set(xlabel="Output position", ylabel="Input cell", title="Rounded permutation")
        ax[2].scatter(d.ZT, pred.pred_hour, s=3, alpha=.3)
        ax[2].set(xlabel="Sample ZT (h)", ylabel="Predicted phase (h)", ylim=(0,24), title="ZT calibration (no truth)")
        fig.tight_layout(); fig.savefig(ct_dir / "ordering_QC.png", dpi=160); plt.close(fig)
        return raw, pred, alignment, qc
except BaseException as exc:
    run_info.update(status="failed", error=repr(exc))
    save_json("run_manifest.json", run_info)
    raise


# ---- notebook code cell 6 ----
try:
    raw_parts, pred_parts, alignment_rows, qc_rows = [], [], [], []
    for ct in CT_ORDER:
        ct_dir = OUT_DIR / ct
        ct_dir.mkdir(parents=True, exist_ok=True)
        run_info["ct_status"][ct] = "started"
        save_json("run_manifest.json", run_info)
        idx = np.flatnonzero(metadata.celltype.to_numpy() == ct)
        print(f"{ct}: {len(idx)} cells; starting scPrisma", flush=True)
        raw_ct, pred_ct, alignment, qc = fit_one_ct(X_counts[idx,:], metadata.iloc[idx][DESIGN_COLS].copy(),
            genes.gene.to_numpy(), run_info["ct_seeds"][ct], ct_dir, REORDER_KWARGS)
        raw_ct["count_column_index"] = idx[raw_ct.ct_original_index.to_numpy()]
        raw_ct.to_csv(ct_dir / "raw_cell_ordering.csv", index=False)
        raw_parts.append(raw_ct); pred_parts.append(pred_ct); alignment_rows.append(alignment); qc_rows.append(qc)
        save_json(ct + "/QC.json", qc)
        run_info["ct_status"][ct] = "complete"
        save_json("run_manifest.json", run_info)
        print(ct, qc, flush=True)
        gc.collect()
    raw = pd.concat(raw_parts, ignore_index=True).sort_values("count_column_index").reset_index(drop=True)
    aligned = pd.concat(pred_parts, ignore_index=True).set_index("cell").loc[metadata.cell].reset_index()
    require(len(aligned) == 9000 and not aligned.cell.duplicated().any(), "Incomplete prediction coverage")
    require(np.array_equal(aligned.cell, metadata.cell) and np.array_equal(raw.cell, metadata.cell), "Cell order mismatch")
    require(np.isfinite(aligned[["pred_hour","pred_rad"]]).all().all(), "Invalid predictions")
    require(aligned.pred_hour.between(0,24,inclusive="left").all(), "Predictions outside 0-24 h")
    raw.to_csv(OUT_DIR / "scPrisma_raw_ordering.csv", index=False)
    aligned.to_csv(OUT_DIR / "scPrisma_aligned_predictions.csv", index=False)
    pd.DataFrame(alignment_rows).to_csv(OUT_DIR / "scPrisma_alignment_summary.csv", index=False)
    pd.DataFrame(qc_rows).to_csv(OUT_DIR / "scPrisma_QC_by_celltype.csv", index=False)
    run_info.update(status="predictions_frozen", frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        frozen_prediction_sha256=sha256(OUT_DIR / "scPrisma_aligned_predictions.csv"),
        alignment=alignment_rows, celltype_qc=qc_rows)
    save_json("run_manifest.json", run_info)
    print("All 9000 predictions frozen. Evaluation metadata can now be read.")
except BaseException as exc:
    run_info.update(status="failed", error=repr(exc))
    save_json("run_manifest.json", run_info)
    raise


# ---- notebook code cell 7 ----
try:
    require(run_info["status"] == "predictions_frozen", "Freeze all CT predictions first")
    for name, path in [("metadata.csv", METADATA_PATH), ("genes.csv", GENE_PATH), ("count_gene_by_cell.mtx", COUNT_PATH), ("source_h5ad", H5AD_PATH)]:
        require(sha256(path) == run_info["input_sha256"][name], "Input changed during run: " + name)
    META_COLS = ["cell", "sample_name", "ZT", "celltype", "sigma_bio_hour", "T_true", "Theta_true"]
    evaluation_metadata = pd.read_csv(METADATA_PATH, usecols=META_COLS, dtype={"cell": str, "sample_name": str, "celltype": str})
    require(evaluation_metadata[DESIGN_COLS].equals(metadata[DESIGN_COLS]), "Design metadata changed")
    require(not evaluation_metadata[META_COLS].isna().any().any(), "Missing evaluation metadata")
    require(np.isfinite(evaluation_metadata[["ZT", "sigma_bio_hour", "T_true", "Theta_true"]].to_numpy(dtype=float)).all(), "Nonfinite evaluation metadata")
    evaluation_metadata = evaluation_metadata[META_COLS].copy()
    def wrap_hour(x):
        return (np.asarray(x) + 12) % 24 - 12


    require(sha256(OUT_DIR / "scPrisma_aligned_predictions.csv") == run_info["frozen_prediction_sha256"], "Frozen file modified")
    frozen = pd.read_csv(OUT_DIR / "scPrisma_aligned_predictions.csv", dtype={"cell":str})
    require(np.max(np.abs(wrap_hour(evaluation_metadata.Theta_true.to_numpy()*24/(2*np.pi)-evaluation_metadata.T_true.to_numpy()))) < 1e-7, "Truth radians/hours inconsistent")
    phase = evaluation_metadata.merge(frozen[["cell","pred_hour","pred_rad"]], on="cell", how="left", validate="one_to_one")
    require(len(phase)==9000 and phase.cell.equals(metadata.cell), "Evaluation merge changed coverage or order")
    phase["wrapped_error_hour"] = wrap_hour(phase.pred_hour-phase.T_true)
    phase["abs_error_hour"] = np.abs(phase.wrapped_error_hour)
    STANDARD_COLS = META_COLS + ["pred_hour", "pred_rad", "wrapped_error_hour", "abs_error_hour"]
    require(np.isfinite(phase[["pred_hour","pred_rad","wrapped_error_hour","abs_error_hour"]]).all().all(), "Nonfinite output")
    phase = phase[STANDARD_COLS]
    phase.to_csv(OUT_DIR / "scPrisma_cell_phase.csv", index=False)
    rows=[]
    for group in ["overall"] + CT_ORDER:
        d = phase if group == "overall" else phase[phase.celltype == group]
        for method, e in [("scPrisma",d.wrapped_error_hour.to_numpy()), ("ZT",wrap_hour(d.ZT-d.T_true))]:
            rows.append(dict(method=method, group=group, n_cells=len(d), circular_MAE_hour=float(np.mean(np.abs(e))),
                circular_MedAE_hour=float(np.median(np.abs(e))), circular_RMSE_hour=float(np.sqrt(np.mean(e**2)))))
    summary = pd.DataFrame(rows)
    summary[summary.method=="scPrisma"].to_csv(OUT_DIR / "scPrisma_circular_error_summary.csv", index=False)
    summary.to_csv(OUT_DIR / "scPrisma_vs_ZT_summary.csv", index=False)
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    for ax,ct in zip(axes.flat,CT_ORDER):
        d=phase[phase.celltype==ct]
        ax.scatter(d.T_true%24,d.pred_hour,s=3,alpha=.3)
        ax.plot([0,24],[0,24],"--",c="gray",lw=1)
        ax.set(title=ct,xlim=(0,24),ylim=(0,24),xlabel="True phase (h)",ylabel="Predicted phase (h)")
    fig.tight_layout(); fig.savefig(OUT_DIR / "true_vs_predicted_by_CT.png",dpi=180); plt.close(fig)
    run_info.update(status="evaluated", finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        prediction_sha256=sha256(OUT_DIR / "scPrisma_cell_phase.csv"), n_predictions=len(phase), output_columns=STANDARD_COLS)
    run_info["output_sha256"] = {str(p.relative_to(OUT_DIR)):sha256(p) for p in sorted(OUT_DIR.rglob("*"))
        if p.is_file() and p.name not in ["run_manifest.json", "run_manifest.json.tmp"]}
    save_json("run_manifest.json",run_info)
    print(summary.to_string(index=False))
    print("R input:",OUT_DIR / "scPrisma_cell_phase.csv")
except BaseException as exc:
    run_info.update(status="failed", error=repr(exc))
    save_json("run_manifest.json", run_info)
    raise


# ---- notebook code cell 8 ----
try:
    audit_rows = []
    for ct in CT_ORDER:
        ct_dir = OUT_DIR / ct
        with np.load(ct_dir / "reconstruction_matrices.npz") as arrays:
            E, E_rec = arrays["E"], arrays["E_rec"]
            q = json.loads((ct_dir / "QC.json").read_text())
            q.update(soft_bistochastic_tolerance=1e-3, rounded_zero_weight_threshold=1e-12,
                soft_bistochastic_ok=bool(max(q["E_row_sum_max_abs_deviation"],q["E_col_sum_max_abs_deviation"]) <= 1e-3),
                rounded_zero_weight_assignments=int(np.sum(E[E_rec == 1] <= 1e-12)),
                zero_soft_columns=int(np.sum(E.sum(axis=0) <= 1e-12)))
            q["numerical_attention"] = bool(not q["soft_bistochastic_ok"] or q["rounded_zero_weight_assignments"] > 0)
            audit_rows.append(q)
            save_json(ct + "/QC.json", q)
    qc_table = pd.DataFrame(audit_rows)
    qc_table.to_csv(OUT_DIR / "scPrisma_QC_by_celltype.csv", index=False)
    run_info.update(status="complete", celltype_qc=audit_rows,
        numerical_qc_status="attention" if qc_table.numerical_attention.any() else "passed",
        numerical_qc_note="Check soft E row/column sums and zero-weight greedy assignments; completed execution does not establish convergence.")
    run_info["output_sha256"] = {str(p.relative_to(OUT_DIR)):sha256(p) for p in sorted(OUT_DIR.rglob("*"))
        if p.is_file() and p.name not in ["run_manifest.json", "run_manifest.json.tmp"]}
    save_json("run_manifest.json", run_info)
    print(qc_table[["celltype","n_cells","E_row_sum_max_abs_deviation","E_col_sum_max_abs_deviation",
        "rounded_zero_weight_assignments","numerical_attention"]].to_string(index=False))
    print("Numerical QC:", run_info["numerical_qc_status"])
except BaseException as exc:
    run_info.update(status="failed", error=repr(exc))
    save_json("run_manifest.json", run_info)
    raise

