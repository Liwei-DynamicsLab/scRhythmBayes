# Generated from the named notebook; original file is unchanged.
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# ---- notebook code cell 0 ----
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

from CoPhaser import utils
from CoPhaser import gene_sets
from CoPhaser.model import CoPhaser
from CoPhaser.model import VAEModelLoader
from CoPhaser.trainer import Trainer
from CoPhaser.loss import Loss

# ---- notebook code cell 1 ----
import sys

print("Python:", sys.version)
print("Executable:", sys.executable)
print("PyTorch:", torch.__version__)
print("MPS:", torch.backends.mps.is_available())

# ---- notebook code cell 2 ----
RUN_DIR = Path("/path/to/analysis_home/Projects/scRhythmBayes_simulation/output/scRITMO_extended_simulation/sigma_gradient_1p5_2p5_3p5_4p5_5p5_6p5_empirical_fixed_v1")

H5AD_PATH = RUN_DIR / "scRITMO_extended_simulation.h5ad"

METADATA_PATH = RUN_DIR / "R_input" / "metadata.csv"

OUT_DIR = Path("/path/to/scRhythmBayes/reproducibility/07_runtime/runs/CoPhaser")

OUT_DIR.mkdir(parents=True, exist_ok=True)

H5AD_PATH, METADATA_PATH, OUT_DIR

# ---- notebook code cell 3 ----
SEED = 1

np.random.seed(SEED)
torch.manual_seed(SEED)

# ---- notebook code cell 4 ----
source = ad.read_h5ad(H5AD_PATH)

print(source.shape)

# ---- notebook code cell 5 ----
adata = ad.AnnData(
    X=source.X.copy(),
    obs=pd.DataFrame(index=source.obs_names.copy()),
    var=pd.DataFrame(index=source.var_names.copy()),
)

adata.layers["total"] = adata.X.copy()

adata

# ---- notebook code cell 6 ----
genes_tempo = [
    "Arntl",
    "Arntl2",
    "Bhlhe41",
    "Ciart",
    "Clock",
    "Cry1",
    "Cry2",
    "Csnk1a1",
    "Csnk1d",
    "Csnk1e",
    "Dbp",
    "Dec1",
    "Dec2",
    "Fbxo21",
    "Fbxo3",
    "Gm129",
    "Hlf",
    "Nampt",
    "Nfil3",
    "Npas2",
    "Nr1d1",
    "Nr1d2",
    "Per1",
    "Per2",
    "Per3",
    "Rora",
    "Rorc",
    "Tef",
]

gene_sets.SMALL_CIRCADIAN_GENE_SET

rhythmic_genes = np.array(
    sorted(
        set(list(gene_sets.SMALL_CIRCADIAN_GENE_SET) + genes_tempo)
        & set(adata.var_names)
    )
)

print("Rhythmic genes:", len(rhythmic_genes))
print(rhythmic_genes)

# ---- notebook code cell 7 ----
context_genes = utils.get_variable_genes(adata, 2000)

print("Context genes:", len(context_genes))

pd.Series(rhythmic_genes, name="gene").to_csv(
    OUT_DIR / "CoPhaser_rhythmic_genes.txt",
    index=False,
    header=False,
)

pd.Series(context_genes, name="gene").to_csv(
    OUT_DIR / "CoPhaser_context_genes.txt",
    index=False,
    header=False,
)



# ---- notebook code cell 8 ----
amp_phase_circadian = {
    "Bmal1": (1, 5.5),
    "Npas2": (0.5, 5.5),
    "Rorc": (0.5, 5.7),
    "Nr1d1": (1, 2.0),
    "Nr1d2": (0.75, 2.09),
    "Tef": (0.75, 2.35),
    "Ciart": (1, 2.35),
    "Dbp": (1, 2.09),
    "Per3": (1, 2.36),
    "Cry1": (0.5, 4.71),
    "Cry2": (0.25, 3),
    "Per2": (0.6, 3.14),
    "Per1": (0.5, 3.14),
    "Hlf": (0.5, 3.14),
}

amp_phase_circadian_used = {
    gene: value
    for gene, value in amp_phase_circadian.items()
    if gene in rhythmic_genes
}

print("Prior genes:", len(amp_phase_circadian_used))
print(amp_phase_circadian_used)

pd.DataFrame(
    [
        {"gene": gene, "amplitude_prior": value[0], "phase_prior": value[1]}
        for gene, value in amp_phase_circadian_used.items()
    ]
).to_csv(
    OUT_DIR / "CoPhaser_circadian_prior.csv",
    index=False,
)

# ---- notebook code cell 9 ----
import inspect

print("CoPhaser:")
print(inspect.signature(CoPhaser.__init__))

print("\nVAEModelLoader.define_decoder_prior:")
print(inspect.signature(VAEModelLoader.define_decoder_prior))

print("\nTrainer:")
print(inspect.signature(Trainer.__init__))

print("\nTrainer.train_model:")
print(inspect.signature(Trainer.train_model))

model = CoPhaser(
    rhythmic_genes,
    context_genes,
    n_latent=10,
    n_harm=1,
    rhythmic_decoder_to_all_genes=True,
    use_mu_z_encoder=True,
    z_range=20,
)

model.load_anndata(
    adata,
    layer_to_use="total",
)

VAEModelLoader.define_decoder_prior(
    amp_phase_prior=amp_phase_circadian_used,
    model=model,
)


# ---- notebook code cell 10 ----
model_parameters = inspect.signature(CoPhaser.__init__).parameters

requested_parameters = [
    "n_latent",
    "n_harm",
    "rhythmic_decoder_to_all_genes",
    "use_mu_z_encoder",
    "z_range",
]

print("Supported:")
for x in requested_parameters:
    print(x, x in model_parameters)

# ---- notebook code cell 11 ----
trainer = Trainer(
    model,
    Loss.compute_loss,
    non_rhythmic_likelihood_weight=2,
    rhythmic_likelihood_weight=10,
    unfreeze_epoch_layer=[(10, "rhythmic_decoder")],
    L2_Z_decoder_loss_weight=0,
    closed_circle_weight=0,
    MI_weight=50,
    entropy_weight_factor=50,
    cycling_status_prior=1,
    MI_detach="f",
)

# ---- notebook code cell 12 ----
trainer.train_model(
    n_epochs=200,
    lr=1e-2,
    device="cpu",
    batch_size=2048,
)

# ---- notebook code cell 13 ----
model.to("cpu")

generative_outputs, space_outputs = model.get_outputs()

mu_theta = space_outputs["mu_theta"].detach().cpu().numpy()

thetas = np.arctan2(
    mu_theta[:, 1],
    mu_theta[:, 0],
)

z_space = space_outputs["z"].detach().cpu().numpy()

print("theta shape:", thetas.shape)
print("z shape:", z_space.shape)
print("theta min:", thetas.min())
print("theta max:", thetas.max())


# ---- notebook code cell 14 ----
plt.figure(figsize=(7, 4))

sns.histplot(
    x=thetas,
    bins=100,
)

plt.xlabel("CoPhaser raw phase (rad)")
plt.ylabel("Cells")
plt.title("CoPhaser inferred circadian phase")

plt.tight_layout()
plt.show()

# ---- notebook code cell 15 ----
theta_rad_raw = np.mod(thetas, 2 * np.pi)

theta_hour_raw = theta_rad_raw * 24 / (2 * np.pi)

print("Raw phase hour min:", theta_hour_raw.min())
print("Raw phase hour max:", theta_hour_raw.max())

CoPhaser_phase = pd.DataFrame({
    "cell": adata.obs_names,
    "cophaser_theta_raw": thetas,
    "cophaser_theta_rad_raw": theta_rad_raw,
    "cophaser_theta_hour_raw": theta_hour_raw,
})
CoPhaser_phase.head()

CoPhaser_phase.to_csv(
    OUT_DIR / "CoPhaser_cell_phase_raw.csv",
    index=False,
)

np.save(
    OUT_DIR / "CoPhaser_context_z.npy",
    z_space,
)

# ---- notebook code cell 16 ----
metadata = pd.read_csv(METADATA_PATH)

metadata = metadata[
    [
        "cell",
        "sample_name",
        "ZT",
        "celltype",
        "sigma_bio_hour",
        "T_true",
        "Theta_true",
    ]
].copy()

CoPhaser_phase_full = metadata.merge(
    CoPhaser_phase,
    on="cell",
    how="left",
)

CoPhaser_phase_full.head()

# ---- notebook code cell 17 ----
plt.figure(figsize=(8, 4))

sns.histplot(
    data=CoPhaser_phase_full,
    x="cophaser_theta_raw",
    hue="ZT",
    bins=100,
    element="step",
    stat="density",
    common_norm=False,
)

plt.xlabel("CoPhaser raw phase (rad)")
plt.ylabel("Density")
plt.title("CoPhaser raw phase by ZT")

plt.tight_layout()
plt.show()


# ---- notebook code cell 18 ----
def optimal_shift(p, p0, n_s=200):
    def circular_deviation2(x, y, period=2 * np.pi):
        x, y = x % period, y % period
        v1 = np.abs(x - y)
        v2 = (period - v1) % period
        return np.minimum(v1, v2)

    Nc = p.shape[0]

    shifts = np.linspace(0, 2 * np.pi, n_s)

    theta_cs = (p.reshape(Nc, 1) - shifts.reshape(1, n_s)) % (2 * np.pi)

    theta_cs_neg = (-p.reshape(Nc, 1) - shifts.reshape(1, n_s)) % (2 * np.pi)

    delta_cs = circular_deviation2(p0[:, None], theta_cs)

    delta_cs_neg = circular_deviation2(p0[:, None], theta_cs_neg)

    v = np.median(delta_cs, axis=0)

    v_neg = np.median(delta_cs_neg, axis=0)

    best_shift_ind = np.argmin(v)

    best_shift_ind_neg = np.argmin(v_neg)

    mad = v[best_shift_ind]

    mad_neg = v_neg[best_shift_ind_neg]

    if mad < mad_neg:
        phi_aligned = theta_cs[:, best_shift_ind]
        best_mad = mad
    else:
        phi_aligned = theta_cs_neg[:, best_shift_ind_neg]
        best_mad = mad_neg

    return phi_aligned, best_mad

# ---- notebook code cell 19 ----
w = 2 * np.pi / 24

ZT_rad = (CoPhaser_phase_full["ZT"].to_numpy(dtype=float) % 24) * w

phi_raw = CoPhaser_phase_full["cophaser_theta_raw"].to_numpy() % (2 * np.pi)

pred_rad, alignment_mad_rad = optimal_shift(
    phi_raw,
    ZT_rad,
)

pred_hour = pred_rad * 24 / (2 * np.pi)

alignment_mad_hour = alignment_mad_rad / w

print("ZT alignment MAD:", alignment_mad_hour, "h")

# ---- notebook code cell 20 ----
CoPhaser_phase_full["pred_rad"] = pred_rad

CoPhaser_phase_full["pred_hour"] = pred_hour

CoPhaser_phase_full[
    [
        "cell",
        "ZT",
        "T_true",
        "cophaser_theta_hour_raw",
        "pred_hour",
    ]
].head(10)

# ---- notebook code cell 21 ----
CoPhaser_phase_full["wrapped_error_hour"] = (
    CoPhaser_phase_full["pred_hour"]
    - CoPhaser_phase_full["T_true"]
    + 12
) % 24 - 12

CoPhaser_phase_full["abs_error_hour"] = np.abs(
    CoPhaser_phase_full["wrapped_error_hour"]
)

MAE = CoPhaser_phase_full["abs_error_hour"].mean()

MedAE = CoPhaser_phase_full["abs_error_hour"].median()

RMSE = np.sqrt(
    np.mean(
        CoPhaser_phase_full["wrapped_error_hour"] ** 2
    )
)

print("Circular MAE:", MAE)
print("Circular MedAE:", MedAE)
print("Circular RMSE:", RMSE)

# ---- notebook code cell 22 ----
summary = []

for group in ["overall", "CT1", "CT2", "CT3", "CT4", "CT5", "CT6"]:

    if group == "overall":
        x = CoPhaser_phase_full
    else:
        x = CoPhaser_phase_full.loc[
            CoPhaser_phase_full["celltype"] == group
        ]

    e = x["wrapped_error_hour"].to_numpy()

    summary.append({
        "group": group,
        "n_cells": len(x),
        "circular_MAE_hour": np.mean(np.abs(e)),
        "circular_MedAE_hour": np.median(np.abs(e)),
        "circular_RMSE_hour": np.sqrt(np.mean(e ** 2)),
    })

# ---- notebook code cell 23 ----
CoPhaser_summary = pd.DataFrame(summary)

CoPhaser_summary

# ---- notebook code cell 24 ----
CoPhaser_phase_full.to_csv(
    OUT_DIR / "CoPhaser_cell_phase.csv",
    index=False,
)

CoPhaser_summary.to_csv(
    OUT_DIR / "CoPhaser_circular_error_summary.csv",
    index=False,
)

pd.DataFrame({
    "metric": ["ZT_alignment_MAD_hour"],
    "value": [alignment_mad_hour],
}).to_csv(
    OUT_DIR / "CoPhaser_alignment_summary.csv",
    index=False,
)

# ---- notebook code cell 25 ----
plt.figure(figsize=(6, 6))

plt.scatter(
    CoPhaser_phase_full["T_true"] % 24,
    CoPhaser_phase_full["pred_hour"] % 24,
    s=6,
    alpha=0.25,
)

plt.plot(
    [0, 24],
    [0, 24],
    linestyle="--",
)

plt.xlim(0, 24)
plt.ylim(0, 24)

plt.xticks(np.arange(0, 25, 4))
plt.yticks(np.arange(0, 25, 4))

plt.xlabel("True Time (h)")
plt.ylabel("Predicted Time (h)")
plt.title("CoPhaser")

plt.tight_layout()
plt.show()

# ---- notebook code cell 26 ----
print("Cells:", len(CoPhaser_phase_full))
print("Rhythmic genes:", len(rhythmic_genes))
print("Prior genes:", len(amp_phase_circadian_used))
print("Context genes:", len(context_genes))
print("ZT alignment MAD:", alignment_mad_hour)
print("Circular MAE:", MAE)
print("Circular MedAE:", MedAE)
print("Circular RMSE:", RMSE)
print("Result:", OUT_DIR / "CoPhaser_cell_phase.csv")
