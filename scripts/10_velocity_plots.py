"""RNA velocity plots for a recorded run (scvelo; run in the velo_py environment, after 09_analysis_plots.py
has written <run>/analysis/{final,pooled}.h5ad with their UMAPs).

For the last snapshot and for all snapshots pooled:
  - scvelo velocity estimated from the observed spliced / unspliced counts (stochastic model)
  - the simulator's TRUE velocity (beta u - gamma s of the true molecules, scaled to normalised counts)
  - streams of both on the same UMAP; accuracy of the estimate (per-cell cosine, per-gene correlation)
  - phase portraits of the most dynamic genes, with the true steady-state line (s = beta / gamma u)
  - velocity pseudotime (estimated and true) vs the true simulated time

Usage: /ihome/djishnu/pac271/.conda/envs/velo_py/bin/python scripts/10_velocity_plots.py --run results/runs/<name>
Output: <run>/plots/02*_velocity_*.png, <run>/plots/03b_velocity_pseudotime.png
"""
import argparse, os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import scvelo as scv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True)
a = ap.parse_args()
RUN = a.run if os.path.isabs(a.run) else f"{ROOT}/{a.run}"
PLOT = f"{RUN}/plots"
scv.settings.verbosity = 0


def dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x)


def prepare(path):
    b = ad.read_h5ad(path)
    types = list(b.obs.cell_type.cat.categories)
    cap = float(b.uns["capture"])
    for k in ("spliced", "unspliced"):
        b.layers[k] = sp.csr_matrix(dense(b.layers[k]))
    b.X = b.layers["spliced"].copy()
    total = np.asarray(b.layers["spliced"].sum(1)).ravel()
    scale = np.median(total) / np.maximum(total, 1)                       # scvelo's per-cell size factor
    v_true = b.var.beta_per_min.to_numpy() * dense(b.layers["U_true"]) - b.var.gamma_per_min.to_numpy() * dense(b.layers["S_true"])
    scv.pp.normalize_per_cell(b, enforce=True)
    b.X = b.X.log1p()                                                     # only for display; moments use X_pca
    scv.pp.moments(b, n_pcs=min(30, b.n_vars - 1), n_neighbors=30, use_rep="X_pca")
    scv.tl.velocity(b, mode="stochastic")
    b.layers["true_velocity"] = (cap * v_true * scale[:, None]).astype(np.float32)
    b.var["true_velocity_genes"] = True
    C = b.obsp["connectivities"].copy(); C.setdiag(1)
    C = sp.diags(1 / np.asarray(C.sum(1)).ravel()) @ C
    b.layers["true_velocity_smoothed"] = np.asarray(C @ b.layers["true_velocity"])
    b.var["true_velocity_smoothed_genes"] = True
    for vk in ("velocity", "true_velocity"):
        scv.tl.velocity_graph(b, vkey=vk, n_jobs=2)
        scv.tl.velocity_embedding(b, basis="umap", vkey=vk)
    return b, types


def accuracy(b):
    g = b.var.velocity_genes.to_numpy() if "velocity_genes" in b.var else np.ones(b.n_vars, bool)
    est = np.nan_to_num(b.layers["velocity"][:, g])
    res = {}
    for k in ("true_velocity", "true_velocity_smoothed"):
        tru = b.layers[k][:, g]
        cos = (est * tru).sum(1) / np.maximum(np.linalg.norm(est, axis=1) * np.linalg.norm(tru, axis=1), 1e-12)
        res[k] = cos
    gcor = np.array([np.corrcoef(est[:, j], b.layers["true_velocity_smoothed"][:, g][:, j])[0, 1]
                     if est[:, j].std() > 0 else np.nan for j in range(g.sum())])
    return res, pd.Series(gcor, index=b.var_names[g]), g


for key in ("final", "pooled"):
    f = f"{RUN}/analysis/{key}.h5ad"
    if not os.path.exists(f):
        print(f"missing {f}: run 09_analysis_plots.py first"); continue
    print(key, flush=True)
    b, types = prepare(f)
    colors = [plt.get_cmap("tab10")(i) for i in range(len(types))]
    b.uns["cell_type_colors"] = [matplotlib.colors.to_hex(c) for c in colors]
    # streams
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))
    for j, (vk, title) in enumerate((("velocity", "estimated by scvelo (stochastic) from counts"),
                                     ("true_velocity", "true velocity from the simulator"))):
        scv.pl.velocity_embedding_stream(b, basis="umap", vkey=vk, color="cell_type", ax=ax[j], show=False,
                                         legend_loc="right margin" if j else "none", title=title, size=8, alpha=0.6)
    save_as = f"{PLOT}/02a_velocity_stream_{key}.png"
    fig.savefig(save_as, dpi=130, bbox_inches="tight"); plt.close(fig); print("  wrote", save_as)
    # accuracy
    cos, gcor, g = accuracy(b)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4))
    for k, ls in (("true_velocity", "-"), ("true_velocity_smoothed", "--")):
        for i, t in enumerate(types):
            m = (b.obs.cell_type == t).to_numpy()
            ax[0].hist(cos[k][m], bins=40, range=(-1, 1), histtype="step", ls=ls, color=colors[i],
                       label=f"{t}{' (kNN-smoothed truth)' if 'smoothed' in k else ''}: median {np.median(cos[k][m]):+.2f}")
    ax[0].set(title="per-cell cosine(estimated, true velocity)", xlabel="cosine"); ax[0].legend(fontsize=6)
    gs = gcor.sort_values()
    ax[1].bar(range(len(gs)), gs.values, color="0.5"); ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set(title=f"per-gene corr(estimated, smoothed true), median {np.nanmedian(gs):+.2f}", xlabel="genes (sorted)", ylabel="Pearson")
    hl = b.var.halflife_h.reindex(gcor.index)
    ax[2].scatter(hl, gcor, s=8); ax[2].set_xscale("log")
    ax[2].set(title="accuracy vs mRNA half-life", xlabel="half-life (h)", ylabel="per-gene corr")
    fig.tight_layout(); fig.savefig(f"{PLOT}/02b_velocity_accuracy_{key}.png", dpi=130); plt.close(fig)
    print("  wrote", f"02b_velocity_accuracy_{key}.png")
    # phase portraits
    var_true = np.var(b.layers["true_velocity"], 0)
    top = np.argsort(-var_true)[:6]
    fig, ax = plt.subplots(1, 6, figsize=(24, 4))
    for j, gi in enumerate(top):
        Ms, Mu = b.layers["Ms"][:, gi], b.layers["Mu"][:, gi]
        for i, t in enumerate(types):
            m = (b.obs.cell_type == t).to_numpy()
            ax[j].scatter(Ms[m], Mu[m], s=2, color=colors[i], alpha=0.5, label=t)
        bt, gm = b.var.beta_per_min.iloc[gi], b.var.gamma_per_min.iloc[gi]
        xs = np.linspace(0, Ms.max(), 10)
        ax[j].plot(xs, gm / bt * xs, "k--", lw=1, label="true steady state")
        if "velocity_gamma" in b.var:
            ax[j].plot(xs, b.var.velocity_gamma.iloc[gi] * xs, "r:", lw=1, label="scvelo fit")
        ax[j].set(title=f"{b.var_names[gi]} (half-life {b.var.halflife_h.iloc[gi]:.1f} h)", xlabel="spliced (Ms)", ylabel="unspliced (Mu)")
    ax[0].legend(fontsize=6, markerscale=4)
    fig.tight_layout(); fig.savefig(f"{PLOT}/02c_phase_portraits_{key}.png", dpi=130); plt.close(fig)
    print("  wrote", f"02c_phase_portraits_{key}.png")
    # velocity pseudotime (pooled only: needs a time axis)
    if key == "pooled":
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for j, vk in enumerate(("velocity", "true_velocity")):
            vpt = np.full(b.n_obs, np.nan)
            for i, t in enumerate(types):
                m = (b.obs.cell_type == t).to_numpy()
                if m.sum() < 50:
                    continue
                c = b[m].copy()
                scv.pp.neighbors(c, n_neighbors=30, use_rep="X_pca")      # subset graph must be rebuilt
                scv.tl.velocity_graph(c, vkey=vk, n_jobs=2)
                scv.tl.velocity_pseudotime(c, vkey=vk)
                vpt[m] = c.obs[f"{vk}_pseudotime"].to_numpy()
                rho = spearmanr(vpt[m], b.obs.time_h[m]).correlation
                mean = pd.Series(vpt[m]).groupby(b.obs.time_h[m].to_numpy()).mean()
                ax[j].plot(mean.index, mean.values, color=colors[i], label=f"{t}: Spearman {rho:+.2f}")
            ax[j].set(title=f"{'estimated' if vk == 'velocity' else 'true'} velocity pseudotime vs true time (per type)",
                      xlabel="simulated time (h)", ylabel="mean velocity pseudotime"); ax[j].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(f"{PLOT}/03b_velocity_pseudotime.png", dpi=130); plt.close(fig)
        print("  wrote 03b_velocity_pseudotime.png")
print("done")
