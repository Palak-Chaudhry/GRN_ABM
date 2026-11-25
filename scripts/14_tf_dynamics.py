"""Full dynamic-GRN view: every TF, every gene and every GRN edge over time, for every cell type.

Unlike plot 04a (which shows the 30 most dynamic TFs per type, a DIFFERENT set in each panel, so the
panels cannot be compared), every panel here uses the SAME rows in the SAME order and a shared colour
scale, so a row can be read across cell types.

  14a  all TFs, mean log2 vs each type's own baseline, one panel per type
  14b  all genes grouped by role (TF / ligand / receptor / module target)
  14c  GRN edge activity (regulator activation x edge weight), all edges, per type
  14d  L-R signalling PER RECEIVING CELL as well as whole-tissue total -- the total rises mostly
       because the tissue grows and shifts towards high-expressing tumour cells, which the per-cell
       panel divides out

  python scripts/14_tf_dynamics.py --run results/runs/<name> [--top-edges 60]

Output: <run>/plots/14*.png, <run>/tf_dynamics.tsv
"""
import argparse, os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import numpy as np
import pandas as pd
pd.options.future.infer_string = False
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True)
ap.add_argument("--top-edges", type=int, default=60, help="edges shown in 14c, by largest activity swing")
a = ap.parse_args()
RUN = a.run
OUT = f"{RUN}/plots"
os.makedirs(OUT, exist_ok=True)

A = ad.read_h5ad(f"{RUN}/snapshots.h5ad")
S = np.asarray(A.layers["S_true"])
TH = np.sort(A.obs.time_h.unique())
TYPES = [str(t) for t in np.load(f"{RUN}/grn_activity.npz")["types"]]
genes = list(A.var_names)
gi = {g: i for i, g in enumerate(genes)}

ann = pd.read_csv(f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/"
                  + str(pd.read_csv(f"{RUN}/config_engine.yaml", sep="\x00", header=None)[0]
                        .str.extract(r"network_dir:\s*(\S+)").dropna().iloc[0, 0]) + "/genes.tsv", sep="\t")
role = dict(zip(ann.gene, ann.role.fillna("")))


def mean_by_type(idx):
    """idx: gene indices. -> dict type -> (time x gene) mean log2 vs that type's t=0 mean."""
    out = {}
    for t in TYPES:
        M = np.array([np.log1p(S[((A.obs.time_h == x) & (A.obs.cell_type == t)).to_numpy()][:, idx].mean(0))
                      for x in TH])
        out[t] = M - M[0]
    return out


def panel(fig, ax, M, rows, title, vmax, cmap="RdBu_r", ylab=True):
    im = ax.imshow(M.T, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax, interpolation="nearest")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows if ylab else [""] * len(rows), fontsize=5.2)
    step = max(1, len(TH) // 6)
    ax.set_xticks(range(0, len(TH), step))
    ax.set_xticklabels([f"{x:g}" for x in TH[::step]], fontsize=8)
    ax.set_xlabel("time (h)", fontsize=9)
    ax.set_title(title, fontsize=10)
    return im


# ---------------------------------------------------------------- 14a: every TF
tf_names = sorted([g for g in genes if "TF" in role.get(g, "")])
tf_idx = [gi[g] for g in tf_names]
D = mean_by_type(tf_idx)
order = np.argsort(-np.abs(D[TYPES[0]]).max(0))                      # rank by swing in the first type
rows = [tf_names[i] for i in order]
vmax = max(np.abs(D[t]).max() for t in TYPES)
fig, axes = plt.subplots(1, len(TYPES), figsize=(4.6 * len(TYPES) + 1.6, max(7, len(rows) * 0.135)))
for k, t in enumerate(TYPES):
    im = panel(fig, axes[k], D[t][:, order], rows, f"{t}  (n TFs = {len(rows)})", vmax, ylab=(k == 0))
fig.suptitle(f"All {len(rows)} transcription factors, mean log2 vs each type's own t = 0 baseline "
             f"(same row order and colour scale in every panel)", fontsize=11)
fig.colorbar(im, ax=axes, shrink=0.55, label="log2 vs baseline")
fig.savefig(f"{OUT}/14a_all_tfs_over_time.png", dpi=140, bbox_inches="tight"); plt.close(fig)
print("wrote 14a_all_tfs_over_time.png")

pd.DataFrame({f"{t}__{TH[i]:g}h": D[t][i, order] for t in TYPES for i in range(len(TH))},
             index=rows).to_csv(f"{RUN}/tf_dynamics.tsv", sep="\t")

# ---------------------------------------------------------------- 14b: every gene, grouped by role
groups = [("TF", [g for g in genes if "TF" in role.get(g, "")]),
          ("ligand", [g for g in genes if "ligand" in role.get(g, "") and "TF" not in role.get(g, "")]),
          ("receptor", [g for g in genes if "receptor" in role.get(g, "")
                        and "TF" not in role.get(g, "") and "ligand" not in role.get(g, "")]),
          ("module target", [g for g in genes if "target" in role.get(g, "")
                             and not any(r in role.get(g, "") for r in ("TF", "ligand", "receptor"))])]
groups = [(n, gs) for n, gs in groups if gs]
hs = [len(gs) for _, gs in groups]
fig, axes = plt.subplots(len(groups), len(TYPES), figsize=(4.6 * len(TYPES) + 1.6, 0.11 * sum(hs) + 2.4),
                         gridspec_kw={"height_ratios": hs}, squeeze=False)
vmax_all = 0
Dg = []
for n, gs in groups:
    d = mean_by_type([gi[g] for g in gs])
    Dg.append((n, gs, d))
    vmax_all = max(vmax_all, max(np.abs(d[t]).max() for t in TYPES))
for r, (n, gs, d) in enumerate(Dg):
    o = np.argsort(-np.abs(d[TYPES[0]]).max(0))
    for c, t in enumerate(TYPES):
        im = panel(fig, axes[r][c], d[t][:, o], [gs[i] for i in o],
                   f"{t} - {n}" if r == 0 else f"{n}", vmax_all, ylab=(c == 0))
        if r < len(Dg) - 1:
            axes[r][c].set_xlabel("")
fig.suptitle(f"All {len(genes)} genes by role, mean log2 vs each type's own t = 0 baseline", fontsize=11)
fig.colorbar(im, ax=axes, shrink=0.4, label="log2 vs baseline")
fig.savefig(f"{OUT}/14b_all_genes_by_role.png", dpi=140, bbox_inches="tight"); plt.close(fig)
print("wrote 14b_all_genes_by_role.png")

# ---------------------------------------------------------------- 14c: GRN edge activity, every type
g_ = np.load(f"{RUN}/grn_activity.npz")
W, act = g_["W"], g_["tf_activity"]
tfs_a, genes_a = [str(x) for x in g_["tfs"]], [str(x) for x in g_["genes"]]
e_tf, e_g = np.nonzero(W)
E = {t: act[:, k, e_tf] * W[e_tf, e_g] for k, t in enumerate(TYPES)}
E = {t: v - v[:1] for t, v in E.items()}                                  # change from t = 0
swing = np.abs(E[TYPES[0]]).max(0)
top = np.argsort(-swing)[:a.top_edges]
lbl = [f"{tfs_a[e_tf[j]]}->{genes_a[e_g[j]]}" for j in top]
vmax = max(np.abs(E[t][:, top]).max() for t in TYPES)
fig, axes = plt.subplots(1, len(TYPES), figsize=(4.6 * len(TYPES) + 1.6, max(7, len(top) * 0.135)))
for k, t in enumerate(TYPES):
    im = panel(fig, axes[k], E[t][:, top], lbl, t, vmax, cmap="PuOr_r", ylab=(k == 0))
fig.suptitle(f"GRN edge activity (regulator activation x weight), change from t = 0 -- "
             f"{a.top_edges} most dynamic of {len(e_tf)} edges", fontsize=11)
fig.colorbar(im, ax=axes, shrink=0.55, label="change in edge input")
fig.savefig(f"{OUT}/14c_grn_edge_activity_all_types.png", dpi=140, bbox_inches="tight"); plt.close(fig)
print("wrote 14c_grn_edge_activity_all_types.png")

# ---------------------------------------------------------------- 14d: L-R per receiving cell
f = np.load(f"{RUN}/lr_flow.npz")
flow, pairs = f["flow"], [str(p) for p in f["pairs"]]
cnt = {t: np.array([max(1, int(((A.obs.time_h == x) & (A.obs.cell_type == t)).sum())) for x in TH]) for t in TYPES}
ncell = sum(cnt[t] for t in TYPES)
tot = flow.sum((2, 3))                                                    # time x pair, whole tissue
per = tot / ncell[:, None]                                                # per cell in the tissue
top = np.argsort(-tot.mean(0))[:20]
rel_tot = np.log2((tot[:, top] + 1e-9) / (tot[:1, top] + 1e-9))
rel_per = np.log2((per[:, top] + 1e-9) / (per[:1, top] + 1e-9))
fig, ax = plt.subplots(1, 3, figsize=(19, 6.4), gridspec_kw={"width_ratios": [1, 1, 1.05]})
for j, (M, ttl) in enumerate(((rel_tot, "whole-tissue TOTAL (as in plot 04d)"),
                              (rel_per, "PER CELL in the tissue"))):
    im = ax[j].imshow(M.T, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax[j].set_yticks(range(len(top))); ax[j].set_yticklabels([pairs[i] for i in top], fontsize=7)
    st = max(1, len(TH) // 6)
    ax[j].set_xticks(range(0, len(TH), st)); ax[j].set_xticklabels([f"{x:g}" for x in TH[::st]])
    ax[j].set(title=f"top 20 L-R pairs: {ttl}\nlog2 vs t=0", xlabel="time (h)")
    plt.colorbar(im, ax=ax[j], shrink=0.6)
for k, t in enumerate(TYPES):
    r = flow[:, :, :, k].sum((1, 2))
    ax[2].plot(TH, np.log2((r / cnt[t]) / (r[0] / cnt[t][0])), label=f"receiver = {t}", lw=2)
ax[2].axhline(0, color="k", lw=0.8)
ax[2].plot(TH, np.log2(ncell / ncell[0]), "k--", lw=1.2, label="log2 cell count (growth alone)")
ax[2].set(xlabel="time (h)", ylabel="log2 vs t = 0",
          title="signalling received PER CELL, by receiver type\n(growth divided out)")
ax[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/14d_lr_per_cell.png", dpi=140, bbox_inches="tight"); plt.close(fig)
print("wrote 14d_lr_per_cell.png")
print(f"wrote {RUN}/tf_dynamics.tsv")
