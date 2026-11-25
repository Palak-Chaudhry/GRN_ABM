"""Export the first (t = 0) and last (t = end) snapshots of a run as standalone AnnData files.

Each keeps everything the snapshot carries: expression (X = observed spliced counts, plus the
S_true / U_true / spliced / unspliced layers), spatial coordinates (obsm['spatial'], and x / y in
obs), cell type, lineage (cell_id, parent_id, birth_h), per-cell module and rule values, signals,
neighbour counts, and the gene annotation in var.

  python scripts/11_export_endpoints.py --run results/runs/<name>

Output: <run>/adata_t0.h5ad, <run>/adata_tfinal.h5ad, <run>/endpoints_summary.md
"""
import argparse, os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import numpy as np
import pandas as pd
pd.options.future.infer_string = False
import anndata as ad
import scipy.sparse as sp

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True)
a = ap.parse_args()
log = []


def say(*x):
    log.append(" ".join(map(str, x))); print(log[-1])


A = ad.read_h5ad(f"{a.run}/snapshots.h5ad")
times = np.sort(A.obs.time_h.unique())
say(f"# Endpoint export: {os.path.basename(a.run.rstrip('/'))}\n")
say(f"- source: `{a.run}/snapshots.h5ad` -- {A.n_obs} cell-snapshots x {A.n_vars} genes, "
    f"{len(times)} snapshots from t = {times[0]:g} h to t = {times[-1]:g} h\n")

for tag, t in (("t0", times[0]), ("tfinal", times[-1])):
    S = A[(A.obs.time_h == t).to_numpy()].copy()
    S.uns = dict(A.uns)
    S.uns["time_h"] = float(t)
    S.uns["source_run"] = os.path.abspath(a.run)
    out = f"{a.run}/adata_{tag}.h5ad"
    S.write_h5ad(out)
    ct = S.obs.cell_type.value_counts()
    xy = S.obsm["spatial"]
    say(f"## t = {t:g} h -> `adata_{tag}.h5ad`")
    say(f"- {S.n_obs} cells x {S.n_vars} genes; " + ", ".join(f"{k} {v}" for k, v in ct.items()))
    say(f"- spatial extent: x [{xy[:, 0].min():.0f}, {xy[:, 0].max():.0f}] um, "
        f"y [{xy[:, 1].min():.0f}, {xy[:, 1].max():.0f}] um")
    tot = np.asarray(S.X.sum(1)).ravel() if sp.issparse(S.X) else np.asarray(S.X).sum(1)
    say(f"- mean observed counts per cell {tot.mean():.1f}; "
        f"dying {int(S.obs.dying.sum())}; mean neighbours {S.obs.neighbours.mean():.1f}")
    say(f"- layers: {', '.join(S.layers)} | obsm: {', '.join(S.obsm)} | obs: {S.obs.shape[1]} columns\n")

open(f"{a.run}/endpoints_summary.md", "w").write("\n".join(log) + "\n")
say(f"wrote {a.run}/adata_t0.h5ad, {a.run}/adata_tfinal.h5ad, {a.run}/endpoints_summary.md")
