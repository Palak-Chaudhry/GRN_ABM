"""Analysis plots for a recorded run (scripts/08_run_simulation.py). Each plot is one function registered in
PLOTS; run all of them or a subset with --plots. Velocity plots are in scripts/10_velocity_plots.py
(scvelo lives in another environment); this script writes the embeddings that script reuses.

  01_umap            UMAP of the last snapshot and of all snapshots pooled (by type, by time)
  03_pseudotime      diffusion pseudotime per type vs the true simulated time
  04_grn_dynamics    TF expression, module scores, GRN edge activity and L-R flow over time
  05_emergent        population, behaviour rates, mixing, infiltration, spread, pressure, motility (MSD),
                     neighbourhood enrichment
  06_spatial_types   tissue snapshots coloured by cell type
  07_movement        cell tracks and an animation (GIF)
  08_spatial_genes   gene and module maps at start / middle / end
  09_states          emergent cell states: clusters within each type, their share over time, behaviour,
                     niche and marker genes
  10_events          where and when divisions and deaths happen
  11_pressure        mechanical pressure maps and its link to division
  12_svg             spatially variable genes (Moran's I) over time
  13_drivers         genes most associated with each behaviour, vs the genes the rules actually read

Usage: python scripts/09_analysis_plots.py --run results/runs/<name> [--plots all|01_umap,05_emergent]
Output: <run>/plots/*.png, <run>/movement.gif, <run>/analysis/{final,pooled}.h5ad
"""
import argparse, os, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import numpy as np
import pandas as pd
pd.options.future.infer_string = False
import anndata as ad
import scanpy as sc
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True)
ap.add_argument("--plots", default="all")
ap.add_argument("--pool-cells", type=int, default=30000, help="cells sampled across all snapshots for pooled UMAP")
a = ap.parse_args()
RUN = a.run if os.path.isabs(a.run) else f"{ROOT}/{a.run}"
PLOT, ANA = f"{RUN}/plots", f"{RUN}/analysis"
os.makedirs(PLOT, exist_ok=True); os.makedirs(ANA, exist_ok=True)
sc.settings.verbosity = 0

A = ad.read_h5ad(f"{RUN}/snapshots.h5ad")
TYPES = list(A.uns["cell_types"])
COL = {t: plt.get_cmap("tab10")(i) for i, t in enumerate(TYPES)}
TIMES = np.sort(A.obs.time_h.unique())
ABM = yaml.safe_load(open(f"{RUN}/config_abm.yaml"))
RAD = {t: float((ABM.get("types") or {}).get(t, {}).get("radius", ABM["defaults"]["radius"])) for t in TYPES}
ADH = {t: float((ABM.get("types") or {}).get(t, {}).get("max_adhesion_distance", ABM["defaults"]["max_adhesion_distance"])) for t in TYPES}
print(f"{A.n_obs} cell-snapshots, {len(TIMES)} snapshots ({TIMES[0]:g}-{TIMES[-1]:g} h), types {TYPES}", flush=True)


def save(fig, name):
    fig.savefig(f"{PLOT}/{name}.png", dpi=130, bbox_inches="tight"); plt.close(fig)
    print("  wrote", name, flush=True)


def snap(t):
    return A[A.obs.time_h == t]


def pick_times(k=4):
    return TIMES[np.unique(np.linspace(0, len(TIMES) - 1, k).round().astype(int))]


def contacts(obs):
    """Contact pairs (i, j) of one snapshot: distance < f_i R_i + f_j R_j, living cells only."""
    xy = obs[["x", "y"]].to_numpy()
    t = obs.cell_type.astype(str).to_numpy()
    r = np.array([RAD[x] * ADH[x] for x in t])
    pr = cKDTree(xy).query_pairs(2 * r.max(), output_type="ndarray")
    if not len(pr):
        return pr
    d = np.linalg.norm(xy[pr[:, 0]] - xy[pr[:, 1]], axis=1)
    alive = ~obs.dying.to_numpy()
    keep = (d < r[pr[:, 0]] + r[pr[:, 1]]) & alive[pr[:, 0]] & alive[pr[:, 1]]
    return pr[keep]


def neighbour_fractions():
    """Adds nb_<type> (fraction of a cell's contacts of that type) and nb_total to A.obs, per snapshot."""
    if "nb_total" in A.obs:
        return
    cols = {f"nb_{t}": np.zeros(A.n_obs) for t in TYPES}
    tot = np.zeros(A.n_obs)
    for t in TIMES:
        idx = np.where(A.obs.time_h.to_numpy() == t)[0]
        o = A.obs.iloc[idx]
        pr = contacts(o)
        ti = o.cell_type.cat.codes.to_numpy()
        cnt = np.zeros((len(o), len(TYPES)))
        if len(pr):
            np.add.at(cnt, (pr[:, 0], ti[pr[:, 1]]), 1); np.add.at(cnt, (pr[:, 1], ti[pr[:, 0]]), 1)
        n = cnt.sum(1)
        tot[idx] = n
        for k, ty in enumerate(TYPES):
            cols[f"nb_{ty}"][idx] = np.where(n > 0, cnt[:, k] / np.maximum(n, 1), np.nan)
    for k, v in cols.items():
        A.obs[k] = v
    A.obs["nb_total"] = tot


def _drop_none(d):
    for k in list(d):
        if d[k] is None:
            del d[k]
        elif isinstance(d[k], dict):
            _drop_none(d[k])


def write(b, f):
    """None values in .uns (scanpy params) cannot be read by older anndata (velo_py env): drop them."""
    _drop_none(b.uns)
    b.write_h5ad(f)


def embed(adata, key):
    """log-normalise counts, PCA, kNN, UMAP (cached in analysis/<key>.h5ad)."""
    f = f"{ANA}/{key}.h5ad"
    if os.path.exists(f):
        return ad.read_h5ad(f)
    b = adata.copy()
    b.layers["counts"] = b.X.copy()
    sc.pp.normalize_total(b, target_sum=float(np.median(np.asarray(b.X.sum(1)).ravel())))
    sc.pp.log1p(b)
    sc.pp.pca(b, n_comps=min(30, b.n_vars - 1))
    sc.pp.neighbors(b, n_neighbors=15)
    sc.tl.umap(b)
    write(b, f)
    return b


def final_emb():
    return embed(snap(TIMES[-1]).copy(), "final")


def pooled_emb():
    f = f"{ANA}/pooled.h5ad"
    if os.path.exists(f):
        return ad.read_h5ad(f)
    rng = np.random.default_rng(0)
    per = max(1, a.pool_cells // len(TIMES))
    idx = np.concatenate([rng.choice(np.where(A.obs.time_h.to_numpy() == t)[0],
                                     min(per, int((A.obs.time_h == t).sum())), replace=False) for t in TIMES])
    return embed(A[np.sort(idx)].copy(), "pooled")


def scatter_type(ax, obs, s=3, alpha=0.8):
    for t in TYPES:
        m = (obs.cell_type == t).to_numpy() & ~obs.dying.to_numpy()
        ax.scatter(obs.x[m], obs.y[m], s=s * (RAD[t] / 8) ** 2, color=COL[t], label=t, alpha=alpha, lw=0)
    m = obs.dying.to_numpy()
    ax.scatter(obs.x[m], obs.y[m], s=s, color="k", marker="x", lw=0.4, label="dying", alpha=0.6)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])


PLOTS = {}


def plot(fn):
    PLOTS[fn.__name__[1:]] = fn
    return fn


# ====================================================================== plots
@plot
def _01_umap():
    f, p = final_emb(), pooled_emb()
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    for t in TYPES:
        m = (f.obs.cell_type == t).to_numpy()
        ax[0].scatter(*f.obsm["X_umap"][m].T, s=2, color=COL[t], label=t)
        m = (p.obs.cell_type == t).to_numpy()
        ax[1].scatter(*p.obsm["X_umap"][m].T, s=1, color=COL[t], label=t)
    o = np.random.default_rng(0).permutation(p.n_obs)
    im = ax[2].scatter(*p.obsm["X_umap"][o].T, s=1, c=p.obs.time_h.to_numpy()[o], cmap="viridis")
    plt.colorbar(im, ax=ax[2], label="time (h)")
    ax[0].set_title(f"last snapshot ({TIMES[-1]:g} h), {f.n_obs} cells"); ax[1].set_title(f"all snapshots pooled ({p.n_obs} cells)")
    ax[2].set_title("pooled, coloured by time"); ax[0].legend(markerscale=5, fontsize=8)
    for x in ax:
        x.set_xticks([]); x.set_yticks([]); x.set_xlabel("UMAP1"); x.set_ylabel("UMAP2")
    save(fig, "01_umap")


@plot
def _03_pseudotime():
    p = pooled_emb()
    dpt = np.full(p.n_obs, np.nan)
    fig, ax = plt.subplots(1, len(TYPES) + 1, figsize=(4.4 * (len(TYPES) + 1), 4))
    for k, t in enumerate(TYPES):
        m = (p.obs.cell_type == t).to_numpy()
        if m.sum() < 50:
            continue
        b = p[m].copy()
        sc.pp.neighbors(b, n_neighbors=15, use_rep="X_pca")
        sc.tl.diffmap(b)
        t0 = np.where(b.obs.time_h.to_numpy() == TIMES[0])[0]
        dc = b.obsm["X_diffmap"][:, 1]
        b.uns["iroot"] = int(t0[np.argmin(np.abs(dc[t0] - np.median(dc[t0])))])   # typical t = 0 cell
        sc.tl.dpt(b)
        dpt[m] = b.obs.dpt_pseudotime.to_numpy()
        rho = spearmanr(b.obs.dpt_pseudotime, b.obs.time_h).correlation
        ax[k + 1].scatter(b.obs.time_h + np.random.default_rng(0).uniform(-0.6, 0.6, b.n_obs), b.obs.dpt_pseudotime,
                          s=1, color=COL[t], alpha=0.4)
        mean = b.obs.groupby("time_h", observed=True).dpt_pseudotime.mean()
        ax[k + 1].plot(mean.index, mean.values, "k-", lw=1.5)
        ax[k + 1].set(title=f"{t}: Spearman(pseudotime, true time) = {rho:+.2f}", xlabel="simulated time (h)", ylabel="pseudotime")
    im = ax[0].scatter(*p.obsm["X_umap"].T, s=1, c=dpt, cmap="magma")
    plt.colorbar(im, ax=ax[0], label="diffusion pseudotime (per type)")
    ax[0].set(title="pseudotime on the pooled UMAP", xticks=[], yticks=[])
    fig.suptitle("Diffusion pseudotime vs true simulated time (cell types are fixed, so pseudotime can only "
                 "track gradual drift of state, not differentiation)", fontsize=9)
    save(fig, "03_pseudotime")


@plot
def _04_grn_dynamics():
    g = np.load(f"{RUN}/grn_activity.npz")
    base = np.log2(A.uns["baseline_counts"] / A.uns["capture"] + 1)                # type x gene, molecules
    S = A.layers["S_true"]
    ti = A.obs.cell_type.cat.codes.to_numpy(); tt = A.obs.time_h.to_numpy()
    tf_mask = A.var.is_tf.to_numpy()
    # (a) TF expression over time, per type
    fig, ax = plt.subplots(1, len(TYPES), figsize=(5.2 * len(TYPES), 7))
    for k, t in enumerate(TYPES):
        M = np.array([np.log2(S[(ti == k) & (tt == x)].mean(0) + 1) - base[k] if ((ti == k) & (tt == x)).any()
                      else np.full(A.n_vars, np.nan) for x in TIMES])          # time x gene
        Mt = M[:, tf_mask]
        top = np.argsort(-np.nanstd(Mt, 0))[:30]
        im = ax[k].imshow(Mt[:, top].T, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
        ax[k].set_yticks(range(len(top))); ax[k].set_yticklabels(A.var_names[tf_mask][top], fontsize=7)
        ax[k].set_xticks(range(0, len(TIMES), max(1, len(TIMES) // 6)))
        ax[k].set_xticklabels([f"{x:g}" for x in TIMES[::max(1, len(TIMES) // 6)]])
        ax[k].set(title=f"{t}: 30 most dynamic TFs", xlabel="time (h)")
    plt.colorbar(im, ax=ax, label="mean log2 vs baseline", shrink=0.6)
    save(fig, "04a_tf_expression_over_time")
    # (b) module scores over time
    mods = [c for c in A.obs if c.startswith("module_")]
    fig, ax = plt.subplots(1, len(mods), figsize=(3.6 * len(mods), 3.2), sharey=True)
    for j, mcol in enumerate(mods):
        for k, t in enumerate(TYPES):
            q = A.obs[A.obs.cell_type == t].groupby("time_h")[mcol].quantile([0.25, 0.5, 0.75]).unstack()
            ax[j].plot(q.index, q[0.5], color=COL[t], label=t)
            ax[j].fill_between(q.index, q[0.25], q[0.75], color=COL[t], alpha=0.15)
        ax[j].axhline(0, color="grey", lw=0.5); ax[j].set(title=mcol[7:], xlabel="time (h)")
    ax[0].set_ylabel("module score (log2 vs baseline; median, IQR)"); ax[0].legend(fontsize=7)
    save(fig, "04b_module_scores_over_time")
    # (c) GRN edge activity (TF activation x signed weight) over time
    W, act, tfs, genes = g["W"], g["tf_activity"], g["tfs"], g["genes"]
    e_tf, e_g = np.nonzero(W)
    fig, ax = plt.subplots(1, len(TYPES), figsize=(5.2 * len(TYPES), 7))
    for k, t in enumerate(TYPES):
        E = act[:, k, e_tf] * W[e_tf, e_g]                                         # time x edge
        ch = np.nan_to_num(E - E[:1])
        top = np.argsort(-np.abs(ch).max(0))[:30]
        im = ax[k].imshow(ch[:, top].T, aspect="auto", cmap="PuOr_r", vmin=-np.abs(ch[:, top]).max(), vmax=np.abs(ch[:, top]).max())
        ax[k].set_yticks(range(len(top))); ax[k].set_yticklabels([f"{tfs[e_tf[i]]}->{genes[e_g[i]]}" for i in top], fontsize=6.5)
        ax[k].set_xticks(range(0, len(g["time_h"]), max(1, len(g["time_h"]) // 6)))
        ax[k].set_xticklabels([f"{x:g}" for x in g["time_h"][::max(1, len(g["time_h"]) // 6)]])
        ax[k].set(title=f"{t}: 30 most changing GRN edges", xlabel="time (h)")
        plt.colorbar(im, ax=ax[k], label="edge input change vs t=0", shrink=0.6)
    save(fig, "04c_grn_edge_activity_over_time")
    # (d) L-R communication over time
    f = np.load(f"{RUN}/lr_flow.npz")
    flow, pairs = f["flow"], f["pairs"]                                            # time x pair x sender x receiver
    tot = flow.sum((2, 3))
    top = np.argsort(-tot.mean(0))[:20]
    fig, ax = plt.subplots(1, 3, figsize=(18, 6), gridspec_kw={"width_ratios": [1.3, 1, 1]})
    rel = np.log2((tot[:, top] + 1e-9) / (tot[:1, top] + 1e-9))
    im = ax[0].imshow(rel.T, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax[0].set_yticks(range(len(top))); ax[0].set_yticklabels(pairs[top], fontsize=7)
    ax[0].set_xticks(range(0, len(f["time_h"]), max(1, len(f["time_h"]) // 6)))
    ax[0].set_xticklabels([f"{x:g}" for x in f["time_h"][::max(1, len(f["time_h"]) // 6)]])
    ax[0].set(title="top 20 L-R pairs: total signalling, log2 vs t=0", xlabel="time (h)")
    plt.colorbar(im, ax=ax[0], shrink=0.6)
    for j, (k, lab) in enumerate(((0, "start"), (-1, "end"))):
        m = flow[k].sum(0)
        im = ax[j + 1].imshow(m, cmap="Greens")
        ax[j + 1].set_xticks(range(len(TYPES))); ax[j + 1].set_xticklabels(TYPES)
        ax[j + 1].set_yticks(range(len(TYPES))); ax[j + 1].set_yticklabels(TYPES)
        ax[j + 1].set(title=f"sender (row) -> receiver (col), all pairs, {lab} ({f['time_h'][k]:g} h)")
        for (r, c), v in np.ndenumerate(m):
            ax[j + 1].text(c, r, f"{v:.0f}", ha="center", va="center", fontsize=8)
    save(fig, "04d_lr_communication_over_time")


@plot
def _05_emergent():
    neighbour_fractions()
    o = A.obs
    tr = np.load(f"{RUN}/tracks.npz")
    fig, ax = plt.subplots(2, 4, figsize=(18, 8))
    ax = ax.ravel()
    living = o[~o.dying]
    for t in TYPES:
        L = living[living.cell_type == t].groupby("time_h")
        ax[0].plot(L.size().index, L.size().values, color=COL[t], label=t)
        ax[1].plot(L.rule_division_rate.mean().index, 60 * L.rule_division_rate.mean().values, color=COL[t])
        ax[2].plot(L.rule_apoptosis_rate.mean().index, 60 * L.rule_apoptosis_rate.mean().values, color=COL[t])
        ax[3].plot(L.rule_speed.mean().index, L.rule_speed.mean().values, color=COL[t])
        ax[4].plot(L[f"nb_{t}"].mean().index, L[f"nb_{t}"].mean().values, color=COL[t])
        r = np.hypot(living.x, living.y)[living.cell_type == t].groupby(living.time_h[living.cell_type == t])
        ax[6].plot(r.median().index, r.median().values, color=COL[t])
    if "tumor" in TYPES:
        for t in [x for x in TYPES if x != "tumor"]:
            L = living[living.cell_type == t].groupby("time_h")
            ax[5].plot(L.nb_tumor.apply(lambda v: np.nanmean(v > 0)).index,
                       L.nb_tumor.apply(lambda v: np.nanmean(v > 0)).values, color=COL[t], label=t)
        ax[5].set(title="infiltration: fraction touching tumour cells", xlabel="time (h)"); ax[5].legend(fontsize=7)
    # MSD from tracks
    for k, t in enumerate(TYPES):
        m = tr["type"] == k
        df = pd.DataFrame({"t": tr["time_h"][m], "id": tr["cell_id"][m], "x": tr["x"][m], "y": tr["y"][m]})
        piv_x, piv_y = df.pivot(index="t", columns="id", values="x"), df.pivot(index="t", columns="id", values="y")
        lags = np.arange(1, min(len(piv_x), 61))
        msd = [np.nanmean((piv_x.values[l:] - piv_x.values[:-l]) ** 2 + (piv_y.values[l:] - piv_y.values[:-l]) ** 2) for l in lags]
        dt = np.diff(piv_x.index[:2])[0] * 60
        ax[7].loglog(lags * dt, msd, color=COL[t])
    ax[0].set(title="living cells", xlabel="time (h)"); ax[0].legend(fontsize=7)
    ax[1].set(title="mean division rate (per h)", xlabel="time (h)")
    ax[2].set(title="mean death rate (per h)", xlabel="time (h)")
    ax[3].set(title="mean speed set by rules (um/min)", xlabel="time (h)")
    ax[4].set(title="segregation: fraction of contacts of own type", xlabel="time (h)")
    ax[6].set(title="median distance from tissue centre (um)", xlabel="time (h)")
    ax[7].set(title="mean squared displacement", xlabel="lag (min)", ylabel="um^2")
    fig.tight_layout(); save(fig, "05a_emergent_behaviour")
    # neighbourhood enrichment: observed / expected contacts, start vs end (expected = shuffled labels)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    rng = np.random.default_rng(0)
    for j, t in enumerate((TIMES[0], TIMES[-1])):
        ob = o[o.time_h == t]
        pr = contacts(ob)
        ti = ob.cell_type.cat.codes.to_numpy()
        def cm(lab):
            c = np.zeros((len(TYPES), len(TYPES)))
            np.add.at(c, (lab[pr[:, 0]], lab[pr[:, 1]]), 1)
            return c + c.T
        obs_c = cm(ti)
        perm = np.array([cm(rng.permutation(ti)) for _ in range(50)])
        z = (obs_c - perm.mean(0)) / np.maximum(perm.std(0), 1e-9)
        im = ax[j].imshow(z, cmap="RdBu_r", vmin=-np.abs(z).max(), vmax=np.abs(z).max())
        for (r, c), v in np.ndenumerate(z):
            ax[j].text(c, r, f"{v:+.0f}", ha="center", va="center", fontsize=9)
        ax[j].set_xticks(range(len(TYPES))); ax[j].set_xticklabels(TYPES); ax[j].set_yticks(range(len(TYPES))); ax[j].set_yticklabels(TYPES)
        ax[j].set(title=f"neighbourhood enrichment z, {t:g} h")
        plt.colorbar(im, ax=ax[j], shrink=0.7)
    save(fig, "05b_neighbourhood_enrichment")


@plot
def _06_spatial_types():
    ts = pick_times(6)
    fig, ax = plt.subplots(1, len(ts), figsize=(4.2 * len(ts), 4.4))
    lim = np.abs(A.obs[["x", "y"]].to_numpy()).max() * 1.02
    for j, t in enumerate(ts):
        ob = A.obs[A.obs.time_h == t]
        scatter_type(ax[j], ob)
        ax[j].set(title=f"{t:g} h ({(~ob.dying).sum()} living)", xlim=(-lim, lim), ylim=(-lim, lim))
    ax[0].legend(markerscale=3, fontsize=7, loc="lower left")
    save(fig, "06_spatial_cell_types")


@plot
def _07_movement():
    tr = np.load(f"{RUN}/tracks.npz")
    df = pd.DataFrame({k: tr[k] for k in ("time_h", "cell_id", "type", "x", "y", "dying")})
    lim = np.abs(df[["x", "y"]].to_numpy()).max() * 1.02
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(1, len(TYPES), figsize=(5.2 * len(TYPES), 5.2))
    last = df[df.time_h == df.time_h.max()]
    for k, t in enumerate(TYPES):
        ax[k].scatter(last.x, last.y, s=1, color="0.85")
        ids0 = df[(df.time_h == 0) & (df.type == k)].cell_id.unique()
        for cid in rng.choice(ids0, min(40, len(ids0)), replace=False):
            d = df[df.cell_id == cid]
            ax[k].plot(d.x, d.y, lw=0.7, color=COL[t], alpha=0.8)
            ax[k].plot(d.x.iloc[-1], d.y.iloc[-1], "o", ms=2, color="k")
        ax[k].set(title=f"{t}: 40 tracks over {df.time_h.max():g} h (dot = end)", xlim=(-lim, lim), ylim=(-lim, lim))
        ax[k].set_aspect("equal"); ax[k].set_xticks([]); ax[k].set_yticks([])
    save(fig, "07a_cell_tracks")
    # animation
    import imageio.v2 as imageio
    times = np.sort(df.time_h.unique())
    times = times[::max(1, len(times) // 120)]
    frames = []
    for t in times:
        d = df[df.time_h == t]
        fig, axx = plt.subplots(figsize=(5, 5), dpi=90)
        for k, ty in enumerate(TYPES):
            m = (d.type == k) & ~d.dying
            axx.scatter(d.x[m], d.y[m], s=2.2 * (RAD[ty] / 8) ** 2, color=COL[ty], lw=0, label=ty)
        axx.scatter(d.x[d.dying], d.y[d.dying], s=2, color="k", marker="x", lw=0.4)
        axx.set(xlim=(-lim, lim), ylim=(-lim, lim), title=f"t = {t:5.1f} h"); axx.set_aspect("equal")
        axx.set_xticks([]); axx.set_yticks([]); axx.legend(loc="lower left", fontsize=6, markerscale=2)
        fig.canvas.draw()
        frames.append(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()); plt.close(fig)
    imageio.mimsave(f"{RUN}/movement.gif", frames, duration=0.08, loop=0)
    print("  wrote movement.gif", flush=True)


def morans_i(o, X):
    xy = o[["x", "y"]].to_numpy()
    k = min(7, len(xy))
    _, nn = cKDTree(xy).query(xy, k)
    Z = (X - X.mean(0)) / np.maximum(X.std(0), 1e-12)
    WZ = Z[nn[:, 1:]].mean(1)                                                     # row-standardised kNN weights
    return (Z * WZ).sum(0) / np.maximum((Z * Z).sum(0), 1e-12)


@plot
def _08_spatial_genes():
    last = snap(TIMES[-1])
    logS = np.log2(np.asarray(last.layers["S_true"]) + 1)
    I = morans_i(last.obs, logS)
    rules = yaml.safe_load(open(f"{RUN}/config_rules.yaml")) if os.path.exists(f"{RUN}/config_rules.yaml") else {"rules": []}
    named = []
    for r in rules.get("rules") or []:
        for x in (r.get("signal", ""), r.get("behaviour", "").split(":")[-1]):
            if x in A.var_names and x not in named:
                named.append(x)
    svg = [g for g in A.var_names[np.argsort(-I)] if g not in named][:3]
    feats = named[:4] + svg + [c for c in A.obs if c.startswith("module_")][:3]
    ts = pick_times(3)
    lim = np.abs(A.obs[["x", "y"]].to_numpy()).max() * 1.02
    fig, ax = plt.subplots(len(feats), len(ts), figsize=(4.0 * len(ts), 3.6 * len(feats)))
    for i, f in enumerate(feats):
        vals = {t: (snap(t).obs[f].to_numpy() if f in A.obs else
                    np.log2(np.asarray(snap(t)[:, f].layers["S_true"]).ravel() + 1)) for t in ts}
        lo, hi = np.percentile(np.concatenate(list(vals.values())), [1, 99])
        for j, t in enumerate(ts):
            ob = snap(t).obs
            im = ax[i, j].scatter(ob.x, ob.y, s=2, c=vals[t], cmap="viridis" if f in A.obs else "magma", vmin=lo, vmax=hi, lw=0)
            ax[i, j].set(xlim=(-lim, lim), ylim=(-lim, lim)); ax[i, j].set_aspect("equal")
            ax[i, j].set_xticks([]); ax[i, j].set_yticks([])
            tag = "rule gene" if f in named else "spatially variable" if f in svg else "module score"
            ax[i, j].set_title(f"{f} ({tag}), {t:g} h", fontsize=9)
        plt.colorbar(im, ax=ax[i, -1], shrink=0.8, label="log2(mRNA+1)" if f not in A.obs else "log2 vs baseline")
    save(fig, "08_spatial_gene_expression")


@plot
def _09_states():
    neighbour_fractions()
    p = pooled_emb()
    for c in [c for c in A.obs if c.startswith("nb_")]:
        p.obs[c] = A.obs.loc[p.obs_names, c].to_numpy()
    state = pd.Series("", index=p.obs_names)
    for t in TYPES:
        m = (p.obs.cell_type == t).to_numpy()
        if m.sum() < 50:
            continue
        b = p[m].copy()
        sc.pp.neighbors(b, n_neighbors=15, use_rep="X_pca")
        sc.tl.leiden(b, resolution=0.4, key_added="st", flavor="igraph", n_iterations=2, directed=False)
        state[b.obs_names] = [f"{t}_s{x}" for x in b.obs.st]
    p.obs["state"] = pd.Categorical(state)
    write(p, f"{ANA}/pooled.h5ad")
    states = list(p.obs.state.cat.categories)
    cmap = plt.get_cmap("tab20")
    scol = {s: cmap(i % 20) for i, s in enumerate(states)}
    fig = plt.figure(figsize=(20, 10))
    gs = fig.add_gridspec(2, 4)
    ax = fig.add_subplot(gs[0, 0])
    for s in states:
        m = (p.obs.state == s).to_numpy()
        ax.scatter(*p.obsm["X_umap"][m].T, s=1, color=scol[s], label=s)
    ax.set(title="emergent states (Leiden within each type)", xticks=[], yticks=[]); ax.legend(markerscale=6, fontsize=6, ncol=2)
    # share over time
    for k, t in enumerate(TYPES[:3]):
        ax = fig.add_subplot(gs[0, 1 + k])
        frac = pd.crosstab(p.obs.time_h[p.obs.cell_type == t], p.obs.state[p.obs.cell_type == t], normalize="index")
        frac = frac.loc[:, frac.sum() > 0]
        ax.stackplot(frac.index, frac.T.values, colors=[scol[s] for s in frac.columns], labels=frac.columns)
        ax.set(title=f"{t}: state share over time", xlabel="time (h)", ylim=(0, 1)); ax.legend(fontsize=6, loc="upper left")
    # behaviour and niche per state
    feats = ["rule_division_rate", "rule_apoptosis_rate", "rule_speed", "signal_pressure", "nb_total"] + \
            [f"nb_{t}" for t in TYPES] + [c for c in p.obs if c.startswith("module_")]
    feats = [f for f in feats if f in p.obs]
    M = p.obs.groupby("state", observed=True)[feats].mean()
    Z = (M - M.mean()) / M.std().replace(0, 1)
    ax = fig.add_subplot(gs[1, :2])
    im = ax.imshow(Z.T.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(feats))); ax.set_yticklabels([f.replace("rule_", "").replace("signal_", "") for f in feats], fontsize=8)
    ax.set_xticks(range(len(M))); ax.set_xticklabels(M.index, rotation=60, fontsize=7)
    ax.set_title("behaviour, niche (nb_* = share of contacts) and modules per state (z across states)")
    plt.colorbar(im, ax=ax, shrink=0.7)
    # marker genes per state, within each type
    ax = fig.add_subplot(gs[1, 2:])
    rows = []
    for t in TYPES:
        b = p[p.obs.cell_type == t].copy()
        if b.obs.state.nunique() < 2:
            continue
        b.obs["state"] = b.obs.state.astype(str).astype("category")
        sc.tl.rank_genes_groups(b, "state", method="wilcoxon")
        for s in b.obs.state.cat.categories:
            rows += [(s, g) for g in sc.get.rank_genes_groups_df(b, s).names[:3]]
    genes = list(dict.fromkeys(g for _, g in rows))
    E = pd.DataFrame(p[:, genes].X.toarray() if hasattr(p.X, "toarray") else p[:, genes].X, columns=genes, index=p.obs_names)
    E = E.groupby(p.obs.state, observed=True).mean()
    E = (E - E.mean()) / E.std().replace(0, 1)
    im = ax.imshow(E.T.values, aspect="auto", cmap="viridis")
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6)
    ax.set_xticks(range(len(E))); ax.set_xticklabels(E.index, rotation=60, fontsize=7)
    ax.set_title("top 3 marker genes per state vs other states of the same type (z of mean log expression)")
    plt.colorbar(im, ax=ax, shrink=0.7)
    fig.tight_layout(); save(fig, "09a_emergent_states")
    # states in space at the last snapshot
    last = p[p.obs.time_h == TIMES[-1]]
    fig, axx = plt.subplots(1, len(TYPES), figsize=(5 * len(TYPES), 5))
    lim = np.abs(A.obs[["x", "y"]].to_numpy()).max() * 1.02
    allo = A.obs[A.obs.time_h == TIMES[-1]]
    for k, t in enumerate(TYPES):
        axx[k].scatter(allo.x, allo.y, s=1, color="0.88")
        for s in [s for s in states if s.startswith(t + "_s")]:
            m = (last.obs.state == s).to_numpy()
            axx[k].scatter(last.obs.x[m], last.obs.y[m], s=4, color=scol[s], label=s)
        axx[k].set(title=f"{t} states in space, {TIMES[-1]:g} h (sampled cells)", xlim=(-lim, lim), ylim=(-lim, lim))
        axx[k].set_aspect("equal"); axx[k].set_xticks([]); axx[k].set_yticks([]); axx[k].legend(fontsize=6, markerscale=3)
    save(fig, "09b_states_in_space")


@plot
def _10_events():
    ev = pd.read_csv(f"{RUN}/events.tsv", sep="\t")
    tr = np.load(f"{RUN}/tracks.npz")
    pos = pd.DataFrame({"time_h": tr["time_h"], "cell_id": tr["cell_id"], "x": tr["x"], "y": tr["y"]})
    ev = ev[ev.kind.isin(["division", "death_start"])].copy()
    pos["tk"] = pos.time_h.round(3)
    ts = np.sort(pos.tk.unique())
    ev["tk"] = ts[np.clip(np.searchsorted(ts, ev.time_h.round(3)), 0, len(ts) - 1)]
    ev = ev.merge(pos[["tk", "cell_id", "x", "y"]], on=["tk", "cell_id"], how="inner")
    lim = np.abs(A.obs[["x", "y"]].to_numpy()).max() * 1.02
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    last = A.obs[A.obs.time_h == TIMES[-1]]
    for j, (kind, cm) in enumerate((("division", "Greens"), ("death_start", "Reds"))):
        e = ev[ev.kind == kind]
        ax[j].scatter(last.x, last.y, s=1, color="0.9")
        im = ax[j].scatter(e.x, e.y, s=3, c=e.time_h, cmap=cm, lw=0)
        plt.colorbar(im, ax=ax[j], label="time (h)", shrink=0.8)
        ax[j].set(title=f"{kind.replace('_start', '')} locations ({len(e)})", xlim=(-lim, lim), ylim=(-lim, lim))
        ax[j].set_aspect("equal"); ax[j].set_xticks([]); ax[j].set_yticks([])
    bins = np.linspace(0, TIMES[-1], 25)
    for kind, ls in (("division", "-"), ("death_start", "--")):
        for t in TYPES:
            e = ev[(ev.kind == kind) & (ev.type == t)]
            h, _ = np.histogram(e.time_h, bins)
            ax[2].plot(bins[1:], h, ls, color=COL[t], label=f"{t} {kind.replace('_start', '')}")
    ax[2].set(title="events per time bin", xlabel="time (h)"); ax[2].legend(fontsize=6)
    save(fig, "10_division_death_events")


@plot
def _11_pressure():
    if "signal_pressure" not in A.obs:
        print("  (no pressure rule in this run: skipped)"); return
    ts = pick_times(3)
    lim = np.abs(A.obs[["x", "y"]].to_numpy()).max() * 1.02
    fig, ax = plt.subplots(1, len(ts) + 1, figsize=(4.4 * (len(ts) + 1), 4.2))
    for j, t in enumerate(ts):
        ob = snap(t).obs
        im = ax[j].scatter(ob.x, ob.y, s=2, c=ob.signal_pressure, cmap="inferno", vmin=0, vmax=3, lw=0)
        ax[j].set(title=f"pressure, {t:g} h", xlim=(-lim, lim), ylim=(-lim, lim)); ax[j].set_aspect("equal")
        ax[j].set_xticks([]); ax[j].set_yticks([])
    plt.colorbar(im, ax=ax[len(ts) - 1], shrink=0.8, label="pressure (1 = start median)")
    o = A.obs.sample(min(20000, A.n_obs), random_state=0)
    for t in TYPES:
        m = o.cell_type == t
        ax[-1].scatter(o.signal_pressure[m], 60 * o.rule_division_rate[m], s=1, color=COL[t], alpha=0.4, label=t)
    ax[-1].set(xlabel="pressure", ylabel="division rate (per h)", title="contact inhibition", xlim=(0, 4)); ax[-1].legend(fontsize=7, markerscale=5)
    save(fig, "11_pressure")


@plot
def _12_svg():
    I = np.array([morans_i(snap(t).obs, np.log2(np.asarray(snap(t).layers["S_true"]) + 1)) for t in TIMES])
    top = np.argsort(-I.mean(0))[:25]
    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(I[:, top].T, aspect="auto", cmap="viridis")
    ax.set_yticks(range(len(top))); ax.set_yticklabels(A.var_names[top], fontsize=7)
    ax.set_xticks(range(0, len(TIMES), max(1, len(TIMES) // 6))); ax.set_xticklabels([f"{x:g}" for x in TIMES[::max(1, len(TIMES) // 6)]])
    ax.set(title="spatially variable genes: Moran's I (7-NN) over time, top 25", xlabel="time (h)")
    plt.colorbar(im, ax=ax, label="Moran's I")
    save(fig, "12_spatially_variable_genes")


@plot
def _13_drivers():
    last = snap(TIMES[-1])
    logS = np.log2(np.asarray(last.layers["S_true"]) + 1)
    rules = yaml.safe_load(open(f"{RUN}/config_rules.yaml")) if os.path.exists(f"{RUN}/config_rules.yaml") else {"rules": []}
    mods = {}
    for g, m in zip(A.var_names, A.var.modules):
        for it in filter(None, str(m).split(";")):
            mods.setdefault(it.split(":")[0], set()).add(g)
    behs = ["rule_division_rate", "rule_apoptosis_rate", "rule_speed"]
    fig, ax = plt.subplots(len(TYPES), len(behs), figsize=(5 * len(behs), 3.6 * len(TYPES)))
    for k, t in enumerate(TYPES):
        m = (last.obs.cell_type == t).to_numpy()
        for j, b in enumerate(behs):
            y = last.obs[b].to_numpy()[m]
            if np.std(y) == 0:
                ax[k, j].set_title(f"{t}: {b[5:]} constant"); ax[k, j].axis("off"); continue
            rho = np.array([spearmanr(logS[m, g], y).correlation if np.std(logS[m, g]) > 0 else 0 for g in range(A.n_vars)])
            rho = np.nan_to_num(rho)
            top = np.argsort(-np.abs(rho))[:12][::-1]
            truth = set()
            for r in rules.get("rules") or []:
                ct = r.get("cell_type", "*")
                applies = ct == "*" or ct == t or (isinstance(ct, list) and t in ct)
                if r.get("enabled", True) and r.get("behaviour") == b[5:] and applies:
                    truth |= mods.get(r["signal"], {r["signal"]} if r["signal"] in A.var_names else set())
            cols = ["C3" if A.var_names[g] in truth else "0.6" for g in top]
            ax[k, j].barh(range(len(top)), rho[top], color=cols)
            ax[k, j].set_yticks(range(len(top))); ax[k, j].set_yticklabels(A.var_names[top], fontsize=7)
            ax[k, j].set(title=f"{t}: genes vs {b[5:]}", xlabel="Spearman")
    fig.suptitle("Genes most associated with each behaviour (red = genes the rules actually read)", fontsize=10)
    fig.tight_layout(); save(fig, "13_behaviour_drivers")


todo = list(PLOTS) if a.plots == "all" else a.plots.split(",")
for k in todo:
    print(k, flush=True)
    try:
        PLOTS[k]()
    except Exception as e:                                       # one failing plot should not stop the rest
        import traceback; traceback.print_exc()
        print(f"  FAILED {k}: {e}", flush=True)
print(f"plots in {PLOT}/")
