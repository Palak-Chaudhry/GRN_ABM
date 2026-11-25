"""GRN_ABM vs scMultiSim (dynamic GRN, spatial CCI) on the same inputs (comparison/run_scmultisim.R).

Both simulators get the same GRN, L-R pairs, L-R -> TF links, 3 cell types, 3000 cells and capture 0.2.
GRN_ABM also uses the reference's per-type expression (scMultiSim has no way to take it), so for realism
scMultiSim's populations are matched to reference types by best expression correlation.

  1. realism        vs the Atera reference, observed counts, shared genes; depth matched by thinning every
                    simulated cell to a library size drawn from the matched reference type:
                    per-type mean expression, zero fraction vs mean, dispersion, gene-gene correlation
                    structure, library sizes
  2. UMAP / space   embeddings of simulated and reference cells; spatial layouts
  3. GRN fidelity   within each type, |Spearman(regulator, target)| of true GRN edges vs random non-edges
                    (AUROC) and sign agreement, from true and observed counts
  4. CCI fidelity   neighbours' ligand vs the receiver's receptor and the receiver's downstream TFs,
                    within-type centred, z vs permuted ligand
  5. dynamics       per-cell memory over time; dynamic GRN: does the target follow its edge?
  6. velocity / cost / features

Usage: python comparison/compare_scmultisim.py [--ours results/runs/main_72h_v2]
           [--sms results/comparison/scmultisim_spatial] [--smv results/comparison/scmultisim_velocity]
           [--out results/comparison]
"""
import argparse, os, re, sys
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
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, rankdata
import scipy.sparse as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--ours", default="results/runs/main_72h_v2")
ap.add_argument("--sms", default="results/comparison/scmultisim_spatial")
ap.add_argument("--smv", default="results/comparison/scmultisim_velocity")
ap.add_argument("--out", default="results/comparison")
a = ap.parse_args()
P = lambda p: p if os.path.isabs(p) else f"{ROOT}/{p}"
OURS, SMS, SMV, OUT = P(a.ours), P(a.sms), P(a.smv), P(a.out)
PLOT = f"{OUT}/plots"; os.makedirs(PLOT, exist_ok=True)
rng = np.random.default_rng(0)
log = []


def say(*x):
    log.append(" ".join(map(str, x))); print(log[-1], flush=True)


def save(fig, name):
    fig.savefig(f"{PLOT}/{name}.png", dpi=130, bbox_inches="tight"); plt.close(fig); print("  wrote", name)


# ====================================================================== load
A = ad.read_h5ad(f"{OURS}/snapshots.h5ad")
TYPES = list(A.uns["cell_types"])
COL = {t: plt.get_cmap("tab10")(i) for i, t in enumerate(TYPES)}
ECFG = yaml.safe_load(open(f"{OURS}/config_engine.yaml"))
ABM = yaml.safe_load(open(f"{OURS}/config_abm.yaml"))
T0, TEND = A.obs.time_h.min(), A.obs.time_h.max()


def ours_at(t):
    m = ((A.obs.time_h == t) & ~A.obs.dying).to_numpy()
    b = A[m]
    return dict(obs=np.asarray(b.layers["spliced"].toarray() if sp.issparse(b.layers["spliced"]) else b.layers["spliced"]),
                true=np.asarray(b.layers["S_true"]), type=b.obs.cell_type.astype(str).to_numpy(),
                xy=b.obs[["x", "y"]].to_numpy(), genes=list(A.var_names), cell_id=b.obs.cell_id.to_numpy())


def read_sms(d):
    rd = lambda f: pd.read_csv(f"{d}/{f}", index_col=0)
    tc, oc = rd("true_counts.csv.gz"), rd("observed_counts.csv.gz")
    meta = pd.read_csv(f"{d}/cell_meta.csv")
    out = dict(true=tc.T.to_numpy(float), obs=oc.T.to_numpy(float), genes=list(tc.index),
               pop=meta["pop"].astype(str).to_numpy(), meta=meta)
    if {"x", "y"} <= set(meta):
        out["xy"] = meta[["x", "y"]].to_numpy(float)
    info = dict(l.split(": ") for l in open(f"{d}/run_info.txt").read().strip().splitlines())
    out["seconds"] = float(info["seconds"])
    return out


O0, OE = ours_at(T0), ours_at(TEND)
S = read_sms(SMS)
say(f"# GRN_ABM vs scMultiSim on the same inputs\n")
say(f"- GRN_ABM: {OURS} ({len(O0['type'])} cells at t = 0, {len(OE['type'])} at {TEND:g} h)")
say(f"- scMultiSim spatial: {SMS} ({len(S['pop'])} cells, populations {sorted(set(S['pop']))})")

# reference: full panel, sampled cells of the configured types
ref_cfg = ECFG["baseline"]["reference"]
lab = ad.read_h5ad(os.path.join(ref_cfg["dir"], ref_cfg["labels_h5ad"]), backed="r").obs[ref_cfg.get("label_key", "cell_type")].astype(str)
if ref_cfg.get("counts_h5"):
    R = sc.read_10x_h5(os.path.join(ref_cfg["dir"], ref_cfg["counts_h5"]), gex_only=True); R.var_names_make_unique()
else:
    R = ad.read_h5ad(os.path.join(ref_cfg["dir"], ref_cfg["labels_h5ad"]))
GENES = [g for g in A.var_names if g in set(S["genes"]) and g in set(R.var_names)]
say(f"- shared genes (network genes present in both simulators and the reference): {len(GENES)} of {A.n_vars}")
ref = {}
for t, label in ECFG["cell_types"].items():
    cells = lab.index[lab == label].intersection(R.obs_names)
    cells = rng.choice(cells, min(1000, len(cells)), replace=False)
    ref[t] = np.asarray(R[cells, GENES].X.toarray(), float)
del R


def sub(d, genes=GENES):
    idx = [d["genes"].index(g) for g in genes]
    return {k: (v[:, idx] if k in ("obs", "true") else v) for k, v in d.items()}


O0, OE, S = sub(O0), sub(OE), sub(S)

# match scMultiSim populations to reference types (it cannot be given type identities)
lcpm = lambda X: np.log1p(X.sum(0) / max(X.sum(), 1) * 1e4)
pops = sorted(set(S["pop"]))
C = np.array([[np.corrcoef(lcpm(S["obs"][S["pop"] == p]), lcpm(ref[t]))[0, 1] for t in TYPES] for p in pops])
r_, c_ = linear_sum_assignment(-C)
POP2T = {pops[i]: TYPES[j] for i, j in zip(r_, c_)}
S["type"] = np.array([POP2T.get(p, "unmatched") for p in S["pop"]])
say(f"- scMultiSim populations matched to reference types by expression: "
    + ", ".join(f"{p} -> {POP2T[p]} (r = {C[pops.index(p), TYPES.index(POP2T[p])]:.2f})" for p in pops if p in POP2T))

SIMS = {"GRN_ABM t=0": O0, f"GRN_ABM {TEND:g} h": OE, "scMultiSim": S}
SCOL = {"reference": "k", "GRN_ABM t=0": "C0", f"GRN_ABM {TEND:g} h": "C9", "scMultiSim": "C3"}


def rank_corr(X):
    """Spearman correlation matrix of the columns (NaN for constant columns), always a matrix."""
    Rk = np.apply_along_axis(rankdata, 0, X)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.corrcoef(Rk.T)


def thin(X, t):
    """Binomial thinning of each cell to a library size drawn from reference type t (never upward)."""
    lib = X.sum(1)
    target = rng.choice(ref[t].sum(1), len(X))
    p = np.clip(target / np.maximum(lib, 1), 0, 1)
    return rng.binomial(X.astype(np.int64), p[:, None]).astype(float)


# ====================================================================== 1. realism
say("\n## 1. Realism vs the reference (observed counts, shared genes, depth-matched)\n")
rows = []
fig, ax = plt.subplots(len(TYPES), 4, figsize=(20, 4.6 * len(TYPES)))
for k, t in enumerate(TYPES):
    Rt = ref[t]
    mr = np.log1p(Rt.mean(0) / Rt.sum(1).mean() * 1e3)
    zr, dr = (Rt == 0).mean(0), np.log10(Rt.var(0) / np.maximum(Rt.mean(0), 1e-9) + 1e-3)
    ax[k, 0].hist(np.log10(Rt.sum(1) + 1), bins=40, histtype="step", color="k", density=True, label="reference")
    ax[k, 1].scatter(np.log10(Rt.mean(0) + 1e-3), zr, s=6, color="k", label="reference")
    ax[k, 2].scatter(np.log10(Rt.mean(0) + 1e-3), dr, s=6, color="k", label="reference")
    top = np.argsort(-Rt.mean(0))[:40]
    cr = rank_corr(Rt[:, top])
    iu = np.triu_indices(len(top), 1)
    for name, d in SIMS.items():
        X = d["obs"][d["type"] == t]
        if len(X) < 20:
            continue
        Xt = thin(X, t)
        ms = np.log1p(X.mean(0) / X.sum(1).mean() * 1e3)
        zs, ds = (Xt == 0).mean(0), np.log10(Xt.var(0) / np.maximum(Xt.mean(0), 1e-9) + 1e-3)
        cs = rank_corr(Xt[:, top])
        ok = np.isfinite(cs[iu]) & np.isfinite(cr[iu])
        row = dict(sim=name, type=t, cells=len(X),
                   mean_r=np.corrcoef(ms, mr)[0, 1],
                   zero_frac_abs_err=np.mean(np.abs(zs - zr)),
                   dispersion_abs_err=np.nanmean(np.abs(ds - dr)),
                   genecorr_r=np.corrcoef(cs[iu][ok], cr[iu][ok])[0, 1] if ok.sum() > 10 else np.nan,
                   median_lib_raw=np.median(X.sum(1)), median_lib_ref=np.median(Rt.sum(1)))
        rows.append(row)
        ax[k, 0].hist(np.log10(X.sum(1) + 1), bins=40, histtype="step", color=SCOL[name], density=True, label=name)
        ax[k, 1].scatter(np.log10(Xt.mean(0) + 1e-3), zs, s=6, color=SCOL[name], alpha=0.6, label=name)
        ax[k, 2].scatter(np.log10(Xt.mean(0) + 1e-3), ds, s=6, color=SCOL[name], alpha=0.6, label=name)
        ax[k, 3].scatter(mr, ms, s=6, color=SCOL[name], alpha=0.6, label=f"{name}: r = {row['mean_r']:.2f}")
    ax[k, 0].set(title=f"{t}: library size (raw, before depth matching)", xlabel="log10 counts per cell")
    ax[k, 1].set(title=f"{t}: zero fraction vs mean (depth-matched)", xlabel="log10 mean", ylabel="fraction of zeros")
    ax[k, 2].set(title=f"{t}: dispersion (depth-matched)", xlabel="log10 mean", ylabel="log10 variance / mean")
    ax[k, 3].set(title=f"{t}: relative mean expression vs reference", xlabel="reference log1p(per-1k)", ylabel="simulated")
    for x in ax[k]:
        x.legend(fontsize=6)
save(fig, "c1_realism")
real = pd.DataFrame(rows)
real.to_csv(f"{OUT}/realism.tsv", sep="\t", index=False)
summ = real.groupby("sim")[["mean_r", "zero_frac_abs_err", "dispersion_abs_err", "genecorr_r"]].mean()
say("Mean over cell types (scMultiSim types = its populations matched to the reference):\n")
say("| simulator | mean expression r vs reference (↑) | zero-fraction error (↓) | dispersion error, log10 (↓) | gene-gene correlation structure r (↑) |\n|---|---|---|---|---|")
for s_, r in summ.iterrows():
    say(f"| {s_} | {r.mean_r:.2f} | {r.zero_frac_abs_err:.3f} | {r.dispersion_abs_err:.2f} | {r.genecorr_r:.2f} |")

# ====================================================================== 2. UMAP / space
say("\n## 2. Embeddings and space\n")
fig, ax = plt.subplots(1, 3, figsize=(17, 5))
for j, (name, X, ty) in enumerate((("reference (sampled)", np.vstack([ref[t] for t in TYPES]), np.repeat(TYPES, [len(ref[t]) for t in TYPES])),
                                   (f"GRN_ABM {TEND:g} h", OE["obs"], OE["type"]), ("scMultiSim", S["obs"], S["type"]))):
    b = ad.AnnData(X.astype(np.float32)); b.obs["t"] = ty
    sc.pp.normalize_total(b, target_sum=float(np.median(X.sum(1)) or 1)); sc.pp.log1p(b)
    sc.pp.pca(b, n_comps=30); sc.pp.neighbors(b); sc.tl.umap(b)
    for t in TYPES:
        m = b.obs.t.to_numpy() == t
        ax[j].scatter(*b.obsm["X_umap"][m].T, s=2, color=COL[t], label=t)
    ax[j].set(title=name, xticks=[], yticks=[])
ax[0].legend(markerscale=5, fontsize=8)
save(fig, "c2_umap")
fig, ax = plt.subplots(1, 2, figsize=(12, 5.5))
for t in TYPES:
    m = OE["type"] == t
    ax[0].scatter(*OE["xy"][m].T, s=2, color=COL[t], label=t)
    if "xy" in S:
        m = S["type"] == t
        ax[1].scatter(*S["xy"][m].T, s=6, marker="s", color=COL[t], label=t)
ax[0].set(title=f"GRN_ABM at {TEND:g} h (continuous space, µm; cells moved, divided, died)"); ax[0].set_aspect("equal")
ax[1].set(title="scMultiSim (grid; one cell per site, static)"); ax[1].set_aspect("equal"); ax[0].legend(markerscale=4)
save(fig, "c2_space")


def contacts(d, kind):
    xy = d["xy"]
    if kind == "grid":
        pr = cKDTree(xy).query_pairs(1.01, output_type="ndarray")
    else:
        rad = {t: float((ABM.get("types") or {}).get(t, {}).get("radius", ABM["defaults"]["radius"])) *
                  float((ABM.get("types") or {}).get(t, {}).get("max_adhesion_distance", ABM["defaults"]["max_adhesion_distance"])) for t in TYPES}
        r = np.array([rad[t] for t in d["type"]])
        pr = cKDTree(xy).query_pairs(2 * r.max(), output_type="ndarray")
        dd = np.linalg.norm(xy[pr[:, 0]] - xy[pr[:, 1]], axis=1)
        pr = pr[dd < r[pr[:, 0]] + r[pr[:, 1]]]
    n = len(xy)
    return sp.csr_matrix((np.ones(2 * len(pr)), (np.r_[pr[:, 0], pr[:, 1]], np.r_[pr[:, 1], pr[:, 0]])), shape=(n, n))


def enrichment(d, Aadj):
    ti = np.array([TYPES.index(t) if t in TYPES else -1 for t in d["type"]])
    i, j = Aadj.nonzero()
    def cm(lab):
        c = np.zeros((len(TYPES), len(TYPES)))
        ok = (lab[i] >= 0) & (lab[j] >= 0)
        np.add.at(c, (lab[i][ok], lab[j][ok]), 1)
        return c
    o = cm(ti); perm = np.array([cm(rng.permutation(ti)) for _ in range(30)])
    return (o - perm.mean(0)) / np.maximum(perm.std(0), 1e-9)


AO = contacts(OE, "abm")
AS = contacts(S, "grid") if "xy" in S else None
for name, d, Ad in ((f"GRN_ABM {TEND:g} h", OE, AO), ("scMultiSim", S, AS)):
    if Ad is None:
        continue
    z = enrichment(d, Ad)
    say(f"- {name}: neighbours per cell {Ad.sum(1).mean():.1f}; same-type contact enrichment z "
        + ", ".join(f"{t} {z[k, k]:+.0f}" for k, t in enumerate(TYPES)))

# ====================================================================== 3. GRN fidelity
say("\n## 3. GRN fidelity: are true edges visible in the data? (within each type)\n")
E = pd.read_csv(f"{ROOT}/{ECFG['network_dir']}/edges.tsv", sep="\t")
E = E[(E.regulator != E.target) & E.regulator.isin(GENES) & E.target.isin(GENES)]
gi = {g: k for k, g in enumerate(GENES)}
tfs = sorted(set(E.regulator))
edge_set = set(zip(E.regulator, E.target))
null = set()
while len(null) < 3 * len(E):
    r, t = rng.choice(tfs), rng.choice(GENES)
    if r != t and (r, t) not in edge_set:
        null.add((r, t))
null = list(null)


def spearman_pairs(X, pairs):
    Rk = np.apply_along_axis(rankdata, 0, X)
    Rk = (Rk - Rk.mean(0)) / np.maximum(Rk.std(0), 1e-12)
    return np.array([(Rk[:, gi[a]] * Rk[:, gi[b]]).mean() for a, b in pairs])


def auroc(pos, neg):
    s = np.r_[pos, neg]; y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    r = rankdata(s)
    return (r[y == 1].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


grn_rows = []
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
for name, d in (("GRN_ABM t=0", O0), (f"GRN_ABM {TEND:g} h", OE), ("scMultiSim", S)):
    for layer in ("true", "obs"):
        au, sg, pe, ne = [], [], [], []
        for t in TYPES:
            X = d[layer][d["type"] == t]
            if len(X) < 50:
                continue
            p = spearman_pairs(X, list(zip(E.regulator, E.target)))
            n = spearman_pairs(X, null)
            ok = np.isfinite(p); okn = np.isfinite(n)
            au.append(auroc(np.abs(p[ok]), np.abs(n[okn])))
            strong = ok & (np.abs(p) > 0.05)
            sg.append(np.mean(np.sign(p[strong]) == E.sign.to_numpy()[strong]) if strong.any() else np.nan)
            pe.append(np.abs(p[ok])); ne.append(np.abs(n[okn]))
        Xp = d[layer][np.isin(d["type"], TYPES)]                      # all types pooled
        pp, nn = spearman_pairs(Xp, list(zip(E.regulator, E.target))), spearman_pairs(Xp, null)
        pooled_auc = auroc(np.abs(pp[np.isfinite(pp)]), np.abs(nn[np.isfinite(nn)]))
        grn_rows.append(dict(sim=name, counts=layer, auroc=np.mean(au), auroc_pooled=pooled_auc, sign_agree=np.nanmean(sg),
                             median_abs_rho_edges=np.median(np.concatenate(pe)), median_abs_rho_null=np.median(np.concatenate(ne))))
        if layer == "true":
            ax[0].hist(np.concatenate(pe), bins=50, histtype="step", density=True, color=SCOL[name], label=f"{name} edges")
            ax[0].hist(np.concatenate(ne), bins=50, histtype="step", ls="--", density=True, color=SCOL[name], label=f"{name} non-edges")
g = pd.DataFrame(grn_rows)
for j, layer in enumerate(("true", "obs")):
    gg = g[g.counts == layer]
    ax[1].bar(np.arange(len(gg)) + 0.4 * j, gg.auroc, width=0.4, label=f"{layer} counts")
ax[1].set_xticks(np.arange(len(gg)) + 0.2); ax[1].set_xticklabels(gg.sim, fontsize=8)
ax[1].axhline(0.5, color="k", lw=0.5); ax[1].set(title="AUROC: |Spearman| separates true edges from non-edges", ylim=(0.4, 1)); ax[1].legend(fontsize=7)
ax[0].set(title="|Spearman(regulator, target)| within type, true counts", xlabel="|rho|"); ax[0].legend(fontsize=6)
save(fig, "c3_grn_fidelity")
g.to_csv(f"{OUT}/grn_fidelity.tsv", sep="\t", index=False)
say("| simulator | counts | AUROC edges vs non-edges, within type (0.5 = invisible) | AUROC, types pooled | sign agreement of correlated edges | median abs rho, edges / non-edges |\n|---|---|---|---|---|---|")
for r in g.itertuples():
    say(f"| {r.sim} | {r.counts} | {r.auroc:.3f} | {r.auroc_pooled:.3f} | {r.sign_agree:.2f} | {r.median_abs_rho_edges:.3f} / {r.median_abs_rho_null:.3f} |")

# ====================================================================== 4. CCI fidelity
say("\n## 4. Cell-cell signalling: does a neighbour's ligand show up in the receiver?\n")
LRd = pd.read_csv(f"{ROOT}/{ECFG['network_dir']}/lr_downstream.tsv", sep="\t")
subu = lambda x: [s for s in str(x).split("_") if s in gi]


def cci_specificity(d, Aadj, readout):
    """Per L-R pair: corr(sum of neighbours' ligand, receiver's readout), both centred within type, for the
    matching (cognate) readout and for the readouts of every other pair (non-cognate), after regressing out own
    type and neighbour counts per type from both. The cognate pair should stand out from the non-cognate ones."""
    X = np.log1p(d["true"]); ty = d["type"]
    has = np.asarray(Aadj.sum(1)).ravel() > 0
    # covariates: own type, and number of neighbours of each type (neighbour sums of every ligand rise and
    # fall together with these); both sides are residualised on them
    oh = np.array([[t == u for u in TYPES] for t in ty], float)
    Zc = np.c_[oh, Aadj @ oh][has]
    resid = lambda M: M - Zc @ np.linalg.lstsq(Zc, M, rcond=None)[0]
    Xc = np.zeros_like(X); Xc[has] = resid(X[has])
    ligs, reads, names = [], [], []
    for (l, r), grp in LRd.groupby(["ligand", "receptor"]):
        ls, rs = subu(l), subu(r)
        tg = [gi[x] for x in grp.tf if x in gi]
        if not ls or not rs or not tg:
            continue
        nb = Aadj @ X[:, [gi[s_] for s_ in ls]].min(1)
        ligs.append(resid(nb[has][:, None])[:, 0])
        reads.append((Xc[:, [gi[s_] for s_ in rs]].mean(1) if readout == "receptor" else Xc[:, tg].mean(1))[has])
        names.append(f"{l}->{r}")
    Lm, Ym = np.array(ligs), np.array(reads)
    Lz = (Lm - Lm.mean(1, keepdims=True)) / np.maximum(Lm.std(1, keepdims=True), 1e-12)
    Yz = (Ym - Ym.mean(1, keepdims=True)) / np.maximum(Ym.std(1, keepdims=True), 1e-12)
    Cm = Lz @ Yz.T / Lz.shape[1]                                     # ligand of pair i vs readout of pair j
    cog = np.diag(Cm)
    non = Cm[~np.eye(len(Cm), dtype=bool)]
    rank_pct = np.array([(Cm[i, np.arange(len(Cm)) != i] < Cm[i, i]).mean() for i in range(len(Cm))])
    return dict(cognate=cog, noncognate=non, auroc=auroc(cog, non), rank_pct=rank_pct, names=names)


fig, ax = plt.subplots(1, 2, figsize=(12, 4))
cci_rows = []
for j, readout in enumerate(("receptor", "downstream TFs")):
    for name, d, Ad in ((f"GRN_ABM {TEND:g} h", OE, AO), ("scMultiSim", S, AS)):
        if Ad is None:
            continue
        r = cci_specificity(d, Ad, readout)
        cci_rows.append(dict(sim=name, readout=readout, auroc=r["auroc"], median_cognate=np.median(r["cognate"]),
                             median_noncognate=np.median(r["noncognate"]), top10pct=np.mean(r["rank_pct"] >= 0.9)))
        ax[j].hist(r["noncognate"], bins=40, histtype="step", ls="--", density=True, color=SCOL[name], label=f"{name} non-cognate")
        ax[j].hist(r["cognate"], bins=20, histtype="step", density=True, color=SCOL[name], lw=1.5, label=f"{name} cognate")
    ax[j].set(title=f"neighbour ligand vs receiver {readout}", xlabel="partial correlation (type, neighbour composition removed)"); ax[j].legend(fontsize=6)
save(fig, "c4_cci_fidelity")
cc = pd.DataFrame(cci_rows); cc.to_csv(f"{OUT}/cci_fidelity.tsv", sep="\t", index=False)
say("| simulator | receiver readout | AUROC cognate vs non-cognate (0.5 = no specificity) | median corr cognate / non-cognate | cognate pair in its ligand's top 10% |\n|---|---|---|---|---|")
for r in cc.itertuples():
    say(f"| {r.sim} | {r.readout} | {r.auroc:.2f} | {r.median_cognate:+.3f} / {r.median_noncognate:+.3f} | {r.top10pct:.0%} |")
say("(GRN_ABM: ligand acts on the receiver's downstream TFs, gated by its receptor; scMultiSim: ligand acts on "
    "the receptor gene itself, which then regulates its TFs through the added receptor -> TF edges.)")

# ====================================================================== 5. dynamics
say("\n## 5. Dynamics\n")
# per-cell memory: the same cell 3 h apart (GRN_ABM only; scMultiSim cells are independent samples)
ts = np.sort(A.obs.time_h.unique())
t1, t2 = ts[-2], ts[-1]
a1, a2 = A[(A.obs.time_h == t1).to_numpy()], A[(A.obs.time_h == t2).to_numpy()]
common = np.intersect1d(a1.obs.cell_id, a2.obs.cell_id)
i1 = pd.Series(np.arange(a1.n_obs), index=a1.obs.cell_id.to_numpy())[common].to_numpy()
i2 = pd.Series(np.arange(a2.n_obs), index=a2.obs.cell_id.to_numpy())[common].to_numpy()
X1, X2 = np.log1p(np.asarray(a1.layers["S_true"])[i1]), np.log1p(np.asarray(a2.layers["S_true"])[i2])
ty = a1.obs.cell_type.astype(str).to_numpy()[i1]
for t in set(ty):
    m = ty == t
    X1[m] -= X1[m].mean(0); X2[m] -= X2[m].mean(0)
mem = np.array([np.corrcoef(X1[:, k], X2[:, k])[0, 1] for k in range(X1.shape[1]) if X1[:, k].std() > 0])
hl = A.var.halflife_h.to_numpy()[[k for k in range(X1.shape[1]) if X1[:, k].std() > 0]]
say(f"- GRN_ABM: the same cell {t2 - t1:g} h apart: median per-gene correlation (within type) {np.median(mem):.2f} "
    f"(genes with half-life > 10 h: {np.median(mem[hl > 10]):.2f}; < 3 h: {np.median(mem[hl < 3]):.2f}); "
    f"{len(common)} cells followed")
say("- scMultiSim: no cell is followed in time (each simulated cell is an independent sample; the probe found "
    "lag-1 within-cell autocorrelation ~0), so this cannot be measured")
# dynamic GRN: does the target follow its edge?
fig, ax = plt.subplots(1, 3, figsize=(19, 4.8))
f = f"{SMS}/dynamic_grn_changing_edges.csv.gz"
sms_track, sms_act = [], []
if os.path.exists(f):
    Wd = pd.read_csv(f, index_col=0)
    order = np.arange(Wd.shape[1])                                  # simulation order of cells
    full = pd.read_csv(f"{SMS}/true_counts.csv.gz", index_col=0)
    Xall = np.log1p(full.T.to_numpy(float)); gl = list(full.index)
    pop = S["pop"]
    im = ax[0].imshow(Wd.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-np.abs(Wd.to_numpy()).max(), vmax=np.abs(Wd.to_numpy()).max())
    ax[0].set_yticks(range(len(Wd))); ax[0].set_yticklabels(Wd.index, fontsize=6)
    ax[0].set(title=f"scMultiSim dynamic GRN: weights of the {len(Wd)} changing edges", xlabel="cell (simulation order)")
    plt.colorbar(im, ax=ax[0], shrink=0.7)
    for e, w in Wd.iterrows():
        reg, tgt = e.split("->")
        if tgt not in gl:
            continue
        y = Xall[:, gl.index(tgt)]
        yc = y.copy()
        for p in set(pop):
            yc[pop == p] -= y[pop == p].mean()
        if np.std(w.to_numpy()) > 0:
            sgn = np.sign(w.to_numpy()[np.argmax(np.abs(w.to_numpy()))])
            sms_track.append(spearmanr(w.to_numpy(), yc).correlation * sgn)
            if reg in gl:                                             # same quantity as GRN_ABM: weight x regulator
                act_c = w.to_numpy() * Xall[:, gl.index(reg)]
                sms_act.append(spearmanr(act_c, yc).correlation * sgn if np.std(act_c) > 0 else np.nan)
    say(f"- scMultiSim dynamic GRN: {len(Wd)} of the GRN's edges change (rewiring along cell order); target follows "
        f"its edge weight (sign-adjusted Spearman, within population): median {np.nanmedian(sms_track):+.2f}; "
        f"follows weight x regulator expression (the same quantity as GRN_ABM's edge activity): median "
        f"{np.nanmedian(sms_act) if sms_act else float('nan'):+.2f}")
g_ = np.load(f"{OURS}/grn_activity.npz")
W, act, tfs_a, genes_a, tt = g_["W"], g_["tf_activity"], list(g_["tfs"]), list(g_["genes"]), g_["time_h"]
e_tf, e_g = np.nonzero(W)
mean_t = {t: np.array([np.log1p(np.asarray(A.layers["S_true"])[((A.obs.time_h == x) & (A.obs.cell_type == t)).to_numpy()].mean(0))
                       for x in tt]) for t in TYPES}
ours_track = []
for k, t in enumerate(TYPES):
    Eact = act[:, k, e_tf] * W[e_tf, e_g]
    for j in np.argsort(-np.abs(Eact - Eact[:1]).max(0))[:30]:
        y = mean_t[t][:, e_g[j]]
        if np.std(Eact[:, j]) > 0 and np.std(y) > 0:
            ours_track.append(spearmanr(Eact[:, j], y).correlation * np.sign(W[e_tf[j], e_g[j]]))
Eact = act[:, 0, e_tf] * W[e_tf, e_g]
top = np.argsort(-np.abs(Eact - Eact[:1]).max(0))[:20]
im = ax[1].imshow((Eact - Eact[:1])[:, top].T, aspect="auto", cmap="PuOr_r")
ax[1].set_yticks(range(len(top))); ax[1].set_yticklabels([f"{tfs_a[e_tf[j]]}->{genes_a[e_g[j]]}" for j in top], fontsize=6)
ax[1].set_xticks(range(0, len(tt), 4)); ax[1].set_xticklabels([f"{x:g}" for x in tt[::4]])
ax[1].set(title=f"GRN_ABM ({TYPES[0]}): change in edge input over real time", xlabel="time (h)")
plt.colorbar(im, ax=ax[1], shrink=0.7)
say(f"- GRN_ABM: edges are fixed; their activity (regulator activation x weight) changes with each cell's state "
    f"and neighbourhood. Target follows its edge activity over time (per type, 30 most dynamic edges per type, "
    f"sign-adjusted Spearman): median {np.nanmedian(ours_track):+.2f}")
ax[2].hist(ours_track, bins=20, range=(-1, 1), histtype="step", color="C0", label="GRN_ABM (over time)")
if sms_track:
    ax[2].hist(sms_track, bins=20, range=(-1, 1), histtype="step", color="C3", label="scMultiSim: edge weight (over cell order)")
    ax[2].hist(sms_act, bins=20, range=(-1, 1), histtype="step", ls="--", color="C3", label="scMultiSim: weight x regulator")
ax[2].set(title="does the target follow its changing edge?", xlabel="sign-adjusted Spearman"); ax[2].legend(fontsize=7)
save(fig, "c5_dynamics")

# ====================================================================== 6. velocity, cost, features
say("\n## 6. Velocity, cost and features\n")
if os.path.exists(f"{SMV}/velocity.csv.gz"):
    V = pd.read_csv(f"{SMV}/velocity.csv.gz", index_col=0)
    say(f"- scMultiSim velocity: available only without spatial CCI (velocity run: {V.shape[1]} cells, "
        f"{np.isnan(V.to_numpy()).mean():.0%} NA). With CCI it returns no unspliced counts or velocity.")
say(f"- GRN_ABM velocity: true velocity for every cell and snapshot, with CCI, movement and divisions "
    f"(see the run's 02* plots for scvelo vs truth)")
ours_sec = None
for fn in sorted(os.listdir(f"{ROOT}/logs")):                   # the 08 log that wrote this run
    if fn.startswith("08_run_"):
        txt = open(f"{ROOT}/logs/{fn}").read()
        m_ = re.search(r"ran [\d.]+ h in (\d+) s", txt)
        if m_ and f"{OURS.rstrip('/')}/" in txt:
            ours_sec = float(m_.group(1))
say(f"- runtime: GRN_ABM {ours_sec or float('nan'):.0f} s for 3000 -> {len(OE['type'])} cells over {TEND:g} h "
    f"(720 coupled steps, 25 recorded snapshots); scMultiSim {S['seconds']:.0f} s for {len(S['pop'])} cells "
    f"(one static snapshot)")
say("""
| feature | GRN_ABM | scMultiSim |
|---|---|---|
| same GRN / L-R / L-R -> TF inputs | yes | yes (receptor -> TF as GRN edges; homophilic pairs dropped) |
| cell types match a reference's expression | yes (calibrated per type) | no (random type profiles) |
| real time; the same cell followed over time | yes | no (cells are independent samples along pseudotime) |
| cells move, divide, die; tissue grows and reorganises | yes | no (static grid) |
| expression changes behaviour (rules) | yes | no |
| per-type GRN activity | yes (TF expression differs by type) | no (one GRN for all types) |
| GRN rewiring (edges created / deleted) | no (fixed edges, activity changes) | yes (dynamic.GRN) |
| spliced / unspliced and velocity with CCI | yes | no (only without CCI) |
| differentiation trajectories | no (fixed types) | yes (tree) |
| ATAC / chromatin | no | yes |
""")
open(f"{OUT}/comparison_report.md", "w").write("\n".join(log) + "\n")
say(f"wrote {OUT}/comparison_report.md and {PLOT}/")
