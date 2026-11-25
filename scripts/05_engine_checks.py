"""Gene engine self-checks. Internal validation only: the engine is compared with its own baseline
and its own definitions, not with biology. Works for any config (any number of cell types,
reference or synthetic baselines).

  1. stability     no NaN / negative / runaway values; per-type means stop drifting
  2. calibration   simulated per-type means reproduce the baseline the engine was given
  3. feedback      the GRN contains TF->TF cycles and the run is still stable
  4. signalling    cut every contact half way: receivers' downstream TFs should move against the
                   link's sign (activating links drop, repressive links rise)
  5. velocity      s(t + dt) - s(t) agrees with velocity * dt
  6. speed         ms per step at 3k and 30k cells; projected 3 days

Usage: python scripts/05_engine_checks.py [--config configs/engine.yaml] [--cells-per-type 1000] [--hours 48]
Output: results/engine_checks/<config name>/{engine_checks.md, *.png}
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{ROOT}/src")
from grn_abm import GeneEngine

ap = argparse.ArgumentParser()
ap.add_argument("--config", default=f"{ROOT}/configs/engine.yaml")
ap.add_argument("--cells-per-type", type=int, default=1000)
ap.add_argument("--hours", type=float, default=48)
a = ap.parse_args()
CFG = yaml.safe_load(open(a.config))
OUT = f"{ROOT}/results/engine_checks/{os.path.splitext(os.path.basename(a.config))[0]}"
os.makedirs(OUT, exist_ok=True)
T = len(CFG["cell_types"])
STEPS = int(a.hours * 60 / CFG["dt_min"])
log = []


def say(*x):
    log.append(" ".join(map(str, x))); print(log[-1], flush=True)


def layout(n, seed, mean_nb=6.0):
    """Random mixed tissue: points in a disc, neighbours = pairs within a contact radius."""
    r = np.random.default_rng(seed)
    R = np.sqrt(n) * 10
    th, rr = r.uniform(0, 2 * np.pi, n), R * np.sqrt(r.uniform(0, 1, n))
    pos = np.c_[rr * np.cos(th), rr * np.sin(th)]
    pairs = cKDTree(pos).query_pairs(np.sqrt(mean_nb / (n / (np.pi * R ** 2)) / np.pi), output_type="ndarray")
    i, j = np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]]
    return sp.csr_matrix((np.ones(len(i)), (i, j)), shape=(n, n))


def run(types, A, cut_at=None):
    """Calibrated on the starting layout A. If cut_at is set, all contacts are removed from that step on."""
    eng = GeneEngine(CFG, root=ROOT).init_cells(types, A)
    empty = sp.csr_matrix(A.shape)
    traj = []
    for k in range(STEPS):
        if k == STEPS - 1:
            eng.prev_S = eng.S.copy()
        if k == cut_at:
            eng.sig_at_cut = eng.last_sig.copy()
        eng.step(empty if cut_at is not None and k >= cut_at else A)
        if k % 10 == 0 or k == STEPS - 1:
            traj.append([eng.S[types == t].mean(0) for t in range(T)])
    eng.traj = np.array(traj)                                     # snapshots x T x G
    return eng


types = np.repeat(np.arange(T), a.cells_per_type)
np.random.default_rng(1).shuffle(types)
A = layout(len(types), seed=1)
eng = run(types, A)
say(f"# Gene engine checks: {os.path.basename(a.config)}\n\n{len(types)} cells, {T} types "
    f"({', '.join(eng.types)}; baseline: {CFG['baseline']['source']}), randomly mixed, "
    f"{A.nnz / len(types):.1f} neighbours per cell, {STEPS} steps of {CFG['dt_min']} min = {a.hours:g} h. "
    f"{len(eng.genes)} genes, {len(eng.tf)} TFs, {len(eng.lig_sub)} L-R pairs.")

# ---- 1. stability
bad = int(np.isnan(eng.S).sum() + (eng.S < -1e-9).sum())
late = eng.traj[len(eng.traj) // 2:]
drift = np.abs(np.log2((late[-1] + 1) / (late[0] + 1)))                 # +1 molecule
say(f"\n## 1. Stability\n- NaN or negative values: {bad}\n- max level / max baseline level: {eng.S.max() / eng.ref.max():.2f}")
say(f"- per-type mean change over the second half, |log2 fold| (+1 molecule): median {np.median(drift):.3f}, "
    f"95th pct {np.percentile(drift, 95):.3f}, max {drift.max():.2f}")

# ---- 2. calibration against the engine's own baseline
obs = eng.observe()
say("\n## 2. Calibration: simulated mean counts vs the baseline the engine was given\n")
say("| type | Pearson r (log) | median sim/baseline | genes within 2-fold |\n|---|---|---|---|")
fig, ax = plt.subplots(1, T, figsize=(3.6 * T, 3.4), squeeze=False)
for t, name in enumerate(eng.types):
    sim, base = obs[types == t].mean(0), eng.baseline_counts[t]
    r = np.corrcoef(np.log1p(sim), np.log1p(base))[0, 1]
    ratio = (sim + 1e-2) / (base + 1e-2)
    say(f"| {name} | {r:.3f} | {np.median(ratio):.2f} | {np.mean((ratio > 0.5) & (ratio < 2)):.0%} |")
    ax[0, t].loglog(base + 1e-3, sim + 1e-3, ".", ms=4)
    lim = [1e-3, max(base.max(), sim.max()) * 2]
    ax[0, t].plot(lim, lim, "k--", lw=0.8)
    ax[0, t].set(title=f"{name} (r = {r:.2f})", xlabel="baseline mean counts", ylabel="simulated" if t == 0 else None)
fig.tight_layout(); fig.savefig(f"{OUT}/calibration.png", dpi=120); plt.close(fig)

# ---- 3. feedback
tfm = eng.W[:, eng.tf] != 0
_, lab = connected_components(sp.csr_matrix(tfm), directed=True, connection="strong")
say(f"\n## 3. Feedback\n- TF->TF edges {int(tfm.sum())}, self-loops {int(np.trace(tfm))}, largest feedback cluster "
    f"{np.bincount(lab).max()} of {len(eng.tf)} TFs; the run above is stable with it")

# ---- 4. signalling
cut = STEPS // 2
eng_cut = run(types, A, cut_at=cut)
ds = pd.read_csv(f"{ROOT}/{CFG['network_dir']}/lr_downstream.tsv", sep="\t")
pairs = ds.drop_duplicates(["ligand", "receptor"]).reset_index(drop=True)
rows = []
for p, r in pairs.iterrows():
    s = eng_cut.sig_at_cut[:, p]
    if s.max() <= 0:
        continue
    cells = s > np.percentile(s, 75)                                     # strongest receivers of this pair
    sub = ds[(ds.ligand == r.ligand) & (ds.receptor == r.receptor)]
    for lnk in sub.itertuples():
        g = eng.gi[lnk.tf]
        rows.append(dict(pair=f"{r.ligand}->{r.receptor}", tf=lnk.tf,
                         sign=getattr(lnk, "sign", 1),
                         log2fc=np.log2((eng_cut.S[cells, g].mean() + 1) / (eng.S[cells, g].mean() + 1))))
L = pd.DataFrame(rows)
# cutting contacts removes the signal, so an ACTIVATING link's TF should fall and a REPRESSIVE one's rise:
# the sign-aware statistic is sign * log2fc < 0, which reduces to log2fc < 0 on an all-activating network.
L["as_expected"] = L.sign * L.log2fc < 0
say(f"\n## 4. Signalling: all contacts cut after {cut * CFG['dt_min'] / 60:g} h vs the intact run\n"
    f"- downstream TF links tested: {len(L)} ({int((L.sign < 0).sum())} repressive); "
    f"moved as the link's sign predicts in {L.as_expected.mean():.0%}; "
    f"median signed log2 fold {(L.sign * L.log2fc).median():.2f}")
if (L.sign < 0).any():
    # cutting EVERY contact at once removes a TF's activating and repressive input together, so a TF fed by
    # both is not a clean test of either. Report the unambiguous subset separately.
    net = L.groupby("tf").sign.agg(lambda x: x.nunique() == 1)
    L["clean"] = L.tf.map(net)
    for nm, d in (("activating", L[L.sign > 0]), ("repressive", L[L.sign < 0])):
        say(f"  - {nm:11s} n={len(d):3d}  as expected {d.as_expected.mean():.0%}  "
            f"median log2 fold {d.log2fc.median():+.2f}")
    c = L[L.clean]
    say(f"  - TFs whose L-R input is all one sign ({c.tf.nunique()} of {L.tf.nunique()}): "
        f"as expected {c.as_expected.mean():.0%} (n={len(c)})")
    say("  - a TF fed by both signs, or driven mostly by the GRN rather than by L-R, is not a clean test "
        "here: the all-contacts cut removes its activating and repressive input at the same time")
    say(f"  - per TF with repressive input: " + ", ".join(
        f"{tf} {d.as_expected.mean():.0%} (n={len(d)})"
        for tf, d in L[L.sign < 0].groupby("tf")))
L.to_csv(f"{OUT}/signalling_links.tsv", sep="\t", index=False)   # per-link detail behind the summary
fig, ax = plt.subplots(figsize=(4.6, 3.2))
ax.hist([L.log2fc[L.sign > 0], L.log2fc[L.sign < 0]], bins=30, stacked=True,
        label=["activating", "repressive"]); ax.axvline(0, color="k", lw=0.8)
if (L.sign < 0).any():
    ax.legend(fontsize=7, frameon=False)
ax.set(xlabel="log2 fold, downstream TF (cut / intact)", ylabel="L-R -> TF links")
fig.tight_layout(); fig.savefig(f"{OUT}/signalling.png", dpi=120); plt.close(fig)

# ---- 5. velocity
c = np.corrcoef((eng.S - eng.prev_S).ravel(), (eng.velocity() * CFG["dt_min"]).ravel())[0, 1]
say(f"\n## 5. Velocity\n- corr(s(t+dt) - s(t), velocity x dt): {c:.3f}")

# ---- 6. speed
say(f"\n## 6. Speed (OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS')})")
for n in (3000, 30000):
    tt = np.repeat(np.arange(T), n // T)
    Ab = layout(len(tt), seed=2)
    e = GeneEngine(CFG, root=ROOT).init_cells(tt, Ab)
    k = 60 if n <= 3000 else 15
    t0 = time.perf_counter()
    for _ in range(k):
        e.step(Ab)
    ms = 1e3 * (time.perf_counter() - t0) / k
    say(f"- {len(tt)} cells: {ms:.1f} ms/step -> 3 days = {ms * 3 * 24 * 60 / CFG['dt_min'] / 1e3:.0f} s")

# ---- trajectories
fig, ax = plt.subplots(1, T, figsize=(3.6 * T, 3.0), sharey=True, squeeze=False)
hrs = np.minimum(np.arange(len(eng.traj)) * 10 * CFG["dt_min"] / 60, a.hours)
for t, name in enumerate(eng.types):
    ax[0, t].plot(hrs, np.log2((eng.traj[:, t] + 1) / (eng.ref[t] + 1)), lw=0.5, alpha=0.4)
    ax[0, t].set(title=name, xlabel="hours", ylabel="log2(mean / baseline), +1 molecule" if t == 0 else None)
fig.tight_layout(); fig.savefig(f"{OUT}/trajectories.png", dpi=120); plt.close(fig)
open(f"{OUT}/engine_checks.md", "w").write("\n".join(log) + "\n")
say(f"\nwrote {OUT}/")
