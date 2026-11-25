"""Tissue (ABM) self-checks. Internal validation: each check compares the code with the model's own
equations, not with biology. Works for any number of cell types.

  1. two-cell equilibrium   final distance vs the analytic balance of repulsion and adhesion
  2. motility               mean squared displacement vs persistent-random-walk theory
  3. events                 division / death counts vs Poisson expectation
  4. packed tissue          relaxes without blow-up; overlap and neighbour counts; speed
  5. reference crop         (if the engine baseline has reference paths) cells and types loaded
  6. coupled run            tissue + gene engine for 24 h: rows in sync, growth vs rates, speed

Usage: python scripts/06_abm_checks.py [--engine configs/engine.yaml] [--abm configs/abm.yaml]
Output: results/abm_checks/<engine config>__<abm config>/{abm_checks.md, *.png}
"""
import argparse, copy, os, sys, time
import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{ROOT}/src")
from grn_abm import Tissue, Simulation

ap = argparse.ArgumentParser()
ap.add_argument("--engine", default=f"{ROOT}/configs/engine.yaml")
ap.add_argument("--abm", default=f"{ROOT}/configs/abm.yaml")
ap.add_argument("--hours", type=float, default=24)
a = ap.parse_args()
ECFG, ACFG = yaml.safe_load(open(a.engine)), yaml.safe_load(open(a.abm))
TYPES = list(ECFG["cell_types"])
DT = ACFG["dt_mech_min"]
stem = lambda x: os.path.splitext(os.path.basename(x))[0]
OUT = f"{ROOT}/results/abm_checks/{stem(a.engine)}__{stem(a.abm)}"
os.makedirs(OUT, exist_ok=True)
log = []


def say(*x):
    log.append(" ".join(map(str, x))); print(log[-1], flush=True)


def tissue(pos, type_idx, seed=0):
    t = Tissue(ACFG, TYPES, seed=seed)
    t._set_cells(np.asarray(pos, float), np.asarray(type_idx))
    return t


say(f"# ABM checks: {os.path.basename(a.abm)} ({len(TYPES)} types: {', '.join(TYPES)})")

# ---- 1. two-cell equilibrium per type
say("\n## 1. Two-cell equilibrium (no motility)\n\n| type | analytic distance (um) | simulated | error |\n|---|---|---|---|")
probe = Tissue(ACFG, TYPES)
for k, name in enumerate(TYPES):
    R, f = probe.p["radius"][k], probe.p["max_adhesion_distance"][k]
    cr, ca = probe.p["repulsion"][k], probe.p["adhesion"][k]
    rs, ra = 2 * R, 2 * f * R
    d_star = (np.sqrt(cr) - np.sqrt(ca)) / (np.sqrt(cr) / rs - np.sqrt(ca) / ra)
    t = tissue([[0, 0], [0.5 * R, 0]], [k, k])
    for _ in range(int(200 / DT)):
        t.mechanics_step(DT, motility=False)
    d = np.linalg.norm(t.pos[0] - t.pos[1])
    say(f"| {name} | {d_star:.3f} | {d:.3f} | {abs(d - d_star) / d_star:.2%} |")

# ---- 2. motility: MSD of isolated walkers vs theory 2 s^2 P (t - P (1 - exp(-t/P)))
say("\n## 2. Motility: mean squared displacement of isolated cells\n")
say("| type | speed (um/min) | persistence (min) | MSD ratio sim/theory at 10 / 60 / 300 min |\n|---|---|---|---|")
fig, ax = plt.subplots(figsize=(4.8, 3.4))
for k, name in enumerate(TYPES):
    s, P = probe.p["speed"][k], probe.p["persistence_min"][k]
    if s <= 0:
        say(f"| {name} | 0 | - | not motile |")
        continue
    n = 400
    g = np.stack(np.meshgrid(np.arange(20), np.arange(20)), -1).reshape(-1, 2) * 2000.0     # far apart
    t = tissue(g, np.full(n, k), seed=k)
    start, times, msd = t.pos.copy(), [], []
    for step in range(1, int(300 / DT) + 1):
        t.mechanics_step(DT)
        if step % int(1 / DT) == 0:
            times.append(step * DT)
            msd.append(np.mean(np.sum((t.pos - start) ** 2, 1)))
    times, msd = np.array(times), np.array(msd)
    theory = 2 * s ** 2 * P * (times - P * (1 - np.exp(-times / P)))
    r = [msd[np.searchsorted(times, x)] / theory[np.searchsorted(times, x)] for x in (10, 60, 300)]
    say(f"| {name} | {s:g} | {P:g} | {' / '.join(f'{x:.2f}' for x in r)} |")
    ax.loglog(times, msd, label=f"{name} sim")
    ax.loglog(times, theory, "k--", lw=0.8)
ax.set(xlabel="minutes", ylabel="MSD (um^2)", title="dashed = theory"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{OUT}/msd.png", dpi=120); plt.close(fig)

# ---- 3. event counts vs Poisson
say("\n## 3. Event counts (fixed population of 2000, 200 steps of 6 min)\n")
t = tissue(np.zeros((2000, 2)), np.zeros(2000, int), seed=3)
div_r, apo_r, steps, dt = 1e-3, 5e-4, 200, 6.0
nd = na = 0
for _ in range(steps):
    t.set_cell_params({"division_rate": np.full(2000, div_r), "apoptosis_rate": np.full(2000, apo_r)})
    p_, d_ = t.draw_events(dt)
    nd, na = nd + len(p_), na + len(d_)
p_apo = 1 - np.exp(-apo_r * dt)
p_div = (1 - np.exp(-div_r * dt)) * (1 - p_apo)          # a cell starting to die that step cannot divide
for lab, obs, p in (("divisions", nd, p_div), ("deaths", na, p_apo)):
    exp_, sd = 2000 * steps * p, np.sqrt(2000 * steps * p * (1 - p))
    say(f"- {lab}: {obs} observed vs {exp_:.0f} expected (z = {(obs - exp_) / sd:+.2f})")

# ---- 4. packed tissue
say("\n## 4. Packed tissue (synthetic layout, mechanics only then with motility)\n")
cfg4 = copy.deepcopy(ACFG)
t = Tissue(cfg4, TYPES, seed=4).init_synthetic()
n = len(t.pos)


def overlap_stats(t):
    R = t.par("radius")
    pr = t._pairs(2 * R.max())
    d = np.linalg.norm(t.pos[pr[:, 0]] - t.pos[pr[:, 1]], axis=1)
    rs = R[pr[:, 0]] + R[pr[:, 1]]
    return np.mean(d < 0.7 * rs), np.median(d[d < rs] / rs[d < rs]) if (d < rs).any() else np.nan


o0 = overlap_stats(t)
t0 = time.perf_counter()
for _ in range(int(60 / DT)):
    t.mechanics_step(DT, motility=False)
relax_ms = 1e3 * (time.perf_counter() - t0) / 10                      # per 6-min phenotype step
o1 = overlap_stats(t)
for _ in range(int(60 / DT)):
    t.mechanics_step(DT)
o2 = overlap_stats(t)
nb = t.adjacency().sum(1).A.ravel()
say(f"- {n} cells; strongly overlapping pairs (d < 0.7 x (R_i + R_j)): start {o0[0]:.1%} -> after 1 h relaxing "
    f"{o1[0]:.1%} -> after 1 h with motility {o2[0]:.1%}")
say(f"- median contact distance / (R_i + R_j): {o2[1]:.2f} (two-cell equilibrium for default cells "
    f"is {16.0 / 16.82:.2f})")
say(f"- NaN positions: {int(np.isnan(t.pos).sum())}; signalling neighbours per cell: median {np.median(nb):.0f}, "
    f"mean {nb.mean():.1f}, isolated cells {np.mean(nb == 0):.1%}")
say(f"- mechanics cost: {relax_ms:.0f} ms per 6-min phenotype step ({int(6 / DT)} substeps)")
fig, ax = plt.subplots(figsize=(4.8, 4.8))
for k, name in enumerate(TYPES):
    m = t.type_idx == k
    ax.scatter(t.pos[m, 0], t.pos[m, 1], s=(t.p["radius"][k] / 4) ** 2, label=name, alpha=0.7)
ax.set_aspect("equal"); ax.legend(fontsize=7, markerscale=2); ax.set(title="packed synthetic tissue after 2 h")
fig.tight_layout(); fig.savefig(f"{OUT}/tissue.png", dpi=120); plt.close(fig)

# ---- 5. reference crop
if ECFG["baseline"].get("reference") and ACFG["init"].get("reference"):
    try:
        tr = Tissue(ACFG, TYPES).init_reference(ECFG)
        nbr = tr.adjacency().sum(1).A.ravel()
        say(f"\n## 5. Reference crop ({ACFG['init']['reference']['crop_um']} um)\n- cells per type: "
            + ", ".join(f"{n_} {c}" for n_, c in zip(TYPES, np.bincount(tr.type_idx, minlength=len(TYPES))))
            + f"; signalling neighbours per cell: median {np.median(nbr):.0f} (only the configured types are kept, "
              f"so other cell types leave gaps)")
    except Exception as e:                                     # reference is optional
        say(f"\n## 5. Reference crop\n- skipped: {e}")

# ---- 6. coupled run with the gene engine
say(f"\n## 6. Coupled run: tissue + gene engine, {a.hours:g} h, base rates\n")
sim = Simulation(ECFG, ACFG, root=ROOT, seed=6)
c0 = sim.counts().astype(float)
times, counts, t0 = [], [], time.perf_counter()
for _ in range(int(a.hours * 60 / sim.dt)):
    sim.step()
    times.append(sim.t / 60); counts.append(sim.counts())
ms = 1e3 * (time.perf_counter() - t0) / len(times)
counts = np.array(counts)
net = sim.tissue.p["division_rate"] - sim.tissue.p["apoptosis_rate"]
say("| type | start | end (living) | expected from rates* |\n|---|---|---|---|")
for k, name in enumerate(TYPES):
    say(f"| {name} | {c0[k]:.0f} | {counts[-1, k]} | {c0[k] * np.exp(net[k] * a.hours * 60):.0f} |")
say("*N0 exp((division - death) t); dying cells are excluded from 'living', so 'end' sits slightly below.")
say(f"- rows in sync throughout (asserted every step); events logged: {len(sim.events)}; "
    f"cost {ms:.0f} ms per 6-min step -> 3 days of this tissue ~ {ms * 720 / 1e3 / 60:.1f} min")
fig, ax = plt.subplots(figsize=(4.8, 3.2))
for k, name in enumerate(TYPES):
    ax.plot(times, counts[:, k], label=name)
    ax.plot(times, c0[k] * np.exp(net[k] * np.array(times) * 60), "k--", lw=0.7)
ax.set(xlabel="hours", ylabel="living cells", title="dashed = expected from rates"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{OUT}/population.png", dpi=120); plt.close(fig)
open(f"{OUT}/abm_checks.md", "w").write("\n".join(log) + "\n")
say(f"\nwrote {OUT}/")
