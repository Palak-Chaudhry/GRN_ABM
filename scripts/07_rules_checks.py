"""Behaviour-rule self-checks (internal validation; any cell types, any rule set).

  1. formula        known signals give the values the PhysiCell rule form predicts
  2. forced module  raise a module that drives apoptosis_rate in half of one type: their death rate rises,
                    and death counts match the per-cell rates (Poisson)
  2b. chemotaxis    gradient estimate on a known field; drift of biased cells
  3. crowding       growth with the pressure rules vs without them
  4. gene rules     each gene:<GENE> rule raises GENE in the cells with a high signal (vs the rule removed)
  5. bad config     typos fail with readable errors
  6. full run       all enabled rules for 48 h: counts, behaviour values over time, per-rule response
                    (does each behaviour move with its signal in the right direction?), cost
  7. extensibility  a rule on a tissue parameter no default rule uses (adhesion) and a plugin's new signal
                    and behaviour (plugins/examples.py) work without code changes

Usage: python scripts/07_rules_checks.py [--engine ...] [--abm ...] [--rules configs/rules.yaml]
Output: results/rules_checks/<engine>__<abm>__<rules>/{rules_checks.md, *.png}
"""
import argparse, os, sys, time
import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{ROOT}/src")
from grn_abm import Simulation, Rules, Tissue, hill, combine

ap = argparse.ArgumentParser()
ap.add_argument("--engine", default=f"{ROOT}/configs/engine.yaml")
ap.add_argument("--abm", default=f"{ROOT}/configs/abm.yaml")
ap.add_argument("--rules", default=f"{ROOT}/configs/rules.yaml")
a = ap.parse_args()
ECFG, ACFG, RCFG = (yaml.safe_load(open(x)) for x in (a.engine, a.abm, a.rules))
stem = lambda x: os.path.splitext(os.path.basename(x))[0]
OUT = f"{ROOT}/results/rules_checks/{stem(a.engine)}__{stem(a.abm)}__{stem(a.rules)}"
os.makedirs(OUT, exist_ok=True)
TYPES = list(ECFG["cell_types"])
ENABLED = [r for r in RCFG.get("rules") or [] if r.get("enabled", True)]
log = []


def say(*x):
    log.append(" ".join(map(str, x))); print(log[-1], flush=True)


def new_sim(seed=0):
    return Simulation(ECFG, ACFG, root=ROOT, seed=seed)


def subset(keep):
    return {**RCFG, "rules": [r for r in ENABLED if keep(r)]}


def per_type(v, ti):
    return [v[ti == k].mean() if (ti == k).any() else np.nan for k in range(len(TYPES))]


say(f"# Behaviour-rule checks ({len(TYPES)} types: {', '.join(TYPES)}; {len(ENABLED)} enabled rules)")

# ---- 1. formula
base, mx, mn = np.array([1e-3]), np.array([5e-3]), np.array([0.0])
cases = [("no signal", combine(base, [(hill(np.array([0.0]), 1, 4), mx)], []), 1e-3),
         ("up rule at half_max", combine(base, [(hill(np.array([1.0]), 1, 4), mx)], []), 3e-3),
         ("up rule saturated", combine(base, [(hill(np.array([1e3]), 1, 4), mx)], []), 5e-3),
         ("down rule saturated", combine(base, [], [(hill(np.array([1e3]), 1, 4), mn)]), 0.0),
         ("up half + down half", combine(base, [(np.array([0.5]), mx)], [(np.array([0.5]), mn)]), 0.5 * 3e-3)]
say("\n## 1. Rule formula\n\n| case | computed | expected |\n|---|---|---|")
for lab, got, exp_ in cases:
    say(f"| {lab} | {got[0]:.3e} | {exp_:.3e} |")
assert all(abs(g[0] - e) < 1e-9 for _, g, e in cases), "rule formula mismatch"

# ---- 2. forced module -> apoptosis_rate
sim = new_sim(seed=2)
probe = Rules({"rules": []}, sim)
death_rules = [r for r in ENABLED if r["behaviour"] == "apoptosis_rate" and r["signal"] in probe.modules]
if death_rules:
    mod = death_rules[0]["signal"]
    rules = Rules(subset(lambda r: r["behaviour"] == "apoptosis_rate"), sim)
    eng, tis = sim.engine, sim.tissue
    cells = np.where(tis.type_idx == 0)[0]
    forced = np.random.default_rng(0).choice(cells, len(cells) // 2, replace=False)
    pro = rules.modules[mod]["pro"]
    eng.S[np.ix_(forced, pro)] = 4 * (eng.ref[0, pro] + 1) - 1                  # +2 log2 on every pro gene
    rates = rules(sim)["params"]["apoptosis_rate"]
    tis.set_cell_params({"apoptosis_rate": rates, "division_rate": np.zeros(len(rates))})
    reps, dt = 100, sim.dt
    obs = np.zeros(len(tis.pos))
    for _ in range(reps):
        obs[tis.draw_events(dt)[1]] += 1
    p = 1 - np.exp(-rates * dt)
    say(f"\n## 2. Forced {mod} module ({TYPES[0]}, +2 log2 on every pro gene in half the cells; {reps} "
        f"independent draws of one {dt:g}-min step)\n\n| group | mean rule rate (per min) | deaths | expected | z |\n|---|---|---|---|---|")
    for lab, idx in (("forced", forced), ("unforced (same type)", np.setdiff1d(cells, forced))):
        e_, v_ = reps * p[idx].sum(), reps * (p[idx] * (1 - p[idx])).sum()
        say(f"| {lab} | {rates[idx].mean():.2e} | {int(obs[idx].sum())} | {e_:.0f} | {(obs[idx].sum() - e_) / np.sqrt(v_):+.2f} |")

# ---- 2b. chemotaxis machinery
chemo_rules = [r for r in ENABLED if r["behaviour"].startswith("chemotaxis:")]
if chemo_rules:
    lig = chemo_rules[0]["behaviour"].split(":", 1)[1]
    s2 = new_sim(seed=21)
    r2 = Rules(RCFG, s2)
    x = s2.tissue.pos[:, 0]
    s2.engine.S[:, s2.engine.gi[lig]] = x - x.min()                     # linear field increasing along +x
    g = r2.gradient(lig, s2)
    say(f"\n## 2b. Chemotaxis\n- gradient estimate on a linear {lig} field along +x: mean cos(angle to +x) = "
        f"{g[:, 0].mean():.3f} (1 = perfect), radius {r2.chemo_radius:g} um")
    k = next((i for i in range(len(TYPES)) if Tissue(ACFG, TYPES).p["speed"][i] > 0), 0)
    grid = np.stack(np.meshgrid(np.arange(10), np.arange(10)), -1).reshape(-1, 2) * 2000.0
    tt = Tissue(ACFG, TYPES, seed=5)
    sp_ = max(tt.p["speed"][k], 1.0)
    for b in (1.0, 0.5, 0.0):
        tt._set_cells(grid, np.full(100, k)); x0 = tt.pos[:, 0].copy()
        tt.set_cell_params({"speed": np.full(100, sp_)})
        tt.set_chemotaxis(np.full(100, b), np.tile([1.0, 0.0], (100, 1)))
        for _ in range(600):
            tt.mechanics_step(0.1)
        say(f"- isolated cells, speed {sp_:g} um/min, bias {b:g}: drift along the bias direction "
            f"{(tt.pos[:, 0] - x0).mean() / 60:.2f} um/min")

# ---- 3. crowding
if any(r["signal"] == "pressure" for r in ENABLED):
    H3 = 36
    say(f"\n## 3. Pressure rules vs none, {H3} h\n")
    runs = {}
    for lab, rc in (("no rules", None), ("pressure rules only", subset(lambda r: r["signal"] == "pressure"))):
        s3 = new_sim(seed=3)
        r3 = Rules(rc, s3) if rc else None
        rr = r3 or Rules({"rules": []}, s3)                                 # for measuring pressure only
        tr, cnt, prs = [], [], []
        for _ in range(int(H3 * 60 / s3.dt)):
            s3.step(r3)
            tr.append(s3.t / 60); cnt.append(s3.counts().sum()); prs.append(np.median(rr.pressure(s3)))
        runs[lab] = (np.array(tr), np.array(cnt), np.array(prs))
        say(f"- {lab}: cells {cnt[0]} -> {cnt[-1]}; median relative pressure at end {prs[-1]:.2f}")
    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.2))
    for lab, (tr, cnt, prs) in runs.items():
        ax[0].plot(tr, cnt, label=lab); ax[1].plot(tr, prs, label=lab)
    ax[0].set(xlabel="hours", ylabel="living cells"); ax[1].set(xlabel="hours", ylabel="median pressure (1 = start)")
    ax[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(f"{OUT}/crowding.png", dpi=120); plt.close(fig)

# ---- 4. gene rules
gene_rules = [r for r in ENABLED if r["behaviour"].startswith("gene:")]
if gene_rules:
    say("\n## 4. Gene rules (12 h, all rules vs all rules except gene:<GENE>)\n")
    res = {}
    for lab, rc in (("with gene rules", subset(lambda r: True)),
                    ("without", subset(lambda r: not r["behaviour"].startswith("gene:")))):
        s4 = new_sim(seed=4)
        r4 = Rules(rc, s4)
        s4.run(12, r4)
        res[lab] = s4
    r_all = Rules(subset(lambda r: True), new_sim(seed=4))              # same start -> same signal scaling
    for r in gene_rules:
        g, sig = r["behaviour"][5:], r["signal"]
        line = []
        for lab, s4 in res.items():
            sc = r_all.signal(sig, s4)
            top = sc > np.percentile(sc, 90)
            e = s4.engine
            lf = np.log2((e.S[:, e.gi[g]] + 1) / (e.ref[e.type_idx, e.gi[g]] + 1))
            line.append(f"{lab}: {lf[top].mean():+.2f}")
        say(f"- {g} (log2 vs baseline) in the top 10% {sig} cells — " + "; ".join(line))

# ---- 5. bad config
try:
    Rules({"rules": [{"cell_type": "tumour", "signal": "apoptsis", "behaviour": "divide", "direction": "up",
                      "halfmax": 2}]}, new_sim())
    say("\n## 5. Bad config\n- NOT caught")
except ValueError as e:
    say("\n## 5. Bad config is rejected with:\n```\n" + str(e) + "\n```")

# ---- 6. full run
H6 = 48
sim = Simulation(ECFG, ACFG, RCFG, root=ROOT, seed=6)
rules = sim.rules
c0 = sim.counts()
SHOW = ["division_rate", "apoptosis_rate", "speed"]
tr, counts, vals, resp, align = [], [], {b: [] for b in SHOW}, {k: [] for k in range(len(rules.rules))}, []
t0 = time.perf_counter()
for _ in range(int(H6 * 60 / sim.dt)):
    before = dict(zip(sim.tissue.cell_id, sim.tissue.pos.copy()))
    sim.step()
    l = sim.last
    ti = l["type_idx"]
    for b in SHOW:                                                     # rule value, or the type's value
        v = l["values"].get(b)
        vals[b].append(per_type(v if v is not None else sim.tissue.p[b][ti], ti))
    sig = rules.last_signals
    for k, r in enumerate(rules.rules):                                # does each behaviour follow its signal?
        m = np.ones(len(ti), bool) if r["types"] is None else np.isin(ti, r["types"])
        s_, v_ = sig[r["signal"]][m], l["values"][r["behaviour"]][m]
        if v_.ndim == 1 and s_.std() > 0 and v_.std() > 0:
            resp[k].append(np.corrcoef(s_, v_)[0, 1])
    if "bias" in l:                                                    # do biased cells move along the gradient?
        ids_then = sorted(before)                                       # rule-time row order == cell_id order
        bd = l["bias_dir"] / np.maximum(np.linalg.norm(l["bias_dir"], axis=1), 1e-12)[:, None]
        now = dict(zip(sim.tissue.cell_id, sim.tissue.pos))
        cos = []
        for i in np.where(l["bias"] > 0.2)[0]:
            cid = ids_then[i]
            if cid in now:
                dsp = now[cid] - before[cid]
                cos.append(dsp @ bd[i] / max(np.linalg.norm(dsp), 1e-9))
        if cos:
            align.append((len(cos), float(np.mean(cos))))
    tr.append(sim.t / 60); counts.append(sim.counts())
ms = 1e3 * (time.perf_counter() - t0) / len(tr)
counts = np.array(counts)
bad = int(np.isnan(sim.engine.S).sum() + np.isnan(sim.tissue.pos).sum())
say(f"\n## 6. All enabled rules, {H6} h\n\n| type | cells start -> end | division rate (per h), first -> last | "
    f"death rate (per h), first -> last | speed (um/min), first -> last |\n|---|---|---|---|---|")
for k, name in enumerate(TYPES):
    d_, a_, s_ = (np.array(vals[b])[:, k] for b in SHOW)
    say(f"| {name} | {c0[k]} -> {counts[-1, k]} | {60 * d_[0]:.3f} -> {60 * d_[-1]:.3f} | "
        f"{60 * a_[0]:.4f} -> {60 * a_[-1]:.4f} | {s_[0]:.2f} -> {s_[-1]:.2f} |")
say("\nPer-rule response: mean over steps of corr(signal, behaviour value) in the cells the rule applies to. "
    "Expected sign: + for increases, - for decreases (other rules on the same behaviour dilute it).\n\n"
    "| rule | cells | corr | sign ok |\n|---|---|---|---|")
for k, r in enumerate(rules.rules):
    c = np.nanmean(resp[k]) if resp[k] else np.nan
    want = 1 if r["direction"] == "increases" else -1
    who = "all" if r["types"] is None else ",".join(TYPES[t] for t in r["types"])
    ok = "-" if np.isnan(c) else ("yes" if np.sign(c) == want else "NO")
    say(f"| {r['signal']} -> {r['behaviour']} ({r['direction']}) | {who} | {c:+.2f} | {ok} |")
if align:
    al = np.array(align)
    say(f"\n- chemotaxis: on average {al[:, 0].mean():.0f} cells per step with bias > 0.2; mean cos(angle between "
        f"their step displacement and the gradient) {al[:, 1].mean():+.2f} (0 = random walk)")
say(f"- NaN: {bad}; rows in sync (asserted every step); cost {ms:.0f} ms per step at the end-size tissue "
    f"(~{ms * 720 / 6e4:.1f} min per 3 days)")
fig, ax = plt.subplots(1, 4, figsize=(15, 3.2))
for k, name in enumerate(TYPES):
    ax[0].plot(tr, counts[:, k], label=name)
    for j, (b, sc) in enumerate(zip(SHOW, (60, 60, 1))):
        ax[j + 1].plot(tr, sc * np.array(vals[b])[:, k], label=name)
ax[0].set(xlabel="hours", ylabel="living cells"); ax[1].set(xlabel="hours", ylabel="division rate (per h)")
ax[2].set(xlabel="hours", ylabel="death rate (per h)"); ax[3].set(xlabel="hours", ylabel="speed (um/min)")
ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{OUT}/full_run.png", dpi=120); plt.close(fig)

# ---- 7. extensibility
say("\n## 7. Extensibility (no code changes)\n")
s7 = new_sim(seed=7)
other = TYPES[1 % len(TYPES)]
r7 = Rules({"rules": [{"cell_type": TYPES[0], "signal": f"contact_{other}", "behaviour": "adhesion",
                       "direction": "decreases", "min": 0.0, "half_max": 1.0}]}, s7)
s7.tissue.set_cell_params(r7(s7)["params"])
c7 = r7.last_signals[f"contact_{other}"]
touch, alone = (c7 > 0) & (s7.tissue.type_idx == 0), (c7 == 0) & (s7.tissue.type_idx == 0)
say(f"- tissue-parameter rule (contact_{other} -> adhesion decreases, {TYPES[0]} only): per-cell adhesion read "
    f"by the mechanics {s7.tissue.par('adhesion')[touch].mean():.3f} in {TYPES[0]} cells touching {other} vs "
    f"{s7.tissue.par('adhesion')[alone].mean():.3f} in the others (type value {s7.tissue.p['adhesion'][0]:g})")
mod = next(iter(r7.modules))
half = float(np.median(np.linalg.norm(s7.tissue.pos, axis=1)))
pcfg = {"plugins": ["plugins/examples.py"],
        "rules": [{"signal": "distance_from_centre", "behaviour": f"module:{mod}", "direction": "increases",
                   "max": 1.5, "half_max": half},
                  {"signal": "dying_neighbours", "behaviour": "speed", "direction": "decreases", "min": 0.0,
                   "half_max": 1.0}]}
res7 = {}
for lab, rc in (("with", pcfg), ("without", {"rules": []})):
    s = new_sim(seed=7)
    rr = Rules(rc, s)
    s.run(12, rr)
    res7[lab] = np.corrcoef(np.linalg.norm(s.tissue.pos, axis=1), rr._raw(mod, s, {}))[0, 1]
say(f"- plugin (plugins/examples.py): new signal distance_from_centre -> new behaviour module:{mod} for 12 h: "
    f"corr(distance from centre, {mod} score) {res7['with']:+.2f} with the rule vs {res7['without']:+.2f} "
    f"without; new signal dying_neighbours -> speed also loaded and ran")
open(f"{OUT}/rules_checks.md", "w").write("\n".join(log) + "\n")
say(f"\nwrote {OUT}/")
