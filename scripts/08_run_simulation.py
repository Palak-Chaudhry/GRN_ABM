"""Run a simulation and record it for analysis (see src/grn_abm/record.py for what is written).

Usage: python scripts/08_run_simulation.py [--engine configs/engine.yaml] [--abm configs/abm.yaml]
           [--rules configs/rules.yaml | --rules none] [--hours 72] [--snapshot-min 180] [--track-min 6]
           [--seed 0] [--name <run name>]
Output: results/runs/<name>/   (default name: <engine>__<abm>__<rules>__<hours>h_seed<seed>)
"""
import argparse, os, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{ROOT}/src")
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import pandas as pd
pd.options.future.infer_string = False
import yaml
from grn_abm import Simulation
from grn_abm.record import Recorder

ap = argparse.ArgumentParser()
ap.add_argument("--engine", default="configs/engine.yaml")
ap.add_argument("--abm", default="configs/abm.yaml")
ap.add_argument("--rules", default="configs/rules.yaml")
ap.add_argument("--hours", type=float, default=72)
ap.add_argument("--snapshot-min", type=float, default=180)
ap.add_argument("--track-min", type=float, default=6)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--name", default=None)
a = ap.parse_args()
path = lambda p: p if os.path.isabs(p) else f"{ROOT}/{p}"
stem = lambda p: os.path.splitext(os.path.basename(p))[0]
load = lambda p: yaml.safe_load(open(path(p)))
rules_cfg = None if a.rules == "none" else load(a.rules)
name = a.name or f"{stem(a.engine)}__{stem(a.abm)}__{stem(a.rules)}__{a.hours:g}h_seed{a.seed}"
out = f"{ROOT}/results/runs/{name}"

t0 = time.perf_counter()
sim = Simulation(load(a.engine), load(a.abm), rules_cfg, root=ROOT, seed=a.seed)
rec = Recorder(sim, snapshot_min=a.snapshot_min, track_min=a.track_min, seed=a.seed)
print(f"start-up {time.perf_counter() - t0:.0f} s; cells {sim.counts()}", flush=True)
sim.run(a.hours, callback=rec)
print(f"ran {a.hours:g} h in {time.perf_counter() - t0:.0f} s; cells {sim.counts()}", flush=True)
adata = rec.save(out)
for k, v in (("engine", a.engine), ("abm", a.abm), ("rules", a.rules)):
    if v != "none":
        os.system(f"cp {path(v)} {out}/config_{k}.yaml")
print(f"wrote {out}/ ({adata.n_obs} cell-snapshots, {len(rec.snaps)} snapshots)")
