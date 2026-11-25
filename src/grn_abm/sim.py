"""Couples the tissue (abm.py), the gene engine (engine.py) and the behaviour rules (rules.py).

Start-up: the layout relaxes and moves for `warmup_h`; the gene engine is calibrated on it (including each type's growth dilution)
and starts from the bursting model's steady-state distribution; genes then run a burn-in (`burn_in_h`,
split into `calibration_rounds`, means nudged back to baseline after each) with the tissue held still and
virtual divisions (mRNA halved at each type's division rate), so t = 0 includes cell-cycle variation.

Every phenotype step (engine dt, e.g. 6 min):
  1. rules read the current state -> per-cell division / apoptosis rates, speed, chemotaxis, gene inputs
  2. mechanics + motility substeps (abm dt_mech) at those speeds
  3. signalling contacts + rule gene inputs -> gene engine step (dying cells neither send nor receive)
  4. events: divisions and new deaths at the rule rates
  5. births appended and finished deaths removed in BOTH layers, so row i is the same cell in both

Without rules, every cell uses its type's parameters (abm.yaml). Rules return per-cell overrides of any
tissue parameter, chemotaxis and gene inputs; this file only hands them to the two layers, so new rules
or new kinds of rule (plugins, see rules.py) need no change here.

Use:  sim = load_simulation("configs/engine.yaml", "configs/abm.yaml", "configs/rules.yaml")
      sim.run(hours=72, every_min=60, callback=lambda s: ...)       # or sim.step() in your own loop
"""
import numpy as np
import yaml

from .abm import Tissue
from .engine import GeneEngine
from .rules import Rules


class Simulation:
    def __init__(self, engine_cfg, abm_cfg, rules_cfg=None, root=".", seed=0):
        self.engine_cfg, self.abm_cfg, self.root = engine_cfg, abm_cfg, root
        self.tissue = Tissue(abm_cfg, list(engine_cfg["cell_types"]), seed=seed)
        if abm_cfg["init"]["source"] == "reference":
            self.tissue.init_reference(engine_cfg)
        else:
            self.tissue.init_synthetic()
        self.engine = GeneEngine(engine_cfg, root=root)
        # gene noise depends on both the engine seed and the run seed, so replicate runs differ in genes too
        self.engine.rng = np.random.default_rng((int(engine_cfg["seed"]), int(seed)))
        self.dt = engine_cfg["dt_min"]
        self.n_sub = max(1, int(round(self.dt / abm_cfg["dt_mech_min"])))
        # let the layout relax, then move for warmup_h (no divisions or deaths) so the packing at t = 0 is
        # that of a motile tissue; otherwise pressure jumps when motility starts and contact inhibition fires
        for _ in range(self.n_sub * 10):
            self.tissue.mechanics_step(abm_cfg["dt_mech_min"], motility=False)
        for _ in range(int(round(abm_cfg.get("warmup_h", 2) * 60 / abm_cfg["dt_mech_min"]))):
            self.tissue.mechanics_step(abm_cfg["dt_mech_min"])
        A = self.tissue.adjacency()
        growth = self.tissue.p["division_rate"]                       # per type, for dilution of mRNA
        self.engine.init_cells(self.tissue.type_idx, A, growth=growth)
        # gene burn-in on the still tissue, in rounds; after each round the engine is nudged back to its
        # baseline means. Cells halve their mRNA at their type's division rate (virtual divisions, cell count
        # unchanged), so t = 0 already contains the cell-cycle variation that real divisions create later.
        rounds = max(1, int(abm_cfg.get("calibration_rounds", 4)))
        per_round = int(abm_cfg.get("burn_in_h", 48) * 60 / self.dt / rounds)
        p_div = 1 - np.exp(-growth[self.tissue.type_idx] * self.dt)
        rng = np.random.default_rng(seed + 1)
        for k in range(rounds):
            for _ in range(per_round):
                self.engine.step(A)
                self.engine.halve(np.where(rng.random(len(p_div)) < p_div)[0])
            if k < rounds - 1:
                self.engine.recalibrate()
        self.t = 0.0
        self.events = []                                   # (time, kind, cell_id, type)
        self.last = {}
        # rules see the burnt-in state (their gene / module signals are scaled to it)
        self.rules = Rules(rules_cfg, self) if rules_cfg is not None else None

    def step(self, rules=None):
        """One phenotype step. rules: a Rules object; default = the one from rules_cfg (if any)."""
        tis, eng = self.tissue, self.engine
        rules = rules if rules is not None else self.rules
        r = rules(self) if rules is not None else {}
        type_at_rules = tis.type_idx.copy()
        tis.set_cell_params(r.get("params"))
        tis.set_chemotaxis(r.get("bias"), r.get("bias_dir"))
        for _ in range(self.n_sub):
            tis.mechanics_step(self.abm_cfg["dt_mech_min"])
        eng.step(tis.adjacency(), mod=r.get("gene_mod"))
        parents, dying = tis.draw_events(self.dt)
        self.last = dict(r, parents=parents, dying=dying, type_idx=type_at_rules)
        tis.start_dying(dying)
        tis.divide(parents)
        eng.divide(parents)
        gone = tis.advance_death(self.dt)
        self.events += [(self.t, "division", int(tis.cell_id[p]), int(tis.type_idx[p])) for p in parents]
        self.events += [(self.t, "death_start", int(tis.cell_id[d]), int(tis.type_idx[d])) for d in dying]
        self.events += [(self.t, "removed", int(tis.cell_id[g]), int(tis.type_idx[g])) for g in gone]
        tis.remove(gone)
        eng.remove(gone)
        self.t += self.dt
        assert len(tis.pos) == len(eng.S), "tissue and gene engine rows out of sync"

    def run(self, hours, rules=None, every_min=None, callback=None):
        """Advance `hours`; call callback(sim) every `every_min` minutes of simulated time (default: every step)."""
        every = max(1, int(round((every_min or self.dt) / self.dt)))
        for k in range(int(round(hours * 60 / self.dt))):
            self.step(rules)
            if callback is not None and (k + 1) % every == 0:
                callback(self)
        return self

    def counts(self):
        return np.bincount(self.tissue.type_idx[~self.tissue.dying], minlength=len(self.tissue.types))


def load_simulation(engine_yaml, abm_yaml, rules_yaml=None, root=".", seed=0):
    load = lambda p: yaml.safe_load(open(p if p.startswith("/") else f"{root}/{p}"))
    return Simulation(load(engine_yaml), load(abm_yaml), load(rules_yaml) if rules_yaml else None, root=root, seed=seed)
