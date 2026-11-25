"""Behaviour rules: any signal -> any behaviour, in the PhysiCell rule form. Configured in configs/rules.yaml.

A rule reads one SIGNAL per cell, passes it through a Hill function and pushes one BEHAVIOUR towards its
max (direction: increases) or min (decreases). Rules on the same behaviour combine (see `combine`).

Built-in signals
  <module>              gene-set score from network.yaml (mean log2 change of pro genes minus anti genes)
  <gene>                any network gene (log2 change vs its type baseline)
                        -> both in baseline SDs: (value - type mean) / type SD after the burn-in
  pressure              mechanical compression (1 = the starting tissue's median)
  neighbours            number of signalling contacts
  contact_<type>        contacts with cells of <type>
Built-in behaviours
  <tissue parameter>    any key of abm.yaml `defaults:` (division_rate, apoptosis_rate, speed, adhesion,
                        persistence_min, ...); base = the cell type's value
  gene:<GENE>           log change of GENE's transcription on-rate (base 0)
  chemotaxis:<LIGAND>   0-1 bias up the gradient of LIGAND expression around the cell (base 0)

Extending
  - more rules: add rows to rules.yaml (checked when loaded; a typo fails with a readable message)
  - new tissue parameter: add it to abm.yaml `defaults:`; rules can then set it per cell
  - new KIND of signal or behaviour: one function in a plugin file listed under `plugins:` in rules.yaml
    (see plugins/examples.py):

        from grn_abm.rules import register_signal, register_behaviour

        @register_signal("distance_from_centre")            # fn(rules, sim, arg, cache) -> per-cell values
        def _(rules, sim, arg, cache):
            return np.linalg.norm(sim.tissue.pos, axis=1)

        @register_behaviour("module:", check=...)           # fn(rules, sim, arg, value, out)
        def _(rules, sim, arg, value, out): ...              # write effects into `out` (see Rules.__call__)

    A registered name ending in "_" or ":" is a prefix; `arg` is the rest of the name (contact_tumor -> "tumor").
"""
import importlib.util
import os

import numpy as np
import pandas as pd


def hill(s, half, n):
    s = np.maximum(s, 0.0)
    return s ** n / (s ** n + half ** n)


def combine(base, up, down):
    """PhysiCell form. up / down: lists of (H, target) arrays. Returns per-cell values:
    (1 - D) [(1 - U) base + U max] + D min, with U / D = 1 - prod(1 - H) over the up / down rules and
    max / min the H-weighted mean of the rules' targets."""
    U, D = np.zeros_like(base), np.zeros_like(base)
    bmax, bmin = base.copy(), base.copy()
    if up:
        H = np.array([h for h, _ in up]); M = np.array([m for _, m in up])
        U = 1 - np.prod(1 - H, 0)
        w = H.sum(0)
        bmax = np.where(w > 0, (H * M).sum(0) / np.maximum(w, 1e-12), base)
    if down:
        H = np.array([h for h, _ in down]); M = np.array([m for _, m in down])
        D = 1 - np.prod(1 - H, 0)
        w = H.sum(0)
        bmin = np.where(w > 0, (H * M).sum(0) / np.maximum(w, 1e-12), base)
    return (1 - D) * ((1 - U) * base + U * bmax) + D * bmin


# ---------------------------------------------------------------------- registries
SIGNALS = {}        # name -> (fn(rules, sim, arg, cache) -> per-cell array, check(rules, sim, arg) -> bool)
BEHAVIOURS = {}     # name -> (fn(rules, sim, arg, value, out), check, base(rules, sim, arg) -> per-cell or None)


def register_signal(name, check=None):
    def deco(fn):
        SIGNALS[name] = (fn, check)
        return fn
    return deco


def register_behaviour(name, check=None, base=None):
    def deco(fn):
        BEHAVIOURS[name] = (fn, check, base)
        return fn
    return deco


def _lookup(table, name, rules, sim):
    """(entry, arg) for an exact or prefix name whose check passes, else (None, None)."""
    if name in table:
        return table[name], ""
    for key, entry in table.items():
        if key[-1] in "_:" and name.startswith(key) and len(name) > len(key):
            arg, check = name[len(key):], entry[1]
            if check is None or check(rules, sim, arg):
                return entry, arg
    return None, None


# ---------------------------------------------------------------------- built-in signals / behaviours
@register_signal("pressure")
def _pressure(rules, sim, arg, cache):
    return rules.pressure(sim)


@register_signal("neighbours")
def _neighbours(rules, sim, arg, cache):
    return rules.adjacency(sim, cache) @ np.ones(len(sim.tissue.pos))


@register_signal("contact_", check=lambda rules, sim, arg: arg in rules.types)
def _contact(rules, sim, arg, cache):
    return rules.adjacency(sim, cache) @ (sim.tissue.type_idx == rules.types.index(arg)).astype(float)


def is_gene(rules, sim, arg):
    return arg in sim.engine.gi


@register_behaviour("gene:", check=is_gene)
def _gene(rules, sim, arg, value, out):
    if out.get("gene_mod") is None:
        out["gene_mod"] = np.zeros((len(value), len(sim.engine.genes)))
    out["gene_mod"][:, sim.engine.gi[arg]] += value


@register_behaviour("chemotaxis:", check=is_gene)
def _chemotaxis(rules, sim, arg, value, out):
    b = np.clip(value, 0, 1)                                       # several ligands: strengths combine like
    out["bias"] = 1 - (1 - out.get("bias", 0.0)) * (1 - b)         # independent rules, directions add up
    out["bias_dir"] = out.get("bias_dir", 0.0) + b[:, None] * rules.gradient(arg, sim)


# ---------------------------------------------------------------------- rules
class Rules:
    KEYS = {"cell_type", "signal", "behaviour", "direction", "max", "min", "half_max", "hill_power",
            "enabled", "name"}
    DEFAULTS = {"cell_type": "*", "hill_power": 2, "enabled": True}

    def __init__(self, cfg, sim):
        cfg = cfg or {}
        eng, tis = sim.engine, sim.tissue
        self.types = list(tis.types)
        self.chemo_radius = float(cfg.get("chemotaxis_radius_um", 40.0))
        for path in cfg.get("plugins") or []:
            load_plugin(path if os.path.isabs(path) else os.path.join(sim.root, path))
        genes = pd.read_csv(f"{sim.root}/{sim.engine_cfg['network_dir']}/genes.tsv", sep="\t")
        self.modules = {}
        for g, mods in zip(genes.gene, genes.modules.fillna("")):
            for item in filter(None, mods.split(";")):
                m, role = item.split(":")
                self.modules.setdefault(m, {"pro": [], "anti": []})[role].append(eng.gi[g])

        self.rules, errors = [], []
        for k, raw in enumerate(cfg.get("rules") or []):
            r, err = self._parse(raw, sim)
            if err:
                errors.append(f"rule {k + 1} {raw}:" + "".join(f"\n      - {e}" for e in err))
            if r is not None and not err:
                self.rules.append(r)
        if errors:
            raise ValueError("rules.yaml problems:\n  " + "\n  ".join(errors))
        self.by_behaviour = {}
        for r in self.rules:
            self.by_behaviour.setdefault(r["behaviour"], []).append(r)

        self.pressure_ref = 1.0
        p = self.pressure(sim)
        self.pressure_ref = float(np.median(p[p > 0])) if (p > 0).any() else 1.0
        # baseline of every gene / module signal used: per-type mean and SD at the current (burnt-in) state,
        # so these signals are in "baseline standard deviations above the type's usual level"
        self.base_stats, cache = {}, {}
        for name in {r["signal"] for r in self.rules if r["signal_kind"] == "expression"}:
            v = self._raw(name, sim, cache)
            mu = np.array([v[tis.type_idx == t].mean() if (tis.type_idx == t).any() else 0 for t in range(len(self.types))])
            sd = np.array([v[tis.type_idx == t].std() if (tis.type_idx == t).any() else 1 for t in range(len(self.types))])
            self.base_stats[name] = (mu, np.maximum(sd, 1e-6))

    def _parse(self, raw, sim):
        """Fill defaults and validate one rule row. Returns (rule or None if disabled, [errors])."""
        r, err = {**self.DEFAULTS, **raw}, []
        unknown = set(r) - self.KEYS
        if unknown:
            err.append(f"unknown key(s) {sorted(unknown)}; allowed: {sorted(self.KEYS)}")
        if not r["enabled"]:
            return None, []
        if "direction" not in r:                                  # infer from which target is given
            r["direction"] = "increases" if "max" in r and "min" not in r else "decreases" if "min" in r and "max" not in r else None
        if r["direction"] not in ("increases", "decreases") or ("max" if r["direction"] == "increases" else "min") not in r:
            err.append("direction must be increases (with max) or decreases (with min)")
        else:
            r["target"] = float(r["max"] if r["direction"] == "increases" else r["min"])
        err += [f"missing {k}" for k in ("signal", "behaviour", "half_max") if k not in r]
        if "half_max" in r and not float(r["half_max"]) > 0:
            err.append("half_max must be > 0")
        ct = r["cell_type"]
        if ct == "*":
            r["types"] = None
        else:
            ct = [ct] if isinstance(ct, str) else list(ct)
            bad = [t for t in ct if t not in self.types]
            if bad:
                err.append(f"unknown cell_type {bad}; types are {self.types} or '*'")
            r["types"] = [self.types.index(t) for t in ct if t in self.types]
        s, b = r.get("signal"), r.get("behaviour")
        if s is None:
            pass
        elif _lookup(SIGNALS, s, self, sim)[0] is not None:
            r["signal_kind"] = "registered"
        elif s in self.modules or s in sim.engine.gi:
            r["signal_kind"] = "expression"
        else:
            reg = sorted(k + ("<...>" if k[-1] in "_:" else "") for k in SIGNALS)
            err.append(f"unknown signal {s!r}; use a module {sorted(self.modules)}, a network gene, or {reg}")
        if b is not None and b not in sim.tissue.p and _lookup(BEHAVIOURS, b, self, sim)[0] is None:
            reg = sorted(k + "<...>" if k[-1] in "_:" else k for k in BEHAVIOURS)
            err.append(f"unknown behaviour {b!r}; use a tissue parameter {sorted(sim.tissue.p)} or {reg} "
                       f"(gene: / chemotaxis: need a network gene)")
        return r, err

    # ------------------------------------------------------------------ signal helpers
    def adjacency(self, sim, cache):
        if "_A" not in cache:
            cache["_A"] = sim.tissue.adjacency()
        return cache["_A"]

    def pressure(self, sim):
        """Sum over overlapping neighbours of (1 - d / (R_i + R_j))^2, relative to the starting tissue."""
        tis = sim.tissue
        R = tis.par("radius")
        pr = tis._pairs(2 * R.max())
        out = np.zeros(len(tis.pos))
        if len(pr):
            i, j = pr[:, 0], pr[:, 1]
            d = np.linalg.norm(tis.pos[i] - tis.pos[j], axis=1)
            x = np.maximum(1 - d / (R[i] + R[j]), 0) ** 2
            out = np.bincount(i, x, len(out)) + np.bincount(j, x, len(out))
        return out / self.pressure_ref

    def gradient(self, gene, sim):
        """Unit direction of increasing `gene` expression, estimated from cells within chemotaxis_radius_um:
        sum over neighbours j of (L_j - L_i) * unit(x_j - x_i). Dying cells count as not expressing."""
        tis, eng = sim.tissue, sim.engine
        L = eng.S[:, eng.gi[gene]] * ~tis.dying
        g = np.zeros_like(tis.pos)
        pr = tis._pairs(self.chemo_radius)
        if len(pr):
            i, j = pr[:, 0], pr[:, 1]
            u = tis.pos[j] - tis.pos[i]
            u /= np.maximum(np.linalg.norm(u, axis=1), 1e-9)[:, None]
            c = (L[j] - L[i])[:, None] * u                     # same vector for both ends of the pair
            for k in (0, 1):
                g[:, k] = np.bincount(i, c[:, k], len(L)) + np.bincount(j, c[:, k], len(L))
        return g / np.maximum(np.linalg.norm(g, axis=1), 1e-12)[:, None]

    def _raw(self, name, sim, cache):
        """Gene / module signal in log2 units vs the type baseline."""
        eng = sim.engine
        if "_log2fc" not in cache:
            cache["_log2fc"] = np.log2((eng.S + 1) / (eng.ref[eng.type_idx] + 1))
        lf = cache["_log2fc"]
        if name in self.modules:
            m = self.modules[name]
            return (lf[:, m["pro"]].mean(1) if m["pro"] else 0) - (lf[:, m["anti"]].mean(1) if m["anti"] else 0)
        return lf[:, eng.gi[name]]

    def signal(self, name, sim, cache=None):
        """Per-cell value of a signal (cached within one rules evaluation)."""
        cache = {} if cache is None else cache
        if name not in cache:
            entry, arg = _lookup(SIGNALS, name, self, sim)
            if entry is not None:
                v = entry[0](self, sim, arg, cache)
            else:
                mu, sd = self.base_stats.get(name, (0.0, 1.0))
                ti = sim.tissue.type_idx
                v = (self._raw(name, sim, cache) - (mu[ti] if np.ndim(mu) else mu)) / (sd[ti] if np.ndim(sd) else sd)
            cache[name] = np.asarray(v, float) * np.ones(len(sim.tissue.pos))
        return cache[name]

    # ------------------------------------------------------------------ evaluate
    def __call__(self, sim):
        """Returns {"params": {tissue parameter: per-cell value}, "gene_mod": (cells x genes) or None,
        "bias"/"bias_dir": chemotaxis (if any), "values": {behaviour: per-cell value}} for sim.step."""
        tis = sim.tissue
        n, cache = len(tis.pos), {}
        out = {"params": {}, "gene_mod": None, "values": {}}
        for b, rs in self.by_behaviour.items():
            up, down = [], []
            for r in rs:
                H = hill(self.signal(r["signal"], sim, cache), float(r["half_max"]), float(r["hill_power"]))
                if r["types"] is not None:
                    H = H * np.isin(tis.type_idx, r["types"])
                (up if r["direction"] == "increases" else down).append((H, np.full(n, r["target"])))
            if b in tis.p:
                val = combine(tis.type_par(b).astype(float), up, down)
                out["params"][b] = val
            else:
                (fn, _, base), arg = _lookup(BEHAVIOURS, b, self, sim)
                val = combine(np.asarray(base(self, sim, arg), float) if base else np.zeros(n), up, down)
                fn(self, sim, arg, val, out)
            out["values"][b] = val
        self.last_signals = {k: v for k, v in cache.items() if not k.startswith("_")}
        return out


_LOADED = set()


def load_plugin(path):
    """Import a plugin file once; its @register_signal / @register_behaviour functions become available."""
    path = os.path.abspath(path)
    if path in _LOADED:
        return
    if not os.path.exists(path):
        raise FileNotFoundError(f"rules plugin not found: {path}")
    spec = importlib.util.spec_from_file_location(f"grn_abm_plugin_{len(_LOADED)}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    _LOADED.add(path)
