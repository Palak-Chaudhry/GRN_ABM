"""Tissue layer: cells as discs with PhysiCell-style mechanics, motility, division and death.

Mechanics (per pair i, j at distance d, overdamped, PhysiCell form):
  repulsion  sqrt(rep_i rep_j) (1 - d / (R_i + R_j))^2                     if d < R_i + R_j
  adhesion   sqrt(adh_i adh_j) (1 - d / (f_i R_i + f_j R_j))^2              if d < f_i R_i + f_j R_j
  velocity_i += (repulsion - adhesion) * (x_i - x_j) / d
  position   += dt (1.5 v - 0.5 v_prev)                                      (Adams-Bashforth 2)
Motility: velocity += speed * direction; direction is redrawn with probability dt / persistence.
Division: daughter placed half a radius away in a random direction (radius kept, no growth).
Death: a dying cell stops moving and is excluded from signalling; it is removed after its duration.

Parameters: every key under abm.yaml `defaults:` is a parameter, set per type under `types:` and readable
per cell with tissue.par(name). Rules can override any of them per cell for one step
(set_cell_params), so adding a key to `defaults:` makes it rule-controllable with no code change; it only
has an effect once code reads it with par(name). The keys in REQUIRED are the ones this file reads.

Rows are cells; births are appended, removals delete rows. The gene engine keeps the same row order
(see sim.py), so the two layers always line up.
"""
import os
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree

REQUIRED = ["radius", "repulsion", "adhesion", "max_adhesion_distance", "speed", "persistence_min",
            "division_rate", "apoptosis_rate", "apoptosis_duration_min"]


class Tissue:
    def __init__(self, cfg, type_names, seed=0):
        self.cfg, self.types = cfg, list(type_names)
        self.rng = np.random.default_rng(seed)
        defaults, over = cfg["defaults"], cfg.get("types") or {}
        errors = [f"defaults: missing {k}" for k in REQUIRED if k not in defaults]
        errors += [f"types: unknown cell type {t!r} (types are {self.types})" for t in over if t not in self.types]
        errors += [f"types.{t}: unknown parameter {k!r} (add it to defaults first)"
                   for t, o in over.items() for k in (o or {}) if k not in defaults]
        if errors:
            raise ValueError("abm.yaml problems:\n  " + "\n  ".join(errors))
        # per-type parameter vectors, e.g. self.p["radius"][type]
        self.p = {k: np.array([float((over.get(t) or {}).get(k, v)) for t in self.types]) for k, v in defaults.items()}
        self.cell = {}                          # per-cell overrides for the current step (from rules)
        self.bias = self.bias_dir = None        # per-cell chemotaxis strength (0-1) and unit direction
        self.domain = None                      # ("disc", r) or ("box", xmin, xmax, ymin, ymax)
        self.next_id = 0

    # ------------------------------------------------------------------ setup
    def _set_cells(self, pos, type_idx):
        n = len(pos)
        self.pos, self.type_idx = np.asarray(pos, float), np.asarray(type_idx, int)
        self.v_prev = np.zeros_like(self.pos)
        self.dir = self._random_dirs(n)
        self.dying = np.zeros(n, bool)
        self.death_timer = np.zeros(n)
        self.cell_id = np.arange(self.next_id, self.next_id + n)
        self.parent_id = np.full(n, -1)
        self.next_id += n
        self.cell, self.bias, self.bias_dir = {}, None, None

    def init_synthetic(self):
        s = self.cfg["init"]["synthetic"]
        counts = [int((s.get("cells") or {}).get(t, s.get("cells_default", 0))) for t in self.types]
        type_idx = np.repeat(np.arange(len(self.types)), counts)
        r = np.random.default_rng(s["seed"])
        r.shuffle(type_idx)
        area = np.sum(np.pi * self.p["radius"][type_idx] ** 2) / s["packing"]
        R = np.sqrt(area / np.pi)
        th, rr = r.uniform(0, 2 * np.pi, len(type_idx)), R * np.sqrt(r.uniform(0, 1, len(type_idx)))
        self.domain = ("disc", R * s.get("domain_scale", 1.0))
        self._set_cells(np.c_[rr * np.cos(th), rr * np.sin(th)], type_idx)
        return self

    def init_reference(self, engine_cfg):
        """Crop of the reference tissue: positions and types of the configured cell types only."""
        os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
        import anndata as ad
        ref, c = engine_cfg["baseline"]["reference"], self.cfg["init"]["reference"]
        a = ad.read_h5ad(f"{ref['dir']}/{ref['labels_h5ad']}", backed="r")
        labels = a.obs[ref.get("label_key", "cell_type")].astype(str).values
        xy = np.asarray(a.obsm[ref.get("spatial_key", "spatial")])[:, :2]
        lab2t = {lab: i for i, lab in enumerate(engine_cfg["cell_types"].values())}
        sel = np.array([l in lab2t for l in labels])
        h = c["crop_um"] / 2
        if c["center_um"] is None:
            # crop-sized grid square that best contains every configured type: max sum_t log(1 + count_t)
            ex, ey = (np.arange(xy[sel, k].min(), xy[sel, k].max() + c["crop_um"], c["crop_um"]) for k in (0, 1))
            tsel = np.array([lab2t[l] for l in labels[sel]])
            score = sum(np.log1p(np.histogram2d(*xy[sel][tsel == t].T, bins=[ex, ey])[0]) for t in range(len(lab2t)))
            ix, iy = np.unravel_index(np.argmax(score), score.shape)
            center = np.array([ex[ix] + h, ey[iy] + h])
        else:
            center = np.asarray(c["center_um"], float)
        win = sel & (np.abs(xy - center) <= h).all(1)
        self.domain = ("box", -h, h, -h, h)
        self._set_cells(xy[win] - center, [lab2t[l] for l in labels[win]])
        return self

    # ------------------------------------------------------------------ helpers
    def _random_dirs(self, n):
        th = self.rng.uniform(0, 2 * np.pi, n)
        return np.c_[np.cos(th), np.sin(th)]

    def par(self, name):
        """Per-cell value of a parameter: the rule override for this step if set, else the type's value."""
        return self.cell[name] if name in self.cell else self.p[name][self.type_idx]

    def type_par(self, name):
        """Per-cell value of the type's parameter (ignores rule overrides)."""
        return self.p[name][self.type_idx]

    def set_cell_params(self, values=None):
        """Per-cell overrides {parameter: array} for the coming step; replaces any previous ones."""
        values = values or {}
        bad = [k for k in values if k not in self.p]
        if bad:
            raise KeyError(f"unknown tissue parameter(s) {bad}; parameters are {sorted(self.p)}")
        self.cell = {k: np.asarray(v, float) for k, v in values.items()}

    def set_chemotaxis(self, bias=None, direction=None):
        """Per-cell chemotaxis strength in [0, 1] and direction (normalised here), or None for none."""
        if bias is None:
            self.bias = self.bias_dir = None
            return
        d = np.asarray(direction, float)
        self.bias = np.clip(np.asarray(bias, float), 0, 1)
        self.bias_dir = d / np.maximum(np.linalg.norm(d, axis=1), 1e-12)[:, None]

    def _per_cell(self):
        """Names of every per-cell array (kept in step through births and removals)."""
        return ["pos", "type_idx", "v_prev", "dir", "dying", "death_timer", "cell_id", "parent_id"]

    def _resize_overrides(self, rows):
        """Overrides / chemotaxis follow the cells: rows = parent rows of the kept + new cells."""
        self.cell = {k: v[rows] for k, v in self.cell.items()}
        if self.bias is not None:
            self.bias, self.bias_dir = self.bias[rows], self.bias_dir[rows]

    def _pairs(self, cutoff):
        return cKDTree(self.pos).query_pairs(cutoff, output_type="ndarray")

    # ------------------------------------------------------------------ dynamics
    def mechanics_step(self, dt, motility=True):
        """Every parameter is read per cell (par), so rule overrides apply. With chemotaxis set, the motile
        direction is normalize(bias * bias_dir + (1 - bias) * random persistent direction)."""
        n = len(self.pos)
        R, f = self.par("radius"), self.par("max_adhesion_distance")
        pr = self._pairs(2 * (R * f).max())
        v = np.zeros((n, 2))
        if len(pr):
            i, j = pr[:, 0], pr[:, 1]
            dvec = self.pos[i] - self.pos[j]
            d = np.maximum(np.linalg.norm(dvec, axis=1), 1e-9)
            rs, ra = R[i] + R[j], f[i] * R[i] + f[j] * R[j]
            rep = np.where(d < rs, np.sqrt(self.par("repulsion")[i] * self.par("repulsion")[j]) * (1 - d / rs) ** 2, 0)
            adh = np.where(d < ra, np.sqrt(self.par("adhesion")[i] * self.par("adhesion")[j]) * (1 - d / ra) ** 2, 0)
            fv = dvec * ((rep - adh) / d)[:, None]
            for k in (0, 1):
                v[:, k] = np.bincount(i, fv[:, k], n) - np.bincount(j, fv[:, k], n)
        if motility:
            speed = np.where(self.dying, 0.0, self.par("speed"))
            redraw = self.rng.random(n) < dt / self.par("persistence_min")
            self.dir[redraw] = self._random_dirs(int(redraw.sum()))
            d = self.dir
            if self.bias is not None:
                d = self.bias[:, None] * self.bias_dir + (1 - self.bias)[:, None] * self.dir
                d = d / np.maximum(np.linalg.norm(d, axis=1), 1e-12)[:, None]
            v += speed[:, None] * d
        self.pos += dt * (1.5 * v - 0.5 * self.v_prev)
        self.v_prev = v
        self._confine()

    def _confine(self):
        if self.domain is None:
            return
        if self.domain[0] == "disc":
            r = np.linalg.norm(self.pos, axis=1)
            out = r > self.domain[1]
            self.pos[out] *= (self.domain[1] / r[out])[:, None]
        else:
            _, x0, x1, y0, y1 = self.domain
            self.pos[:, 0] = np.clip(self.pos[:, 0], x0, x1)
            self.pos[:, 1] = np.clip(self.pos[:, 1], y0, y1)

    def adjacency(self):
        """Signalling contacts (both directions) between living, non-dying cells."""
        n = len(self.pos)
        R, f = self.par("radius"), self.par("max_adhesion_distance")
        pr = self._pairs(self.cfg["contact_factor"] * 2 * (R * f).max())
        if len(pr):
            i, j = pr[:, 0], pr[:, 1]
            d = np.linalg.norm(self.pos[i] - self.pos[j], axis=1)
            keep = (d < self.cfg["contact_factor"] * (f[i] * R[i] + f[j] * R[j])) & ~self.dying[i] & ~self.dying[j]
            i, j = i[keep], j[keep]
        else:
            i = j = np.array([], int)
        return sp.csr_matrix((np.ones(2 * len(i)), (np.r_[i, j], np.r_[j, i])), shape=(n, n))

    def draw_events(self, dt):
        """Poisson events over dt at the per-cell division_rate and apoptosis_rate (per min)."""
        div, apo = self.par("division_rate"), self.par("apoptosis_rate")
        alive = ~self.dying
        start_dying = np.where(alive & (self.rng.random(len(self.pos)) < 1 - np.exp(-apo * dt)))[0]
        parents = np.where(alive & (self.rng.random(len(self.pos)) < 1 - np.exp(-div * dt)))[0]
        parents = np.setdiff1d(parents, start_dying)
        return parents, start_dying

    def divide(self, parents):
        k = len(parents)
        if k == 0:
            return
        u = self._random_dirs(k)
        off = 0.5 * self.par("radius")[parents][:, None] * u
        self._resize_overrides(np.r_[np.arange(len(self.pos)), parents])     # daughters inherit overrides
        kids = self.pos[parents] + off
        self.pos[parents] -= off
        self.pos = np.vstack([self.pos, kids])
        self.type_idx = np.r_[self.type_idx, self.type_idx[parents]]
        self.v_prev = np.vstack([self.v_prev, np.zeros((k, 2))])
        self.dir = np.vstack([self.dir, self._random_dirs(k)])
        self.dying = np.r_[self.dying, np.zeros(k, bool)]
        self.death_timer = np.r_[self.death_timer, np.zeros(k)]
        self.parent_id = np.r_[self.parent_id, self.cell_id[parents]]
        self.cell_id = np.r_[self.cell_id, np.arange(self.next_id, self.next_id + k)]
        self.next_id += k
        self._confine()

    def start_dying(self, idx):
        self.dying[idx] = True
        self.death_timer[idx] = self.par("apoptosis_duration_min")[idx]

    def advance_death(self, dt):
        """Counts down dying cells; returns the rows whose time is up (to be removed)."""
        self.death_timer[self.dying] -= dt
        return np.where(self.dying & (self.death_timer <= 0))[0]

    def remove(self, idx):
        keep = np.ones(len(self.pos), bool)
        keep[idx] = False
        for k in self._per_cell():
            setattr(self, k, getattr(self, k)[keep])
        self._resize_overrides(np.where(keep)[0])
