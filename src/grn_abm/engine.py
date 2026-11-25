"""Gene engine: per-cell gene state stepped forward in time, driven by a GRN and ligand-receptor signals.

State per cell and gene: promoter ON/OFF, unspliced u, spliced s (molecules).

One step of dt minutes:
  x        = s of every TF (mRNA as a proxy for TF activity)
  H        = x^n / (x^n + h^n), h = the TF's mean reference level across cell types   (SERGIO-style)
  grn_in   = H @ W, with W[tf, gene] = sign * weight from edges.tsv
  lig      = sum over neighbours of the ligand level (complex = weakest subunit)
  sig_p    = lig_p / (lig_p + h_lig) * rec_p / (rec_p + h_rec)    receptor-gated, per L-R pair
             h_lig = a typical neighbourhood's ligand input (ligand_half_response: neighbourhood), h_rec =
             the receptor's mean level across types
  lr_in    = sig @ M, with M[pair, tf] = lr_gain * sign * weight from lr_downstream.tsv
             (sign -1 = engagement represses the TF, e.g. cadherin contact -> TCF7/LEF1)
  k_on     = k_on0[type] * exp(grn_in + lr_in + mod)              mod = ABM rule input (log scale)
  promoter switches with the exact two-state transition probabilities over dt
  du/dt = alpha * ON - beta u,   ds/dt = beta u - gamma s         solved exactly over dt
Velocity = beta u - gamma s. Observed counts = Poisson(capture * s), sampled only when asked for.

Cell types and their baseline expression come from the config: any number of types, with baselines
either from a reference table (`baseline.source: reference`) or generated (`synthetic`).

Calibration: alpha (transcription rate while ON) is set per gene so the highest-expressing type sits at
`max_on_fraction`. k_on0[type, gene] is then chosen so that, with every regulator at that type's
reference level (and, if the starting layout is given, that type's average ligand input in it), the
mean of s equals the reference level. Cell types therefore start at their reference expression;
changes in neighbourhood and ABM rules move them away from it.
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


class GeneEngine:
    def __init__(self, cfg, root="."):
        self.cfg = cfg
        nd = os.path.join(root, cfg["network_dir"])
        genes = pd.read_csv(f"{nd}/genes.tsv", sep="\t")
        edges = pd.read_csv(f"{nd}/edges.tsv", sep="\t")
        lr = pd.read_csv(f"{nd}/lr_pairs.tsv", sep="\t")
        ds = pd.read_csv(f"{nd}/lr_downstream.tsv", sep="\t")

        self.genes = genes.gene.tolist()
        self.types = list(cfg["cell_types"])
        G, gi = len(self.genes), {g: i for i, g in enumerate(self.genes)}
        self.gi = gi
        self.rng = np.random.default_rng(cfg["seed"])
        self.dt = cfg["dt_min"]

        # ---- kinetics
        self.gamma = np.log(2) / (genes.halflife_h.values * 60.0)                  # decay, per min
        self.beta = np.full(G, np.log(2) / cfg["splice_halftime_min"])            # splicing, per min
        self.k_off = cfg["burst_off_rate"]

        # ---- baseline levels (molecules) per type: the calibration target
        self.baseline_counts = baseline_counts(cfg, self.genes, root)             # T x G, observed-count scale
        self.ref = np.maximum(self.baseline_counts / cfg["capture"], cfg["min_level"])

        # ---- GRN
        self.tf = np.array(sorted({gi[r] for r in edges.regulator}))
        tf_pos = {g: k for k, g in enumerate(self.tf)}
        W = np.zeros((len(self.tf), G))
        for r in edges.itertuples():
            W[tf_pos[gi[r.regulator]], gi[r.target]] += r.sign * r.weight
        self.W = W
        self.h_tf = self.ref[:, self.tf].mean(0)                                 # Hill half-response
        self.n = cfg["hill_n"]

        # ---- ligand-receptor
        sub = lambda x: [gi[s] for s in str(x).split("_") if s in gi]
        self.lig_sub = [sub(x) for x in lr.ligand]
        self.rec_sub = [sub(x) for x in lr.receptor]
        self.h_lig_cell = np.array([self.ref[:, s].min(1).mean() for s in self.lig_sub])   # one cell's ligand
        self.h_lig = self.h_lig_cell.copy()
        self.h_rec = np.array([self.ref[:, s].min(1).mean() for s in self.rec_sub])
        pair = {(a, b): p for p, (a, b) in enumerate(zip(lr.ligand, lr.receptor))}
        M = np.zeros((len(lr), G))
        signed = "sign" in ds.columns            # older networks have no sign column: treat every link as +1
        for r in ds.itertuples():
            M[pair[(r.ligand, r.receptor)], gi[r.tf]] += (r.sign if signed else 1) * r.weight
        self.M = M * cfg["lr_gain"]

        # ---- calibration (no growth dilution until init_cells says otherwise)
        self.grn_ref = self._grn_input(self.ref)                                  # T x G
        self._calibrate(np.zeros(len(self.types)))

    def _calibrate(self, mu):
        """Set alpha, the ON fraction and k_on0 so each type's mean s equals its baseline, given a per-type
        growth dilution rate mu (per min): mean u = p alpha / (beta + mu), mean s = beta u / (gamma + mu)."""
        self.mu = np.asarray(mu, float)
        a = self.beta[None, :] + self.mu[:, None]                                 # T x G, u loss
        c = self.gamma[None, :] + self.mu[:, None]                                # T x G, s loss
        need = self.ref * a * c / self.beta                                       # = p_on * alpha
        self.alpha = need.max(0) / self.cfg["max_on_fraction"]                    # per gene
        p_on = np.clip(need / self.alpha, 1e-6, 1 - 1e-6)
        self.p_on_ref = p_on
        self.k_on_star = p_on * self.k_off / (1 - p_on)                           # on-rate giving the baseline mean
        self.k_on0 = self.k_on_star / np.exp(self.grn_ref)

    # ------------------------------------------------------------------ pieces
    def _grn_input(self, S):
        x = S[:, self.tf]
        H = x ** self.n / (x ** self.n + self.h_tf ** self.n)
        return H @ self.W

    def _lr_input(self, S, A):
        lig = np.stack([S[:, s].min(1) for s in self.lig_sub], 1)               # N x P
        rec = np.stack([S[:, s].min(1) for s in self.rec_sub], 1)
        lig_nb = A @ lig                                                         # summed over neighbours
        sig = lig_nb / (lig_nb + self.h_lig) * rec / (rec + self.h_rec)
        return sig @ self.M, sig

    # ------------------------------------------------------------------ API
    def init_cells(self, type_idx, A=None, growth=None):
        """type_idx: int array (N,) indexing cfg['cell_types']. Starts at each type's reference state.

        A: starting adjacency. Each type's average ligand input in this layout is folded into the
        calibration, so the starting tissue sits at baseline and only CHANGES in a cell's neighbourhood
        move its genes. Without A, calibration assumes no signalling.
        growth: per-type division rate (per min). Dividing halves a cell's mRNA, so a growing population
        sits below its no-growth level for long-lived genes; calibration compensates for this.

        Each cell's u and s are drawn from the steady-state distribution of the bursting model (gamma with
        its mean and variance), so genes start fully dispersed whatever their half-life; observed counts
        at t = 0 are then Poisson-gamma = negative binomial.
        """
        self.type_idx = np.asarray(type_idx)
        if growth is not None:
            self._calibrate(growth)
        if A is not None:
            # ligand half-response: a typical NEIGHBOURHOOD's ligand input (one cell's level x mean contacts), so
            # the signal sits mid-range and responds to ligand dose; "cell" = one cell's level (saturates)
            if self.cfg.get("ligand_half_response", "neighbourhood") == "neighbourhood":
                nb = np.asarray(A.sum(1)).ravel()
                self.h_lig = self.h_lig_cell * max(nb[nb > 0].mean() if (nb > 0).any() else 1.0, 1.0)
            lr_ref, _ = self._lr_input(self.ref[self.type_idx], A)
            for t in range(len(self.types)):
                # mean of exp(input), not exp(mean input): ligand exposure is uneven across cells and
                # the rate is exponential in it, so averaging the input first under-corrects
                boost = np.exp(lr_ref[self.type_idx == t]).mean(0) if (self.type_idx == t).any() else 1.0
                self.k_on0[t] = self.k_on_star[t] / (np.exp(self.grn_ref[t]) * boost)
        p = self.p_on_ref[self.type_idx]
        self.on = self.rng.random(p.shape) < p
        mu = self.mu[self.type_idx][:, None]
        a, c = self.beta + mu, self.gamma + mu
        k = self.k_on_star[self.type_idx] + self.k_off
        var_on = self.alpha ** 2 * p * (1 - p)
        self.U = self._gamma_draw(p * self.alpha / a, var_on / (a * (a + k)))
        prod_s = self.alpha * self.beta / a                                       # effective production into s
        self.S = self._gamma_draw(p * prod_s / c, prod_s ** 2 * p * (1 - p) / (c * (c + k)))
        return self

    def _gamma_draw(self, mean, var):
        shape = mean ** 2 / np.maximum(var, 1e-12)
        return self.rng.gamma(np.maximum(shape, 1e-6), np.maximum(var, 1e-12) / np.maximum(mean, 1e-12))

    def recalibrate(self, max_change=4.0):
        """Scale each type's per-gene on-rate so its current mean s returns to the baseline. Used during the
        burn-in to remove the bias from fluctuating TF levels (the rate is exponential in TF input)."""
        for t in range(len(self.types)):
            cells = self.type_idx == t
            if cells.any():
                r = (self.ref[t] * 1.0) / np.maximum(self.S[cells].mean(0), 1e-9)
                self.k_on0[t] *= np.clip(r, 1 / max_change, max_change)

    def step(self, A, mod=None):
        """A: N x N sparse adjacency from the ABM. mod: optional N x G log-rate input from ABM rules."""
        S, U = self.S, self.U
        lr_in, self.last_sig = self._lr_input(S, A)
        logk = self._grn_input(S) + lr_in
        if mod is not None:
            logk = logk + mod
        k_on = self.k_on0[self.type_idx] * np.exp(np.clip(logk, -self.cfg["max_log_rate_change"],
                                                          self.cfg["max_log_rate_change"]))
        # exact two-state promoter switching over dt
        tot = k_on + self.k_off
        relax = np.exp(-tot * self.dt)
        p_on_next = np.where(self.on, 1 - (self.k_off / tot) * (1 - relax), (k_on / tot) * (1 - relax))
        on_prev = self.on
        self.on = self.rng.random(p_on_next.shape) < p_on_next
        # mRNA: production uses the promoter state averaged over the step (start and end)
        prod = self.alpha * 0.5 * (on_prev.astype(float) + self.on)
        # du/dt = prod - beta u, ds/dt = beta u - gamma s, solved exactly over dt (growth dilution happens
        # through explicit divisions, which halve the mRNA)
        a, c = self.beta, self.gamma
        u_inf = prod / a
        s_inf = self.beta * u_inf / c
        ea, ec = np.exp(-a * self.dt), np.exp(-c * self.dt)
        du = U - u_inf
        self.U = u_inf + du * ea
        self.S = s_inf + du * (self.beta / (c - a)) * (ea - ec) + (S - s_inf) * ec

    def velocity(self):
        return self.beta * self.U - self.gamma * self.S

    def observe(self):
        return self.rng.poisson(self.cfg["capture"] * np.maximum(self.S, 0))

    def _split(self, X):
        n = np.floor(X)
        return self.rng.binomial(n.astype(np.int64), 0.5) + (X - n) / 2

    def halve(self, idx):
        """Virtual division (burn-in): the cell keeps one daughter's share of its mRNA, cell count unchanged."""
        self.S[idx] -= self._split(self.S[idx])
        self.U[idx] -= self._split(self.U[idx])

    def divide(self, parents):
        """Daughters are appended; mRNA is split binomially, promoter state copied."""
        split = self._split
        ks, ku = split(self.S[parents]), split(self.U[parents])
        self.S[parents] -= ks
        self.U[parents] -= ku
        self.S, self.U = np.vstack([self.S, ks]), np.vstack([self.U, ku])
        self.on = np.vstack([self.on, self.on[parents]])
        self.type_idx = np.concatenate([self.type_idx, self.type_idx[parents]])

    def remove(self, idx):
        keep = np.ones(len(self.type_idx), bool)
        keep[idx] = False
        self.S, self.U, self.on, self.type_idx = self.S[keep], self.U[keep], self.on[keep], self.type_idx[keep]


def baseline_counts(cfg, genes, root="."):
    """Each cell type's mean observed counts per gene (T x G), from the reference or generated."""
    types, b = list(cfg["cell_types"]), cfg["baseline"]
    if b["source"] == "reference":
        expr = pd.read_csv(os.path.join(root, b["reference"]["table"]), sep="\t", index_col=0)
        missing = [t for t in types if f"mean_{t}" not in expr]
        if missing:
            raise KeyError(f"no mean_<type> column for {missing} in {b['reference']['table']}; rerun scripts/03")
        return expr.reindex(genes)[[f"mean_{t}" for t in types]].fillna(0).values.T
    if b["source"] == "synthetic":
        s = b["synthetic"]
        r = np.random.default_rng(s["seed"])
        gene_log = r.normal(np.log(s["median_counts"]), s["gene_sd_log"], len(genes))
        return np.exp(gene_log + r.normal(0, s["type_sd_log"], (len(types), len(genes))))
    raise ValueError(f"baseline.source must be reference or synthetic, not {b['source']!r}")


def load_engine(config_path, root):
    return GeneEngine(yaml.safe_load(open(config_path)), root=root)
