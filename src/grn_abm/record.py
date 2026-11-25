"""Records a running simulation for analysis. Use as the run callback:

    rec = Recorder(sim, snapshot_min=180, track_min=6)      # records t = 0 immediately
    sim.run(72, callback=rec)                                # called every step
    rec.save("results/runs/my_run")

Writes
  snapshots.h5ad   every snapshot_min: one row per living-or-dying cell per snapshot
                     X / layers['spliced'], layers['unspliced']  observed counts (Poisson(capture x mRNA))
                     layers['S_true'], layers['U_true']          true mRNA molecules (velocity = beta U - gamma S,
                                                                  beta / gamma per minute in var)
                     obs: time_h, cell_id, parent_id, birth_h, cell_type, dying, x, y, module_<m> (log2 vs
                          baseline), signal_<s> and rule_<behaviour> (what the rules saw and set), lr_* summary
                     obsm['spatial'], obsm['lr_received'] (cells x L-R pairs: receptor-gated signal received)
  tracks.npz       every track_min: time_h, cell_id, type, x, y, dying (for movement plots / animation)
  events.tsv       time_h, kind (division / death_start / removed), cell_id, type
  lr_flow.npz      per snapshot: flow[time, pair, sender type, receiver type] = sum over contacts of
                   ligand(sender) x receptor gate(receiver)
  grn_activity.npz per snapshot: tf_activity[time, type, TF] = mean Hill activation of each TF; W (TF x gene)
                   so edge activity = tf_activity x W
Observed counts use their own random stream, so recording never changes the simulation.
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


class Recorder:
    def __init__(self, sim, snapshot_min=180, track_min=None, seed=0):
        self.sim = sim
        self.rng = np.random.default_rng(seed + 12345)
        self.snap_k = max(1, int(round(snapshot_min / sim.dt)))
        self.track_k = max(1, int(round((track_min or sim.dt) / sim.dt)))
        eng = sim.engine
        nd = os.path.join(sim.root, sim.engine_cfg["network_dir"])
        self.genes_tab = pd.read_csv(f"{nd}/genes.tsv", sep="\t")
        lr = pd.read_csv(f"{nd}/lr_pairs.tsv", sep="\t")
        self.pairs = [f"{a}->{b}" for a, b in zip(lr.ligand, lr.receptor)]
        self.modules = {}
        for g, mods in zip(self.genes_tab.gene, self.genes_tab.modules.fillna("")):
            for item in filter(None, mods.split(";")):
                m, role = item.split(":")
                self.modules.setdefault(m, {"pro": [], "anti": []})[role].append(eng.gi[g])
        self.birth = {int(c): 0.0 for c in sim.tissue.cell_id}
        self.snaps, self.tracks, self.flow, self.tf_act, self.snap_t = [], [], [], [], []
        self.k = 0
        self._track()
        self._snapshot()

    def __call__(self, sim):
        self.k += 1
        for c in sim.tissue.cell_id:
            c = int(c)
            if c not in self.birth:
                self.birth[c] = (sim.t - sim.dt) / 60          # born during the step that just ran
        if self.k % self.track_k == 0:
            self._track()
        if self.k % self.snap_k == 0:
            self._snapshot()

    # ------------------------------------------------------------------ pieces
    def _track(self):
        t = self.sim.tissue
        self.tracks.append((np.full(len(t.pos), self.sim.t / 60), t.cell_id.copy(), t.type_idx.copy(),
                            t.pos[:, 0].copy(), t.pos[:, 1].copy(), t.dying.copy()))

    def _snapshot(self):
        sim, eng, tis = self.sim, self.sim.engine, self.sim.tissue
        n, cap = len(tis.pos), sim.engine_cfg["capture"]
        types = list(tis.types)
        obs = pd.DataFrame({
            "time_h": sim.t / 60, "cell_id": tis.cell_id, "parent_id": tis.parent_id,
            "birth_h": [self.birth[int(c)] for c in tis.cell_id],
            "cell_type": pd.Categorical([types[k] for k in tis.type_idx], categories=types),
            "dying": tis.dying, "x": tis.pos[:, 0], "y": tis.pos[:, 1]})
        lf = np.log2((eng.S + 1) / (eng.ref[eng.type_idx] + 1))
        for m, d in self.modules.items():
            obs[f"module_{m}"] = (lf[:, d["pro"]].mean(1) if d["pro"] else 0) - (lf[:, d["anti"]].mean(1) if d["anti"] else 0)
        if sim.rules is not None:                                   # pure evaluation: no random draws
            out = sim.rules(sim)
            for b, v in out["values"].items():
                if np.ndim(v) == 1:
                    obs[f"rule_{b}"] = v
            for s, v in sim.rules.last_signals.items():
                obs[f"signal_{s}"] = v
        for b in ("division_rate", "apoptosis_rate", "speed"):      # value in force (rule or type)
            if f"rule_{b}" not in obs:
                obs[f"rule_{b}"] = tis.type_par(b)
        A = tis.adjacency()
        obs["neighbours"] = np.asarray(A.sum(1)).ravel()
        # L-R: received signal per cell and pair; flow between types
        lig = np.stack([eng.S[:, s].min(1) for s in eng.lig_sub], 1)
        rec = np.stack([eng.S[:, s].min(1) for s in eng.rec_sub], 1)
        gate = rec / (rec + eng.h_rec)
        _, received = eng._lr_input(eng.S, A)
        flow = np.zeros((len(self.pairs), len(types), len(types)))
        onehot = np.eye(len(types))[tis.type_idx]
        for b in range(len(types)):
            Z = A @ (gate * onehot[:, [b]])                          # receptor gate of type-b neighbours
            flow[:, :, b] = (onehot.T @ (lig * Z)).T
        self.flow.append(flow)
        x = eng.S[:, eng.tf]
        H = x ** eng.n / (x ** eng.n + eng.h_tf ** eng.n)
        self.tf_act.append(np.array([H[tis.type_idx == t].mean(0) if (tis.type_idx == t).any()
                                     else np.full(len(eng.tf), np.nan) for t in range(len(types))]))
        self.snap_t.append(sim.t / 60)
        a = ad.AnnData(X=sp.csr_matrix(self.rng.poisson(cap * np.maximum(eng.S, 0)).astype(np.float32)), obs=obs)
        a.layers["spliced"] = a.X.copy()
        a.layers["unspliced"] = sp.csr_matrix(self.rng.poisson(cap * np.maximum(eng.U, 0)).astype(np.float32))
        a.layers["S_true"] = eng.S.astype(np.float32)
        a.layers["U_true"] = eng.U.astype(np.float32)
        a.obsm["spatial"] = tis.pos.copy()
        a.obsm["lr_received"] = received.astype(np.float32)
        self.snaps.append(a)

    # ------------------------------------------------------------------ save
    def save(self, out):
        os.makedirs(out, exist_ok=True)
        sim, eng = self.sim, self.sim.engine
        a = ad.concat(self.snaps, index_unique=None)
        a.obs_names = [f"t{t:.1f}_c{c}" for t, c in zip(a.obs.time_h, a.obs.cell_id)]
        a.obs["cell_type"] = pd.Categorical(a.obs["cell_type"].astype(str), categories=list(sim.tissue.types))
        var = self.genes_tab.set_index("gene").loc[eng.genes].copy()
        var["modules"] = var["modules"].fillna("").astype(str)
        var["role"] = var["role"].astype(str)
        var["beta_per_min"], var["gamma_per_min"] = eng.beta, eng.gamma
        var["is_tf"] = np.isin(np.arange(len(eng.genes)), eng.tf)
        a.var = var
        a.var_names = list(eng.genes)
        a.uns["lr_pairs"] = np.array(self.pairs)
        a.uns["cell_types"] = np.array(sim.tissue.types)
        a.uns["capture"] = sim.engine_cfg["capture"]
        a.uns["baseline_counts"] = eng.baseline_counts.astype(np.float32)   # type x gene, observed scale
        a.write_h5ad(f"{out}/snapshots.h5ad", compression="gzip")
        tr = [np.concatenate(z) for z in zip(*self.tracks)]
        np.savez_compressed(f"{out}/tracks.npz", time_h=tr[0], cell_id=tr[1], type=tr[2], x=tr[3], y=tr[4],
                            dying=tr[5], types=np.array(sim.tissue.types))
        ev = pd.DataFrame(sim.events, columns=["time_min", "kind", "cell_id", "type"])
        ev["time_h"] = ev.pop("time_min") / 60
        ev["type"] = [sim.tissue.types[k] for k in ev["type"]]
        ev.to_csv(f"{out}/events.tsv", sep="\t", index=False)
        np.savez_compressed(f"{out}/lr_flow.npz", time_h=np.array(self.snap_t), flow=np.array(self.flow),
                            pairs=np.array(self.pairs), types=np.array(sim.tissue.types))
        np.savez_compressed(f"{out}/grn_activity.npz", time_h=np.array(self.snap_t), tf_activity=np.array(self.tf_act),
                            tfs=np.array([eng.genes[k] for k in eng.tf]), genes=np.array(eng.genes), W=eng.W,
                            types=np.array(sim.tissue.types))
        with open(f"{out}/run_info.yaml", "w") as f:
            yaml.safe_dump({"engine": sim.engine_cfg, "abm": sim.abm_cfg, "hours": float(sim.t / 60),
                            "snapshots": len(self.snaps), "cells_recorded": int(a.n_obs)}, f, sort_keys=False)
        return a
