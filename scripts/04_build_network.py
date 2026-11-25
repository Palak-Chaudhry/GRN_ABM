"""Build the simulator's network from any GRN + an L-R list (configs/network.yaml).

To use a new GRN: point `grn.path` at it and set `grn.format` (celloracle_wide or edge_list).
Nothing else changes. Missing signs / weights / receptor->TF links are filled with seeded defaults.

Steps:
  1. targets = all behaviour-module genes
  2. TFs     = module genes that are regulators, then the regulators with the most edges into targets
  3. edges   = GRN edges from the selected TFs to every selected gene (random cap per target if dense)
  4. L-R     = pairs involving a module gene first, then by ligand x receptor expression (if given)
  5. L-R -> TF links from `lr.downstream` if given, else random selected TFs

Cell-type differences are not decided here: in the gene engine a TF only acts in cells that express it.

Outputs (data/processed/network/): genes.tsv, edges.tsv, lr_pairs.tsv, lr_downstream.tsv, report.md
"""
import argparse
import os
import numpy as np
import pandas as pd
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--config", default=f"{ROOT}/configs/network.yaml")
ap.add_argument("--out", default=f"{ROOT}/data/processed/network")
ARGS = ap.parse_args()
CFG = yaml.safe_load(open(ARGS.config))
SZ, DEF = CFG["sizes"], CFG["defaults"]
OUT = ARGS.out
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(DEF["seed"])
log = []


def say(*a):
    log.append(" ".join(map(str, a))); print(log[-1])


def read_table(path):
    return pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path, sep=None, engine="python")


def load_grn(path, fmt):
    """Return an edge list with columns regulator, target (+ sign / weight if the source has them)."""
    df = read_table(path)
    if fmt == "celloracle_wide":
        tfs = [c for c in df.columns if c not in ("peak_id", "gene_short_name")]
        wide = df.groupby("gene_short_name")[tfs].max()
        e = wide.stack()
        e = e[e > 0].reset_index()
        e.columns = ["target", "regulator", "value"]
        return e[["regulator", "target"]]
    if fmt == "edge_list":
        return df
    raise ValueError(f"unknown grn.format {fmt!r}")


subunits = lambda x: str(x).split("_")

# ---- inputs
grn = load_grn(CFG["grn"]["path"], CFG["grn"]["format"]).drop_duplicates(["regulator", "target"])
_lr_cols = ["ligand", "receptor"] + [c for c in ("pathway", "signaling") if c in read_table(CFG["lr"]["path"])]
lr = read_table(CFG["lr"]["path"])[_lr_cols].drop_duplicates(["ligand", "receptor"]).reset_index(drop=True)
expr = None
if CFG.get("expression"):
    expr = pd.read_csv(os.path.join(ROOT, CFG["expression"]), sep="\t", index_col=0).filter(like="det_").max(axis=1)
modules = CFG["modules"]
targets = sorted({g for m in modules.values() for g in m["pro"] + m["anti"]})
say(f"GRN: {len(grn)} edges, {grn.regulator.nunique()} regulators, {grn.target.nunique()} targets "
    f"({CFG['grn']['format']}: {CFG['grn']['path']})")

# ---- 1. L-R pairs (before TFs: the curated pathway effectors of the SELECTED pairs get
#         reserved TF slots, so a curated receptor->TF link is never dropped by the degree ranking)
lr["genes"] = lr.ligand.map(subunits) + lr.receptor.map(subunits)
lr["module_gene"] = lr.genes.map(lambda gs: any(g in targets for g in gs))
if expr is not None:
    lr["score"] = lr.genes.map(lambda gs: min(expr.get(g, 0.0) for g in gs))   # weakest subunit
    lr = lr[lr.score > 0]
else:
    lr["score"] = 0.0
# module-gene pairs first, but at most `pairs_per_module_gene` per gene (ITGB1 alone matches ~40 ECM pairs)
lr = lr.sort_values("score", ascending=False)
first = set()
for g in targets:
    first |= set(lr.index[lr.genes.map(lambda gs: g in gs)][: SZ["pairs_per_module_gene"]])
lr["module_gene"] = lr.index.isin(first)
lr = lr.sort_values(["module_gene", "score"], ascending=False).head(SZ["n_lr_pairs"]).reset_index(drop=True)
say(f"L-R pairs: {len(lr)} ({int(lr.module_gene.sum())} chosen for a module gene, max {SZ['pairs_per_module_gene']} each; rest ranked by "
    f"{'expression' if expr is not None else 'file order'})")

# ---- 2. TFs: module genes, then curated L-R effectors, then the regulators with the most edges in
regulators = set(grn.regulator)
into_targets = grn[grn.target.isin(targets)].regulator.value_counts()
tfs = [g for g in targets if g in regulators]
n_module_tfs = len(tfs)
effectors = []
if CFG["lr"].get("downstream"):
    _D = read_table(os.path.join(ROOT, CFG["lr"]["downstream"])).merge(lr[["ligand", "receptor"]])
    want = [t for t in dict.fromkeys(_D.tf) if t in regulators and t not in tfs]
    dropped = sorted({t for t in _D.tf} - regulators)
    effectors = want[: max(0, SZ["n_tf"] - len(tfs))]
    tfs += effectors
    if dropped:
        say(f"  effector TFs absent from the GRN, skipped: {', '.join(dropped)}")
    if len(want) > len(effectors):
        say(f"  WARNING: n_tf too small; {len(want) - len(effectors)} curated effectors did not fit: "
            f"{', '.join(want[len(effectors):])}")
n_degree = max(0, SZ["n_tf"] - len(tfs))
tfs += [r for r in into_targets.index if r not in tfs][:n_degree]
say(f"TFs: {len(tfs)} ({n_module_tfs} module genes, {len(effectors)} curated L-R effectors, "
    f"{len(tfs) - n_module_tfs - len(effectors)} by most edges into module genes)")

# ---- 3. edges into every selected gene
genes = sorted(set(targets) | set(tfs) | {g for gs in lr.genes for g in gs})
E = grn[grn.regulator.isin(tfs) & grn.target.isin(genes)].copy()
cap = SZ["max_regulators_per_target"]
if cap:
    E = E.sample(frac=1, random_state=DEF["seed"]).groupby("target").head(cap)
if "sign" not in E:
    E["sign"] = np.where(rng.random(len(E)) < DEF["activator_fraction"], 1, -1)
if "weight" not in E:
    E["weight"] = rng.lognormal(*DEF["weight_lognormal"], len(E)).round(4)
E = E[["regulator", "target", "sign", "weight"]].sort_values(["target", "regulator"])
say(f"edges: {len(E)} into {E.target.nunique()} of {len(genes)} genes (cap {cap or 'none'} per target)")
say(f"module genes with no regulator: {sorted(set(targets) - set(E.target)) or 'none'}")

# ---- 5. L-R -> TF links
if CFG["lr"].get("downstream"):
    D = read_table(os.path.join(ROOT, CFG["lr"]["downstream"])).merge(lr[["ligand", "receptor"]])
    n_all = len(D)
    D = D[D.tf.isin(tfs)].copy()
    if "sign" not in D:
        D["sign"] = 1
    if "weight" not in D:
        D["weight"] = DEF["lr_weight"]
    D["weight"] = D.weight * DEF["lr_weight"]
    if n_all > len(D):
        say(f"  {n_all - len(D)} of {n_all} curated links dropped (TF not among the selected {len(tfs)})")
    pairs_hit = D[["ligand", "receptor"]].drop_duplicates().shape[0]
    say(f"L-R -> TF links: {len(D)} from file ({(D.sign < 0).sum()} inhibitory); "
        f"{pairs_hit}/{len(lr)} pairs have at least one")
else:
    D = pd.DataFrame([dict(ligand=r.ligand, receptor=r.receptor, tf=tf, sign=1, weight=DEF["lr_weight"])
                      for r in lr.itertuples() for tf in rng.choice(tfs, SZ["lr_downstream_tfs"], replace=False)])
    say(f"L-R -> TF links: {len(D)} (random, seeded)")
D = D[["ligand", "receptor"] + [c for c in ("pathway",) if c in D] + ["tf", "sign", "weight"]]

# ---- check: every module reachable from at least one ligand (TF -> TF -> ... -> gene)
reach, frontier = set(D.tf), set(D.tf)   # sign-agnostic: a repressive link still couples
while frontier:
    nxt = set(E[E.regulator.isin(frontier)].target) - reach
    reach |= nxt
    frontier = nxt & set(tfs)
for m, d in modules.items():
    gs = set(d["pro"] + d["anti"])
    say(f"  {m:14s} {len(gs & reach):2d}/{len(gs)} genes reachable from a ligand")

# ---- write
hl_path = os.path.join(ROOT, CFG.get("halflives", "data/processed/mrna_halflife_human.tsv"))
hl = pd.read_csv(hl_path, sep="\t").drop_duplicates("gene").set_index("gene").halflife_h   # columns gene, halflife_h
ligs = {g for x in lr.ligand for g in subunits(x)}
recs = {g for x in lr.receptor for g in subunits(x)}
role = lambda g: ",".join(r for r, s in (("TF", tfs), ("ligand", ligs), ("receptor", recs), ("target", targets)) if g in s)
pd.DataFrame(dict(
    gene=genes, role=[role(g) for g in genes],
    modules=[";".join(f"{m}:{k}" for m, d in modules.items() for k in ("pro", "anti") if g in d[k]) for g in genes],
    halflife_h=[hl.get(g, hl.median()) for g in genes],
)).to_csv(f"{OUT}/genes.tsv", sep="\t", index=False)
E.to_csv(f"{OUT}/edges.tsv", sep="\t", index=False)
lr.drop(columns="genes").to_csv(f"{OUT}/lr_pairs.tsv", sep="\t", index=False)
D.to_csv(f"{OUT}/lr_downstream.tsv", sep="\t", index=False)
say(f"genes: {len(genes)} ({len(tfs)} TFs, {len(ligs)} ligands, {len(recs)} receptor subunits, {len(targets)} module genes; "
    f"roles overlap) -> {OUT}/")
open(f"{OUT}/report.md", "w").write("\n".join(log) + "\n")
