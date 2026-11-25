"""Expand a curated pathway -> TF table into the ligand,receptor,tf,sign,weight table that
`network.yaml: lr.downstream` expects.

Replaces the random receptor->TF assignment in scripts/04 (which drew `lr_downstream_tfs` uniform
TFs per pair) with canonical pathway effectors, so what a signal DOES is curated rather than seeded
noise. Pairs whose pathway is not in the curation table fall back to the immediate-early set and are
listed in the report, so the gap is visible instead of silent.

  python scripts/04b_build_lr_downstream.py [--config configs/network.yaml]

Outputs: data/processed/lr_downstream_curated.tsv (+ .report.md)
"""
import argparse
import os
import pandas as pd
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--config", default=f"{ROOT}/configs/network.yaml")
ap.add_argument("--pathway-tfs", default=f"{ROOT}/data/external/pathway_tfs.tsv")
ap.add_argument("--out", default=f"{ROOT}/data/processed/lr_downstream_curated.tsv")
ap.add_argument("--fallback", default="FOS,JUN,EGR1",
                help="TFs used for pairs whose pathway is not curated (immediate-early response)")
ARGS = ap.parse_args()
CFG = yaml.safe_load(open(ARGS.config))
log = []


def say(*a):
    log.append(" ".join(map(str, a))); print(log[-1])


def read_table(path):
    return pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path, sep=None, engine="python")


# ---- inputs
lr = read_table(CFG["lr"]["path"])
if "pathway" not in lr.columns:
    raise SystemExit(f"{CFG['lr']['path']} has no `pathway` column; cannot map pathways to TFs")
lr = lr.drop_duplicates(["ligand", "receptor"])
P = pd.read_csv(ARGS.pathway_tfs, sep="\t", comment="#")
say(f"L-R list: {len(lr)} pairs, {lr.pathway.nunique()} pathways ({CFG['lr']['path']})")
say(f"curation: {len(P)} pathway->TF links over {P.pathway.nunique()} pathways, {P.tf.nunique()} distinct TFs "
    f"({os.path.relpath(ARGS.pathway_tfs, ROOT)})")

# ---- expand: every pair inherits its pathway's effectors
D = lr.merge(P, on="pathway", how="left")
missing = sorted(set(lr.pathway) - set(P.pathway))
hit = D.tf.notna()
say(f"pairs with a curated pathway: {D.loc[hit, ['ligand', 'receptor']].drop_duplicates().shape[0]}/{len(lr)}")

# ---- fallback for uncurated pathways
fb = [t.strip() for t in ARGS.fallback.split(",") if t.strip()]
gap = lr[lr.pathway.isin(missing)]
if len(gap) and fb:
    F = pd.DataFrame([dict(ligand=r.ligand, receptor=r.receptor, pathway=r.pathway, tf=t, sign=1, weight=0.5)
                      for r in gap.itertuples() for t in fb])
    D = pd.concat([D[hit], F], ignore_index=True)
    say(f"fallback ({'+'.join(fb)}, weight 0.5) applied to {len(gap)} pairs in {len(missing)} uncurated pathways")
else:
    D = D[hit]
if missing:
    say(f"uncurated pathways ({len(missing)}): {', '.join(missing[:40])}{' ...' if len(missing) > 40 else ''}")

D = (D[["ligand", "receptor", "pathway", "tf", "sign", "weight"]]
     .astype({"sign": int}).drop_duplicates().sort_values(["ligand", "receptor", "tf"]))
say(f"L-R -> TF links: {len(D)} ({(D.sign < 0).sum()} inhibitory, {(D.sign > 0).sum()} activating)")
say(f"median TFs per pair: {D.groupby(['ligand', 'receptor']).size().median():.0f}")

os.makedirs(os.path.dirname(ARGS.out), exist_ok=True)
D.to_csv(ARGS.out, sep="\t", index=False)
say(f"-> {os.path.relpath(ARGS.out, ROOT)}")
open(ARGS.out.replace(".tsv", ".report.md"), "w").write("\n".join(log) + "\n")
