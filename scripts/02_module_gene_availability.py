"""Phase 0a: which behaviour-module genes can the simulator actually use?

A gene is usable as a GRN target only if it is in the Atera breast panel AND SpaceTravLR trained a
model for it (a <gene>_betadata.feather exists), because that is where the signed, cell-type-
specific weights come from. TFs must be columns of the SpaceTravLR base GRN; ligands/receptors must
appear in CellChat (receptor complexes such as TGFBR1_TGFBR2 count via their subunits).

Inputs (read-only):
  SpaceTravLR base GRN + CellChat   Consensus_signaling/SpaceTravLR/data/
  Atera breast reference            Consensus_signaling/data/atera_breast/atera_breast_subset_hvg3500_biomarkers.h5ad
  Atera breast SpaceTravLR run      Consensus_signaling/atera_breast_subset_hvg3500_biomarkers_10xcelltype_2026-09-25/
  Gene sets (offline)               ~/.cache/gseapy/Enrichr.{MSigDB_Hallmark_2020,GO_Biological_Process_2023}.gmt
Output: data/processed/module_gene_availability.tsv (+ printed summary)
"""
import os, re
import pandas as pd

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
pd.options.future.infer_string = False
import anndata as ad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CS = "/vast/djishnu/Palak/Collabs/Consensus_signaling"
RUN = f"{CS}/atera_breast_subset_hvg3500_biomarkers_10xcelltype_2026-09-25"
REF = f"{CS}/data/atera_breast/atera_breast_subset_hvg3500_biomarkers.h5ad"
GMT = os.path.expanduser("~/.cache/gseapy")

CORE = {  # from the project brief: pro (raises the behaviour) / anti (lowers it)
    "apoptosis": (["CASP3", "CASP7", "CASP8", "CASP9", "BAX", "BAK1", "BID", "BBC3", "PMAIP1", "APAF1", "CYCS", "FAS"],
                  ["BCL2", "BCL2L1", "MCL1", "XIAP", "BIRC5"]),
    "proliferation": (["CCND1", "CCNE1", "E2F1", "MYC", "CDK4", "CDK2", "PCNA", "MKI67", "TOP2A", "CCNB1", "CDK1", "MCM2"],
                      ["CDKN1A", "CDKN1B", "CDKN2A", "RB1", "TP53"]),
    "migration": (["VIM", "SNAI1", "SNAI2", "ZEB1", "TWIST1", "CDH2", "FN1", "MMP2", "MMP9", "ITGB1", "RAC1", "CDC42"],
                  ["CDH1", "EPCAM", "CLDN4"]),
    "chemotaxis": (["CXCR4", "CCR7", "ITGAL"], []),
    "necrosis": (["RIPK1", "RIPK3", "MLKL"], ["CASP8"]),
}
REFSETS = {  # (gmt file, term regex)
    "apoptosis": [("MSigDB_Hallmark_2020", r"^Apoptosis$")],
    "proliferation": [("MSigDB_Hallmark_2020", r"^G2-M Checkpoint$"), ("MSigDB_Hallmark_2020", r"^E2F Targets$")],
    "migration": [("MSigDB_Hallmark_2020", r"^Epithelial Mesenchymal Transition$"),
                  ("GO_Biological_Process_2023", r"\(GO:0016477\)$")],
    "chemotaxis": [("GO_Biological_Process_2023", r"\(GO:0006935\)$")],
    "necrosis": [("GO_Biological_Process_2023", r"\(GO:0070266\)$")],
}


def read_gmt(name):
    sets = {}
    with open(f"{GMT}/Enrichr.{name}.gmt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            sets[p[0]] = [g for g in p[2:] if g]
    return sets


panel = set(ad.read_h5ad(REF, backed="r").var["gene_name"].astype(str))
modelled = {f[:-len("_betadata.feather")] for f in os.listdir(RUN) if f.endswith("_betadata.feather")}
grn = pd.read_parquet(f"{CS}/SpaceTravLR/data/human_base_grn.parquet")
tfs = set(grn.columns) - {"peak_id", "gene_short_name"}
cc = pd.read_csv(f"{CS}/SpaceTravLR/data/cellchat_human.csv")
ligands = set(cc["ligand"])
receptor_subunits = {s for r in cc["receptor"] for s in str(r).split("_")}
print(f"panel genes {len(panel)} | modelled targets (betadata) {len(modelled)} | base-GRN TFs {len(tfs)} | "
      f"CellChat ligands {len(ligands)} / receptor subunits {len(receptor_subunits)}")
print(f"in panel: TFs {len(tfs & panel)}, ligands {len(ligands & panel)}, receptor subunits {len(receptor_subunits & panel)}\n")

gmts = {n: read_gmt(n) for n in {g for v in REFSETS.values() for g, _ in v}}
rows = []
for mod, (pro, anti) in CORE.items():
    refgenes, terms = set(), []
    for gname, rx in REFSETS[mod]:
        hit = [t for t in gmts[gname] if re.search(rx, t)]
        terms += hit
        for t in hit:
            refgenes |= set(gmts[gname][t])
    for role, genes in (("pro", pro), ("anti", anti)):
        for g in genes:
            rows.append(dict(module=mod, role=role, gene=g, source="core", in_panel=g in panel,
                             modelled=g in modelled, is_tf=g in tfs))
    for g in sorted(refgenes - set(pro) - set(anti)):
        rows.append(dict(module=mod, role="reference_set", gene=g, source="; ".join(terms), in_panel=g in panel,
                         modelled=g in modelled, is_tf=g in tfs))
    core = pd.DataFrame([r for r in rows if r["module"] == mod and r["source"] == "core"])
    usable_ref = sum(1 for r in rows if r["module"] == mod and r["role"] == "reference_set" and r["modelled"])
    print(f"{mod:13s} core usable {int(core.modelled.sum()):2d}/{len(core):2d}  "
          f"(pro {int(core[core.role=='pro'].modelled.sum())}/{len(pro)}, anti {int(core[core.role=='anti'].modelled.sum())}/{len(anti)})  "
          f"| reference-set genes usable {usable_ref} of {len(refgenes)}  terms: {terms}")
    missing = core.loc[~core.modelled, "gene"].tolist()
    if missing:
        print(f"{'':13s} not usable: {', '.join(missing)}")

out = pd.DataFrame(rows)
os.makedirs(f"{ROOT}/data/processed", exist_ok=True)
out.to_csv(f"{ROOT}/data/processed/module_gene_availability.tsv", sep="\t", index=False)
print(f"\nwrote data/processed/module_gene_availability.tsv ({len(out)} rows)")
