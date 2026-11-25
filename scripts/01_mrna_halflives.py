"""Per-gene human mRNA half-lives (hours) for the gene engine's decay rates.

Source: Agarwal & Kelley 2022, Genome Biology 23:245 (doi 10.1186/s13059-022-02811-x),
github.com/vagarwal87/saluki_paper. Files in data/external/saluki_halflives/:
  all_HLs_human_PC1.txt   consensus stability score per gene = PC1 over 39 human datasets
                          (unitless; ranks genes robustly but has no time scale)
  all_HLs_human.txt.gz    the per-study values, each transformed by the authors
                          (Fig1_S1_S2/master_comparison_human.R)

Absolute scale: the authors' script stores the Dieterich (4sU, HEK293 + MCF7) and Mortazavi (4sU,
GM12878) studies as log10(x + 1). Back-transformed, only Mortazavi GM12878 gives a plausible scale
in hours (median ~4.7 h, in line with the several-hour medians of metabolic-labelling studies).
Dieterich back-transforms to medians of 86-507 "hours", i.e. not hours (most likely minutes), so it
is reported for transparency but NOT used. GM12878 is a lymphoblastoid (immune) line.
Caveat: quantile mapping onto one study keeps the consensus ranking but inherits that study's
range, so the fastest transcripts (FOS, MYC: ~15-30 min in the literature) come out longer.

Method: each gene's consensus PC1 quantile is mapped onto the pooled distribution of absolute
half-lives from those calibration studies (quantile mapping). Genes without a PC1 value get the
pooled median and are flagged.

Output: data/processed/mrna_halflife_human.tsv
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXT = os.path.join(ROOT, "data", "external", "saluki_halflives")
OUT = os.path.join(ROOT, "data", "processed", "mrna_halflife_human.tsv")
CHECKED = ["Dieterich_4sU_MCF7_1", "Dieterich_4sU_MCF7_2", "Dieterich_4sU_HEK293_1",
           "Dieterich_4sU_HEK293_2", "Mortazavi_4sU_GM12878_1"]   # log10(x + 1) per the authors' script
CALIB = ["Mortazavi_4sU_GM12878_1"]                                  # the only one whose x is in hours

pc1 = pd.read_csv(os.path.join(EXT, "all_HLs_human_PC1.txt"), sep="\t")
per = pd.read_csv(os.path.join(EXT, "all_HLs_human.txt.gz"), sep="\t")

hours = {c: 10 ** per[c].dropna() - 1 for c in CHECKED}
print("back-transformed studies, median value, n genes (only CALIB used):")
for c, h in hours.items():
    print(f"  {c:26s} {h.median():6.2f}  n={len(h)}")

# PC1 sign: positive correlation with the calibration studies means higher PC1 = more stable
m = pc1.merge(per[["ensembl_gene_id"] + CHECKED], on="ensembl_gene_id", how="left")
rho = {c: m["human_PC1"].corr(m[c], method="spearman") for c in CHECKED}
print("Spearman(PC1, study):", {k: round(v, 2) for k, v in rho.items()})
sign = 1.0 if np.nanmean(list(rho.values())) > 0 else -1.0

pooled = np.sort(np.concatenate([hours[c].values for c in CALIB]))
pooled = pooled[pooled > 0]
q = (sign * pc1["human_PC1"]).rank(pct=True)
pc1["halflife_h"] = np.quantile(pooled, q.clip(1e-6, 1 - 1e-6))
pc1["source"] = "consensus_PC1_quantile_mapped"
pc1 = pc1.rename(columns={"external_gene_name": "gene", "human_PC1": "consensus_pc1"})
pc1 = pc1.dropna(subset=["gene"]).drop_duplicates("gene")
pc1.attrs["pooled_median_h"] = float(np.median(pooled))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
pc1[["gene", "ensembl_gene_id", "consensus_pc1", "halflife_h", "source"]].to_csv(OUT, sep="\t", index=False)
print(f"\nwrote {OUT}: {len(pc1)} genes; half-life median {pc1.halflife_h.median():.2f} h, "
      f"5-95% {pc1.halflife_h.quantile(.05):.2f}-{pc1.halflife_h.quantile(.95):.2f} h; "
      f"pooled calibration median {np.median(pooled):.2f} h (default for genes not in the table)")
