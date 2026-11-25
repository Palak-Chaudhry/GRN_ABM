"""Phase 0b: per-cell-type expression of every gene in the full Atera breast panel.

Usage: python scripts/03_reference_type_expression.py [--config configs/engine.yaml]
Reads baseline.reference in the engine config: labels_h5ad (obs[label_key] = cell type labels), and raw
counts from counts_h5 (10x .h5) if given, else from labels_h5ad (layers[counts_layer], or X).

Used for (a) expression filters when choosing TFs / ligands / receptors and (b) later, the gene
engine's per-type basal rates. Reads the full 18k-gene Xenium matrix (the SpaceTravLR subset h5ad
only has 3,652 genes) and takes cell types from the subset h5ad (same barcodes).

Output: data/processed/reference_type_expression.tsv, one row per gene, with for each type:
  det_<type>   fraction of cells with >= 1 count
  mean_<type>  mean raw counts per cell
  cpm_<type>   mean counts per 10k total counts (library-size normalised)
"""
import os, time
import numpy as np
import pandas as pd
import yaml
import scipy.sparse as sp

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
pd.options.future.infer_string = False
import anndata as ad
import scanpy as sc

import argparse
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--config", default=f"{ROOT}/configs/engine.yaml")
CFG = yaml.safe_load(open(ap.parse_args().config))        # cell types + reference paths
REF = CFG["baseline"]["reference"]
REF_DIR = REF["dir"]
TYPES = CFG["cell_types"]                                  # short name -> reference annotation label
path = lambda p: os.path.join(REF_DIR, p)

t0 = time.time()
cells_h5ad = ad.read_h5ad(path(REF["labels_h5ad"]), backed="r")
labels = cells_h5ad.obs[REF.get("label_key", "cell_type")].astype(str)
if REF.get("counts_h5"):                                   # a separate full-panel 10x matrix (same barcodes)
    full = sc.read_10x_h5(path(REF["counts_h5"]), gex_only=True)
    full.var_names_make_unique()
else:                                                      # raw counts inside the labels h5ad
    full = ad.read_h5ad(path(REF["labels_h5ad"]))
    layer = REF.get("counts_layer")
    if layer:
        full.X = full.layers[layer]
full.X = sp.csr_matrix(full.X)
print(f"count matrix {full.shape} loaded in {time.time() - t0:.0f}s")

rows = {}
for short, label in TYPES.items():
    cells = labels.index[labels == label].intersection(full.obs_names)
    X = full[cells].X.tocsr()
    lib = np.asarray(X.sum(1)).ravel()
    rows[f"det_{short}"] = np.asarray((X > 0).mean(0)).ravel()
    rows[f"mean_{short}"] = np.asarray(X.mean(0)).ravel()
    rows[f"cpm_{short}"] = np.asarray((X.multiply(1e4 / np.maximum(lib, 1)[:, None])).mean(0)).ravel()
    print(f"{short:10s} {label:32s} cells {len(cells):6d}  median counts/cell {np.median(lib):.0f}")

out = pd.DataFrame(rows, index=full.var_names)
out.index.name = "gene"
out.to_csv(os.path.join(ROOT, REF["table"]), sep="\t", float_format="%.5g")
print(f"wrote {REF['table']} {out.shape} in {time.time() - t0:.0f}s")
