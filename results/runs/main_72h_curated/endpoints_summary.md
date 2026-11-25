# Endpoint export: main_72h_curated

- source: `results/runs/main_72h_curated/snapshots.h5ad` -- 132548 cell-snapshots x 167 genes, 25 snapshots from t = 0 h to t = 72 h

## t = 0 h -> `adata_t0.h5ad`
- 3000 cells x 167 genes; tumor 1000, tcell 1000, macrophage 1000
- spatial extent: x [-533, 524] um, y [-502, 487] um
- mean observed counts per cell 60.7; dying 0; mean neighbours 4.5
- layers: S_true, U_true, spliced, unspliced | obsm: lr_received, spatial | obs: 28 columns

## t = 72 h -> `adata_tfinal.h5ad`
- 7970 cells x 167 genes; tumor 4749, tcell 1999, macrophage 1222
- spatial extent: x [-1128, 1275] um, y [-1336, 1370] um
- mean observed counts per cell 79.4; dying 745; mean neighbours 4.4
- layers: S_true, U_true, spliced, unspliced | obsm: lr_received, spatial | obs: 28 columns

