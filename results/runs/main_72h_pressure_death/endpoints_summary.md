# Endpoint export: main_72h_pressure_death

- source: `results/runs/main_72h_pressure_death/snapshots.h5ad` -- 116129 cell-snapshots x 167 genes, 25 snapshots from t = 0 h to t = 72 h

## t = 0 h -> `adata_t0.h5ad`
- 3000 cells x 167 genes; tumor 1000, tcell 1000, macrophage 1000
- spatial extent: x [-533, 524] um, y [-502, 487] um
- mean observed counts per cell 60.7; dying 0; mean neighbours 4.5
- layers: S_true, U_true, spliced, unspliced | obsm: lr_received, spatial | obs: 28 columns

## t = 72 h -> `adata_tfinal.h5ad`
- 6653 cells x 167 genes; tumor 4323, tcell 1588, macrophage 742
- spatial extent: x [-1251, 1271] um, y [-1148, 1108] um
- mean observed counts per cell 75.6; dying 1002; mean neighbours 3.3
- layers: S_true, U_true, spliced, unspliced | obsm: lr_received, spatial | obs: 28 columns

