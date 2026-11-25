# GRN_ABM vs scMultiSim on the same inputs

- GRN_ABM: /vast/djishnu/Palak/GRN_ABM/results/runs/main_72h_v2 (3000 cells at t = 0, 8846 at 72 h)
- scMultiSim spatial: /vast/djishnu/Palak/GRN_ABM/results/_tests/sms_sp (300 cells, populations ['1', '2', '3'])
- shared genes (network genes present in both simulators and the reference): 147 of 147
- scMultiSim populations matched to reference types by expression: 1 -> tumor (r = -0.05), 2 -> tcell (r = -0.02), 3 -> macrophage (r = -0.01)

## 1. Realism vs the reference (observed counts, shared genes, depth-matched)

Mean over cell types (scMultiSim types = its populations matched to the reference):

| simulator | mean expression r vs reference (↑) | zero-fraction error (↓) | dispersion error, log10 (↓) | gene-gene correlation structure r (↑) |
|---|---|---|---|---|
| GRN_ABM 72 h | 0.97 | 0.040 | 0.12 | 0.35 |
| GRN_ABM t=0 | 0.99 | 0.024 | 0.13 | 0.35 |
| scMultiSim | -0.03 | 0.187 | 0.30 | 0.06 |

## 2. Embeddings and space

- GRN_ABM 72 h: neighbours per cell 5.1; same-type contact enrichment z tumor +48, tcell -15, macrophage +15
- scMultiSim: neighbours per cell 3.5; same-type contact enrichment z tumor +20, tcell +6, macrophage +2

## 3. GRN fidelity: are true edges visible in the data? (within each type)

| simulator | counts | AUROC edges vs non-edges, within type (0.5 = invisible) | AUROC, types pooled | sign agreement of correlated edges | median abs rho, edges / non-edges |
|---|---|---|---|---|---|
| GRN_ABM t=0 | true | 0.593 | 0.539 | 0.82 | 0.220 / 0.155 |
| GRN_ABM t=0 | obs | 0.512 | 0.526 | 0.76 | 0.026 / 0.025 |
| GRN_ABM 72 h | true | 0.546 | 0.528 | 0.74 | 0.295 / 0.257 |
| GRN_ABM 72 h | obs | 0.522 | 0.521 | 0.81 | 0.021 / 0.020 |
| scMultiSim | true | 0.496 | 0.509 | 0.56 | 0.070 / 0.072 |
| scMultiSim | obs | 0.499 | 0.522 | 0.57 | 0.074 / 0.074 |

## 4. Cell-cell signalling: does a neighbour's ligand show up in the receiver?

| simulator | receiver readout | AUROC cognate vs non-cognate (0.5 = no specificity) | median corr cognate / non-cognate | cognate pair in its ligand's top 10% |
|---|---|---|---|---|
| GRN_ABM 72 h | receptor | 0.52 | +0.197 / +0.185 | 4% |
| scMultiSim | receptor | 0.64 | +0.046 / +0.012 | 28% |
| GRN_ABM 72 h | downstream TFs | 0.49 | +0.302 / +0.312 | 6% |
| scMultiSim | downstream TFs | 0.54 | +0.021 / +0.008 | 12% |
(GRN_ABM: ligand acts on the receiver's downstream TFs, gated by its receptor; scMultiSim: ligand acts on the receptor gene itself, which then regulates its TFs through the added receptor -> TF edges.)

## 5. Dynamics

- GRN_ABM: the same cell 3 h apart: median per-gene correlation (within type) 0.83 (genes with half-life > 10 h: 0.93; < 3 h: 0.73); 8998 cells followed
- scMultiSim: no cell is followed in time (each simulated cell is an independent sample; the probe found lag-1 within-cell autocorrelation ~0), so this cannot be measured
- scMultiSim dynamic GRN: 20 of the GRN's edges change (rewiring along cell order); target follows its edge weight (sign-adjusted Spearman, within population): median -0.00; follows weight x regulator expression (the same quantity as GRN_ABM's edge activity): median +0.00
- GRN_ABM: edges are fixed; their activity (regulator activation x weight) changes with each cell's state and neighbourhood. Target follows its edge activity over time (per type, 30 most dynamic edges per type, sign-adjusted Spearman): median +0.86

## 6. Velocity, cost and features

- GRN_ABM velocity: true velocity for every cell and snapshot, with CCI, movement and divisions (see the run's 02* plots for scvelo vs truth)
- runtime: GRN_ABM 385 s for 3000 -> 8846 cells over 72 h (720 coupled steps, 25 recorded snapshots); scMultiSim 19 s for 300 cells (one static snapshot)

| feature | GRN_ABM | scMultiSim |
|---|---|---|
| same GRN / L-R / L-R -> TF inputs | yes | yes (receptor -> TF as GRN edges; homophilic pairs dropped) |
| cell types match a reference's expression | yes (calibrated per type) | no (random type profiles) |
| real time; the same cell followed over time | yes | no (cells are independent samples along pseudotime) |
| cells move, divide, die; tissue grows and reorganises | yes | no (static grid) |
| expression changes behaviour (rules) | yes | no |
| per-type GRN activity | yes (TF expression differs by type) | no (one GRN for all types) |
| GRN rewiring (edges created / deleted) | no (fixed edges, activity changes) | yes (dynamic.GRN) |
| spliced / unspliced and velocity with CCI | yes | no (only without CCI) |
| differentiation trajectories | no (fixed types) | yes (tree) |
| ATAC / chromatin | no | yes |

