# Gene engine checks

3000 cells (1000 per type, randomly mixed), 5.9 neighbours per cell, 480 steps of 6.0 min = 48 h
147 genes, 50 TFs, 50 L-R pairs

## 1. Stability
- NaN or negative values: 0
- max level / max reference level: 1.58
- per-type mean drift over day 2, |log fold|: median 0.035, 95th pct 0.163

## 2. Calibration against the reference (per type, across genes)
Calibrated on the starting layout: each type's average ligand input is part of its baseline.

| run | type | Pearson r of log mean counts | median sim/ref ratio | zero fraction sim | zero fraction ref |
|---|---|---|---|---|---|
| 48 h, signalling on | tumor | 0.996 | 0.94 | 0.71 | 0.74 |
| 48 h, signalling on | tcell | 0.987 | 1.04 | 0.86 | 0.88 |
| 48 h, signalling on | macrophage | 0.979 | 1.13 | 0.88 | 0.91 |

## 3. Feedback
- TF->TF edges 294, self-loops 7, largest feedback cluster (strongly connected) 48 of 50 TFs — run above is stable with it

## 4. Signalling: all contacts removed after 24 h, compared with the intact run
- downstream TF links tested (strongest-receiving quarter of cells): 150; went DOWN in 97%; median log2 fold -0.44; largest drop {'pair': 'OCLN->OCLN', 'tf': 'FOXA1', 'log2fc': -4.406780659921636}

## 5. Velocity
- corr(s(t+dt) - s(t), velocity x dt) over all cells x genes: 0.988 (not 1: velocity is instantaneous at t+dt, and promoters switch within the step)

## 6. Speed (1 thread: OPENBLAS_NUM_THREADS=1)
- 3000 cells: 23.6 ms/step -> 3 days (720 steps) = 17 s
- 30000 cells: 350.1 ms/step -> 3 days (720 steps) = 252 s
