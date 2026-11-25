# Gene engine checks: engine_synthetic5.yaml

5000 cells, 5 types (A, B, C, D, E; baseline: synthetic), randomly mixed, 5.8 neighbours per cell, 480 steps of 6.0 min = 48 h. 147 genes, 50 TFs, 50 L-R pairs.

## 1. Stability
- NaN or negative values: 0
- max level / max baseline level: 1.84
- per-type mean change over the second half, |log2 fold| (+1 molecule): median 0.019, 95th pct 0.104, max 0.27

## 2. Calibration: simulated mean counts vs the baseline the engine was given

| type | Pearson r (log) | median sim/baseline | genes within 2-fold |
|---|---|---|---|
| A | 0.997 | 1.11 | 99% |
| B | 0.994 | 1.05 | 99% |
| C | 0.996 | 1.00 | 99% |
| D | 0.997 | 1.06 | 99% |
| E | 0.992 | 1.09 | 99% |

## 3. Feedback
- TF->TF edges 294, self-loops 7, largest feedback cluster 48 of 50 TFs; the run above is stable with it

## 4. Signalling: all contacts cut after 24 h vs the intact run
- downstream TF links tested: 150; went DOWN in 97%; median log2 fold -0.68

## 5. Velocity
- corr(s(t+dt) - s(t), velocity x dt): 0.988

## 6. Speed (OPENBLAS_NUM_THREADS=2)
- 3000 cells: 9.1 ms/step -> 3 days = 7 s
- 30000 cells: 184.3 ms/step -> 3 days = 133 s
