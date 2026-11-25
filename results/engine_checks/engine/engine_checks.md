# Gene engine checks: engine.yaml

3000 cells, 3 types (tumor, tcell, macrophage; baseline: reference), randomly mixed, 5.9 neighbours per cell, 480 steps of 6.0 min = 48 h. 147 genes, 50 TFs, 50 L-R pairs.

## 1. Stability
- NaN or negative values: 0
- max level / max baseline level: 1.62
- per-type mean change over the second half, |log2 fold| (+1 molecule): median 0.014, 95th pct 0.083, max 0.34

## 2. Calibration: simulated mean counts vs the baseline the engine was given

| type | Pearson r (log) | median sim/baseline | genes within 2-fold |
|---|---|---|---|
| tumor | 0.995 | 0.95 | 100% |
| tcell | 0.985 | 1.04 | 99% |
| macrophage | 0.978 | 1.11 | 99% |

## 3. Feedback
- TF->TF edges 294, self-loops 7, largest feedback cluster 48 of 50 TFs; the run above is stable with it

## 4. Signalling: all contacts cut after 24 h vs the intact run
- downstream TF links tested: 150 (0 repressive); moved as the link's sign predicts in 96%; median signed log2 fold -0.31

## 5. Velocity
- corr(s(t+dt) - s(t), velocity x dt): 0.988

## 6. Speed (OPENBLAS_NUM_THREADS=None)
- 3000 cells: 21.7 ms/step -> 3 days = 16 s
- 30000 cells: 335.9 ms/step -> 3 days = 242 s
