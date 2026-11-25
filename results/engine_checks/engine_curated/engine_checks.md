# Gene engine checks: engine_curated.yaml

3000 cells, 3 types (tumor, tcell, macrophage; baseline: reference), randomly mixed, 5.9 neighbours per cell, 480 steps of 6.0 min = 48 h. 167 genes, 69 TFs, 50 L-R pairs.

## 1. Stability
- NaN or negative values: 0
- max level / max baseline level: 1.53
- per-type mean change over the second half, |log2 fold| (+1 molecule): median 0.015, 95th pct 0.084, max 0.36

## 2. Calibration: simulated mean counts vs the baseline the engine was given

| type | Pearson r (log) | median sim/baseline | genes within 2-fold |
|---|---|---|---|
| tumor | 0.994 | 0.93 | 100% |
| tcell | 0.989 | 1.04 | 100% |
| macrophage | 0.983 | 1.15 | 98% |

## 3. Feedback
- TF->TF edges 406, self-loops 6, largest feedback cluster 64 of 69 TFs; the run above is stable with it

## 4. Signalling: all contacts cut after 24 h vs the intact run
- downstream TF links tested: 171 (22 repressive); moved as the link's sign predicts in 92%; median signed log2 fold -0.64
  - activating  n=149  as expected 99%  median log2 fold -0.68
  - repressive  n= 22  as expected 45%  median log2 fold -0.08
  - TFs whose L-R input is all one sign (17 of 19): as expected 98% (n=139)
  - a TF fed by both signs, or driven mostly by the GRN rather than by L-R, is not a clean test here: the all-contacts cut removes its activating and repressive input at the same time
  - per TF with repressive input: LEF1 100% (n=3), MYC 27% (n=11), NFKB1 0% (n=1), SNAI1 0% (n=3), TCF7 100% (n=4)

## 5. Velocity
- corr(s(t+dt) - s(t), velocity x dt): 0.988

## 6. Speed (OPENBLAS_NUM_THREADS=None)
- 3000 cells: 25.4 ms/step -> 3 days = 18 s
- 30000 cells: 393.6 ms/step -> 3 days = 283 s
