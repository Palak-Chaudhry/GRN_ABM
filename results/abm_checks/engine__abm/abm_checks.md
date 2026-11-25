# ABM checks: abm.yaml (3 types: tumor, tcell, macrophage)

## 1. Two-cell equilibrium (no motility)

| type | analytic distance (um) | simulated | error |
|---|---|---|---|
| tumor | 16.019 | 16.019 | 0.00% |
| tcell | 9.524 | 9.524 | 0.00% |
| macrophage | 19.048 | 19.048 | 0.00% |

## 2. Motility: mean squared displacement of isolated cells

| type | speed (um/min) | persistence (min) | MSD ratio sim/theory at 10 / 60 / 300 min |
|---|---|---|---|
| tumor | 0 | - | not motile |
| tcell | 2 | 5 | 1.00 / 1.00 / 0.99 |
| macrophage | 0.5 | 10 | 1.05 / 1.04 / 0.94 |

## 3. Event counts (fixed population of 2000, 200 steps of 6 min)

- divisions: 2398 observed vs 2386 expected (z = +0.25)
- deaths: 1189 observed vs 1198 expected (z = -0.27)

## 4. Packed tissue (synthetic layout, mechanics only then with motility)

- 3000 cells; strongly overlapping pairs (d < 0.7 x (R_i + R_j)): start 31.1% -> after 1 h relaxing 0.0% -> after 1 h with motility 5.6%
- median contact distance / (R_i + R_j): 0.84 (two-cell equilibrium for default cells is 0.95)
- NaN positions: 0; signalling neighbours per cell: median 5, mean 4.6, isolated cells 1.1%
- mechanics cost: 111 ms per 6-min phenotype step (60 substeps)

## 5. Reference crop (500 um)
- cells per type: tumor 636, tcell 729, macrophage 31; signalling neighbours per cell: median 6 (only the configured types are kept, so other cell types leave gaps)

## 6. Coupled run: tissue + gene engine, 24 h, base rates

| type | start | end (living) | expected from rates* |
|---|---|---|---|
| tumor | 1000 | 2541 | 2612 |
| tcell | 1000 | 1241 | 1235 |
| macrophage | 1000 | 929 | 926 |
*N0 exp((division - death) t); dying cells are excluded from 'living', so 'end' sits slightly below.
- rows in sync throughout (asserted every step); events logged: 2401; cost 178 ms per 6-min step -> 3 days of this tissue ~ 2.1 min
