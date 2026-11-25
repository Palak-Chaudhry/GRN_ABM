# ABM checks: abm.yaml (5 types: A, B, C, D, E)

## 1. Two-cell equilibrium (no motility)

| type | analytic distance (um) | simulated | error |
|---|---|---|---|
| A | 16.019 | 16.019 | 0.00% |
| B | 16.019 | 16.019 | 0.00% |
| C | 16.019 | 16.019 | 0.00% |
| D | 16.019 | 16.019 | 0.00% |
| E | 16.019 | 16.019 | 0.00% |

## 2. Motility: mean squared displacement of isolated cells

| type | speed (um/min) | persistence (min) | MSD ratio sim/theory at 10 / 60 / 300 min |
|---|---|---|---|
| A | 0 | - | not motile |
| B | 0 | - | not motile |
| C | 0 | - | not motile |
| D | 0 | - | not motile |
| E | 0 | - | not motile |

## 3. Event counts (fixed population of 2000, 200 steps of 6 min)

- divisions: 2398 observed vs 2386 expected (z = +0.25)
- deaths: 1189 observed vs 1198 expected (z = -0.27)

## 4. Packed tissue (synthetic layout, mechanics only then with motility)

- 5000 cells; strongly overlapping pairs (d < 0.7 x (R_i + R_j)): start 49.9% -> after 1 h relaxing 0.0% -> after 1 h with motility 0.0%
- median contact distance / (R_i + R_j): 0.93 (two-cell equilibrium for default cells is 0.95)
- NaN positions: 0; signalling neighbours per cell: median 5, mean 5.0, isolated cells 0.1%
- mechanics cost: 203 ms per 6-min phenotype step (60 substeps)

## 6. Coupled run: tissue + gene engine, 24 h, base rates

| type | start | end (living) | expected from rates* |
|---|---|---|---|
| A | 1000 | 2530 | 2612 |
| B | 1000 | 2651 | 2612 |
| C | 1000 | 2651 | 2612 |
| D | 1000 | 2518 | 2612 |
| E | 1000 | 2606 | 2612 |
*N0 exp((division - death) t); dying cells are excluded from 'living', so 'end' sits slightly below.
- rows in sync throughout (asserted every step); events logged: 9496; cost 527 ms per 6-min step -> 3 days of this tissue ~ 6.3 min
