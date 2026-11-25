# Behaviour-rule checks (3 types: tumor, tcell, macrophage; 9 enabled rules)

## 1. Rule formula

| case | computed | expected |
|---|---|---|
| no signal | 1.000e-03 | 1.000e-03 |
| up rule at half_max | 3.000e-03 | 3.000e-03 |
| up rule saturated | 5.000e-03 | 5.000e-03 |
| down rule saturated | 1.000e-15 | 0.000e+00 |
| up half + down half | 1.500e-03 | 1.500e-03 |

## 2. Forced apoptosis module (tumor, +2 log2 on every pro gene in half the cells; 100 independent draws of one 6-min step)

| group | mean rule rate (per min) | deaths | expected | z |
|---|---|---|---|---|
| forced | 4.91e-03 | 1473 | 1452 | +0.55 |
| unforced (same type) | 1.82e-04 | 54 | 54 | -0.07 |

## 2b. Chemotaxis
- gradient estimate on a linear CXCL12 field along +x: mean cos(angle to +x) = 0.990 (1 = perfect), radius 40 um
- isolated cells, speed 2 um/min, bias 1: drift along the bias direction 2.00 um/min
- isolated cells, speed 2 um/min, bias 0.5: drift along the bias direction 1.30 um/min
- isolated cells, speed 2 um/min, bias 0: drift along the bias direction 0.01 um/min

## 3. Pressure rules vs none, 36 h

- no rules: cells 3003 -> 6563; median relative pressure at end 2.64
- pressure rules only: cells 3001 -> 4849; median relative pressure at end 1.10

## 4. Gene rules (12 h, all rules vs all rules except gene:<GENE>)

- CXCR4 (log2 vs baseline) in the top 10% migration cells — with gene rules: -0.07; without: -0.11
- CCR7 (log2 vs baseline) in the top 10% migration cells — with gene rules: -0.03; without: -0.05
- CXCL12 (log2 vs baseline) in the top 10% neighbours cells — with gene rules: +0.11; without: -0.03

## 5. Bad config is rejected with:
```
rules.yaml problems:
  rule 1 {'cell_type': 'tumour', 'signal': 'apoptsis', 'behaviour': 'divide', 'direction': 'up', 'halfmax': 2}:
      - unknown key(s) ['halfmax']; allowed: ['behaviour', 'cell_type', 'direction', 'enabled', 'half_max', 'hill_power', 'max', 'min', 'name', 'signal']
      - direction must be increases (with max) or decreases (with min)
      - missing half_max
      - unknown cell_type ['tumour']; types are ['tumor', 'tcell', 'macrophage'] or '*'
      - unknown signal 'apoptsis'; use a module ['apoptosis', 'chemotaxis', 'migration', 'necrosis', 'proliferation'], a network gene, or ['contact_<...>', 'neighbours', 'pressure']
      - unknown behaviour 'divide'; use a tissue parameter ['adhesion', 'apoptosis_duration_min', 'apoptosis_rate', 'division_rate', 'max_adhesion_distance', 'persistence_min', 'radius', 'repulsion', 'speed'] or ['chemotaxis:<...>', 'gene:<...>'] (gene: / chemotaxis: need a network gene)
```

## 6. All enabled rules, 48 h

| type | cells start -> end | division rate (per h), first -> last | death rate (per h), first -> last | speed (um/min), first -> last |
|---|---|---|---|---|
| tumor | 1000 -> 3277 | 0.042 -> 0.036 | 0.0106 -> 0.0152 | 0.17 -> 0.49 |
| tcell | 1000 -> 1785 | 0.018 -> 0.025 | 0.0091 -> 0.0125 | 2.00 -> 2.00 |
| macrophage | 1000 -> 1146 | 0.009 -> 0.012 | 0.0076 -> 0.0085 | 0.59 -> 0.54 |

Per-rule response: mean over steps of corr(signal, behaviour value) in the cells the rule applies to. Expected sign: + for increases, - for decreases (other rules on the same behaviour dilute it).

| rule | cells | corr | sign ok |
|---|---|---|---|
| apoptosis -> apoptosis_rate (increases) | all | +0.62 | yes |
| necrosis -> apoptosis_rate (increases) | all | +0.04 | yes |
| proliferation -> division_rate (increases) | all | +0.57 | yes |
| migration -> speed (increases) | all | +0.10 | yes |
| pressure -> division_rate (decreases) | all | -0.56 | yes |
| migration -> gene:CXCR4 (increases) | all | +0.86 | yes |
| migration -> gene:CCR7 (increases) | all | +0.86 | yes |
| neighbours -> gene:CXCL12 (increases) | all | +0.94 | yes |
| CXCR4 -> chemotaxis:CXCL12 (increases) | tcell | +0.78 | yes |

- chemotaxis: on average 494 cells per step with bias > 0.2; mean cos(angle between their step displacement and the gradient) +0.33 (0 = random walk)
- NaN: 0; rows in sync (asserted every step); cost 245 ms per step at the end-size tissue (~2.9 min per 3 days)

## 7. Extensibility (no code changes)

- tissue-parameter rule (contact_tcell -> adhesion decreases, tumor only): per-cell adhesion read by the mechanics 0.140 in tumor cells touching tcell vs 0.400 in the others (type value 0.4)
- plugin (plugins/examples.py): new signal distance_from_centre -> new behaviour module:apoptosis for 12 h: corr(distance from centre, apoptosis score) +0.56 with the rule vs +0.05 without; new signal dying_neighbours -> speed also loaded and ran
