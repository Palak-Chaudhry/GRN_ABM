# GRN_ABM user guide

GRN_ABM simulates a 2D tissue in which every cell carries a gene regulatory network (GRN). Genes and
behaviour drive each other:

- **genes drive behaviour:** expression sets division, death, speed and chemotaxis;
- **the tissue drives genes:** crowding, contacts and ligands from neighbours change expression.

The output looks like spatial single-cell data:

- spliced and unspliced counts with x, y positions over time;
- every ground truth behind them: true mRNA, true RNA velocity, GRN and signalling activity, rule
  states, cell tracks, divisions and deaths.

This guide covers:

- what you need;
- how to run it;
- every configuration parameter;
- where the results go;
- what the model does and does not do.

Units throughout are micrometres (µm), minutes and mRNA molecules. Rates are per minute unless the
name says otherwise.

---

## 1. How the model works

```
 start-up (once)                                    one phenotype step (dt_min = 6 min), repeated
 ───────────────                                    ──────────────────────────────────────────────
 place cells (synthetic disc or reference crop)     1. RULES read each cell's state (gene modules, genes,
 relax + move for warmup_h                             pressure, neighbours, contacts) and set its
 calibrate genes to each type's baseline               behaviour: division / death rates, speed, adhesion,
 draw t = 0 expression (steady state, NB counts)       chemotaxis, transcription changes, ...
 48 h gene burn-in on the still tissue              2. TISSUE: 60 mechanics + motility substeps (0.1 min)
   (virtual divisions, re-centred to baseline)      3. GENES: GRN + ligand-receptor input from contacting
                                                       neighbours + rule input -> promoter bursts -> u, s
                                                    4. EVENTS: divisions (mRNA split between daughters)
                                                       and deaths (dying cells stop moving and signalling,
                                                       removed after apoptosis_duration_min)
```

**Gene engine.** For each cell and gene, a promoter switches between ON and OFF. The ON rate is:

```
k_on = k_on0[type, gene] · exp( GRN input + L-R input + rule input )
```

- **GRN input** is a sum over the gene's regulators of `sign · weight · x^n/(x^n + h^n)`, where `x`
  is the regulator's mRNA. This is SERGIO-style Hill regulation.
- **L-R input** comes from each L-R pair. Its signal is the summed ligand of the cell's contacting
  neighbours, gated by the cell's own receptor level. The pair's downstream TFs receive that signal,
  multiplied by the link's sign: a repressive link (e.g. cadherin contact -> TCF7/LEF1, contact
  inhibition) lowers its TF instead of raising it.
- **mRNA kinetics.** While the promoter is ON, unspliced mRNA `u` is made at rate `α`. It is spliced
  into `s` at rate `β` (15 min half-time) and `s` decays at rate `γ`, which comes from literature
  mRNA half-lives.
- **Velocity** is `β·u − γ·s`.
- **Observed counts** are `Poisson(capture · s)`.

**Calibration.** `k_on0` is solved so that each cell type, in its starting neighbourhood, reproduces
its baseline expression. The baseline comes from a reference dataset or a synthetic generator. After
that, only **changes** in neighbourhood, GRN state and rules move expression.

**Tissue.**
- Cells are discs with PhysiCell-form repulsion and adhesion.
- Movement is a persistent random walk, optionally biased up a ligand gradient.
- Division and death are Poisson events.

**Rules.** Rules follow the PhysiCell grammar: *signal → Hill response → behaviour up or down*. Any
signal can drive any behaviour, and several rules on the same behaviour combine. This is where
mechanics, genes and behaviour connect, for example crowding → gene → death, or receptor → migration
up a gradient.

---

## 2. Requirements

| what | where / version |
|---|---|
| Python env for simulation and most plots | `/ihome/djishnu/pac271/.conda/envs/spatial` (numpy, scipy, pandas, anndata 0.12, scanpy, leidenalg, matplotlib, imageio, pyyaml) |
| Python env for RNA-velocity plots | `/ihome/djishnu/pac271/.conda/envs/velo_py` (scvelo 0.3.4) |
| Cluster | Slurm; `scripts/run_py.sbatch` and `scripts/run_analysis.sbatch` use partition `htc` |
| Environment variables | `HDF5_USE_FILE_LOCKING=FALSE` (set by the sbatch scripts) |

All commands below are run from the project folder, `/vast/djishnu/Palak/GRN_ABM`. Every script also
accepts absolute paths.

---

## 3. Quick start

The default inputs are already built: a 147-gene network and Atera breast baselines for tumour, T cell
and macrophage. A 72 h run with 3,000 cells takes about 5 minutes.

```bash
cd /vast/djishnu/Palak/GRN_ABM
sbatch -J run      scripts/run_py.sbatch scripts/08_run_simulation.py --hours 72 --name my_first_run
# when it has finished:
sbatch -J analysis scripts/run_analysis.sbatch results/runs/my_first_run
```

Results:
- data: `results/runs/my_first_run/`
- plots: `results/runs/my_first_run/plots/`
- animation: `results/runs/my_first_run/movement.gif`

From Python:

```python
import sys; sys.path.insert(0, "src")
from grn_abm import load_simulation
from grn_abm.record import Recorder

sim = load_simulation("configs/engine.yaml", "configs/abm.yaml", "configs/rules.yaml", root=".", seed=0)
rec = Recorder(sim, snapshot_min=180)          # optional: records everything for the analysis scripts
sim.run(hours=72, callback=rec)                # or: for _ in range(n): sim.step()
rec.save("results/runs/my_run")
print(sim.counts())                            # living cells per type
```

---

## 4. Workflow: the scripts, their inputs and outputs

| step | script | input | output | when to (re)run |
|---|---|---|---|---|
| 01 | `01_mrna_halflives.py` | `data/external/saluki_halflives/` (Agarwal & Kelley 2022) | `data/processed/mrna_halflife_human.tsv` | once (human); replace the table for another species |
| 02 | `02_module_gene_availability.py` | SpaceTravLR run + Atera panel | `data/processed/module_gene_availability.tsv` | optional diagnostic for the default data |
| 03 | `03_reference_type_expression.py --config <engine.yaml>` | reference h5ad (+ optional 10x .h5) | `baseline.reference.table` (per-type `det_`, `mean_`, `cpm_`) | when you change cell types or the reference |
| 04 | `04_build_network.py --config <network.yaml> --out <dir>` | GRN, L-R list, modules, half-lives, optional expression table | `<dir>/genes.tsv, edges.tsv, lr_pairs.tsv, lr_downstream.tsv, report.md` | when you change the GRN, L-R list, modules or sizes |
| 04b | `04b_build_lr_downstream.py [--config <network.yaml>] [--pathway-tfs <tsv>]` | CellChat L-R list with a `pathway` column + `data/external/pathway_tfs.tsv` | `data/processed/lr_downstream_curated.tsv` (+ `.report.md`) | when you edit the pathway->TF curation, or move to a new L-R list |
| 05 | `05_engine_checks.py --config <engine.yaml> [--cells-per-type 1000] [--hours 48]` | engine config + network | `results/engine_checks/<config>/` | after changing the network or engine config |
| 06 | `06_abm_checks.py --engine <..> --abm <..> [--hours 24]` | engine + ABM configs | `results/abm_checks/<engine>__<abm>/` | after changing the ABM config |
| 07 | `07_rules_checks.py --engine <..> --abm <..> --rules <..>` | all three configs | `results/rules_checks/<e>__<a>__<r>/` | after changing rules |
| 08 | `08_run_simulation.py` (options below) | all three configs | `results/runs/<name>/` | each simulation |
| 09 | `09_analysis_plots.py --run <run dir> [--plots all\|names] [--pool-cells 30000]` | a run folder | `<run>/plots/*.png`, `<run>/movement.gif`, `<run>/analysis/{final,pooled}.h5ad` | after each run (spatial env) |
| 10 | `10_velocity_plots.py --run <run dir>` | a run folder after 09 | `<run>/plots/02*_velocity*.png`, `03b_velocity_pseudotime.png` | after 09 (velo_py env) |

**Running a script on the cluster.** Use `sbatch -J <name> scripts/run_py.sbatch scripts/<script>.py
[args]`. The job gets 48 GB, 2 CPUs and 1 h; logs go to `logs/<name>_<jobid>.log`.

**Analysis.** `scripts/run_analysis.sbatch <run dir> [plot names]` runs 09 and then 10 in their two
environments.

### `08_run_simulation.py` options

| option | default | meaning |
|---|---|---|
| `--engine` | `configs/engine.yaml` | gene engine config |
| `--abm` | `configs/abm.yaml` | tissue config |
| `--rules` | `configs/rules.yaml` | rules config; `none` runs without rules, so every cell keeps its type's parameters |
| `--hours` | 72 | simulated time after t = 0 |
| `--snapshot-min` | 180 | minutes between full snapshots (expression, positions, rules, signalling) |
| `--track-min` | 6 | minutes between position records (tracks and animation); 6 = every step |
| `--seed` | 0 | run seed. It changes the layout, motility, events **and** gene noise, so use different seeds for replicates |
| `--name` | `<engine>__<abm>__<rules>__<hours>h_seed<seed>` | output folder name under `results/runs/` |

---

## 5. Inputs you provide

You need at least a GRN, an L-R list and behaviour gene sets. A reference dataset is optional.

| input | format | required? |
|---|---|---|
| **GRN** | `edge_list`: a table (.tsv/.csv/.parquet) with columns `regulator`, `target`, optionally `sign` (+1/−1) and `weight` (>0). Or `celloracle_wide`: rows = targets in column `gene_short_name`, one column per TF, value > 0 = edge (CellOracle base GRN) | yes |
| **L-R list** | table with columns `ligand`, `receptor`. Complexes are joined by `_` (e.g. `TGFBR1_TGFBR2`); a complex's level is its weakest subunit | yes |
| **L-R → TF links** | table with columns `ligand`, `receptor`, `tf`, and optionally `sign` (±1) and `weight` (`lr.downstream`). If not given, each pair activates `lr_downstream_tfs` random selected TFs — arbitrary downstream biology. Build a curated one with script 04b | optional but recommended |
| **behaviour modules** | gene lists in `network.yaml` (`pro` / `anti` per module) | yes (≥ 1 module) |
| **mRNA half-lives** | table with columns `gene`, `halflife_h`. Missing genes get the median | provided for human |
| **reference** (only for `baseline.source: reference` or `init.source: reference`) | h5ad with cell type labels in `obs[label_key]`, positions in µm in `obsm[spatial_key]`, raw counts in a layer / X or in a separate 10x .h5 with the same barcodes | optional |

**Without a reference**, use `baseline.source: synthetic`. Each type then gets a random but
reproducible baseline profile (see `configs/engine_synthetic5.yaml`). The tissue uses
`init.source: synthetic`.

### Built network files (`network_dir`, written by 04; the simulator reads only these)

| file | columns |
|---|---|
| `genes.tsv` | `gene`; `role` (TF / ligand / receptor / target, comma-joined); `modules` (`module:pro\|anti`, `;`-joined); `halflife_h` |
| `edges.tsv` | `regulator`, `target`, `sign` (±1), `weight` |
| `lr_pairs.tsv` | `ligand`, `receptor`, `module_gene` (picked because it involves a module gene), `score` |
| `lr_downstream.tsv` | `ligand`, `receptor`, [`pathway`,] `tf`, `sign`, `weight` |
| `report.md` | build summary: sizes, module genes without regulators, reachability of each module from ligands |

You can also write these four tables by hand, from any source, and point `network_dir` at them.

---

## 6. Configuration reference

There are four YAML files. The simulator reads `engine.yaml`, `abm.yaml` and `rules.yaml`; the network
builder reads `network.yaml`. Copy a file to make a variant (e.g. `configs/engine_5types.yaml`) and
pass it with `--engine`. Output folders are named after the config files, so variants never overwrite
each other.

### 6.1 `configs/network.yaml` (network builder, step 04)

| key | default | meaning / options |
|---|---|---|
| `grn.path` | SpaceTravLR `human_base_grn.parquet` | GRN file |
| `grn.format` | `celloracle_wide` | `celloracle_wide` or `edge_list` (see §5) |
| `lr.path` | SpaceTravLR `cellchat_human.csv` | L-R list (`ligand`, `receptor`) |
| `lr.downstream` | `null` | optional L-R → TF table (signed); `null` = random links. `configs/network_curated.yaml` points at the curated one from script 04b |
| `expression` | `data/processed/reference_type_expression.tsv` | optional; used only to rank L-R pairs by expression (max `det_*`); `null` = file order |
| `halflives` | `data/processed/mrna_halflife_human.tsv` | half-life table (`gene`, `halflife_h`) |
| `modules.<name>.pro` | apoptosis, proliferation, migration, chemotaxis, necrosis lists | genes whose expression **raises** the module score |
| `modules.<name>.anti` | (as above) | genes whose expression **lowers** it; may be `[]` |
| `sizes.n_tf` | 50 | number of TFs: module genes that are TFs first, then the TFs with most edges into module genes |
| `sizes.n_lr_pairs` | 50 | number of L-R pairs |
| `sizes.pairs_per_module_gene` | 2 | L-R pairs involving a module gene are taken first, at most this many per gene |
| `sizes.max_regulators_per_target` | 6 | random cap on regulators per gene; `0` = keep all (dense priors need a cap) |
| `sizes.lr_downstream_tfs` | 3 | TFs per L-R pair when `lr.downstream` is null |
| `defaults.seed` | 0 | seed for every random choice in the build |
| `defaults.activator_fraction` | 0.7 | P(edge is activating) when the GRN has no `sign` |
| `defaults.weight_lognormal` | `[0.0, 0.5]` | \|weight\| ~ LogNormal(μ, σ) when the GRN has no `weight` |
| `defaults.lr_weight` | 1.0 | weight of every L-R → TF link |

**Module names are free.** Rules refer to modules by name. The module score is the mean log2 change of
the `pro` genes minus the mean of the `anti` genes, each measured against the cell type's baseline. A
gene can be in several modules.

### 6.2 `configs/engine.yaml` (gene engine)

| key | default | meaning / options |
|---|---|---|
| `network_dir` | `data/processed/network` | folder with the four network tables |
| `cell_types` | `tumor`, `tcell`, `macrophage` | **the** list of cell types, as `short_name: reference label`. Any number. The label is used only with a reference; use `null` otherwise. Every other config refers to these short names |
| `baseline.source` | `reference` | `reference` (per-type mean counts from data) or `synthetic` (generated) |
| `baseline.reference.table` | `data/processed/reference_type_expression.tsv` | per-gene table with a `mean_<type>` column per type (written by 03) |
| `baseline.reference.dir` | Atera breast folder | folder of the reference files |
| `baseline.reference.labels_h5ad` | `atera_breast_subset_hvg3500_biomarkers.h5ad` | cells, labels and positions |
| `baseline.reference.label_key` | `cell_type` | `obs` column holding the labels named in `cell_types` |
| `baseline.reference.spatial_key` | `spatial` | `obsm` key with x, y in µm (for `init.source: reference`) |
| `baseline.reference.counts_h5` | `cell_feature_matrix.h5` | optional full-panel 10x .h5 with the same barcodes; `null` = counts from `labels_h5ad` |
| `baseline.reference.counts_layer` | `null` | when `counts_h5` is null: the layer of `labels_h5ad` holding raw counts (`null` = `X`) |
| `baseline.synthetic.median_counts` | 0.2 | median observed counts per gene and cell |
| `baseline.synthetic.gene_sd_log` | 1.5 | spread of gene levels (natural-log SD) |
| `baseline.synthetic.type_sd_log` | 1.0 | how different the types are from each other (natural-log SD per gene and type) |
| `baseline.synthetic.seed` | 0 | seed for the synthetic profiles |
| `dt_min` | 6.0 | length of one phenotype step (minutes) |
| `capture` | 0.2 | observed counts = Poisson(capture × molecules); the detection efficiency |
| `splice_halftime_min` | 15.0 | unspliced → spliced half-time (all genes) |
| `burst_off_rate` | 0.05 | promoter ON → OFF rate (per min); mean burst = 1/0.05 = 20 min |
| `max_on_fraction` | 0.5 | each gene's highest-expressing type sits at this promoter ON fraction (sets burstiness and room to go up) |
| `hill_n` | 2 | Hill coefficient of TF regulation; half-response = the regulator's mean level across types |
| `lr_gain` | 1.0 | scales all L-R → TF effects; `0` switches cell-cell signalling off |
| `ligand_half_response` | `neighbourhood` | ligand half-response = a typical neighbourhood's summed ligand input (one cell's level × mean contacts), so signalling responds to ligand **dose**; `cell` = one cell's level, which saturates once a cell has several neighbours (only contact matters) |
| `min_level` | 0.01 | floor (molecules) for genes at ~0 in a type |
| `max_log_rate_change` | 8.0 | clip on the exponent of `k_on`, for numerical safety |
| `seed` | 0 | gene-noise seed (combined with the run `--seed`) |

### 6.3 `configs/abm.yaml` (tissue)

| key | default | meaning / options |
|---|---|---|
| `dt_mech_min` | 0.1 | mechanics / motility substep (minutes); `dt_min / dt_mech_min` substeps per phenotype step |
| `warmup_h` | 2 | before calibration, cells move (no divisions or deaths) this long, so t = 0 has a motile tissue's packing |
| `burn_in_h` | 48 | gene burn-in on the still tissue before t = 0, with virtual divisions |
| `calibration_rounds` | 4 | the burn-in is split into this many rounds; after each but the last, gene means are nudged back to baseline |
| `init.source` | `synthetic` | `synthetic` (random mix in a disc) or `reference` (crop of the reference tissue) |
| `init.synthetic.cells` | `{tumor: 1000, tcell: 1000, macrophage: 1000}` | starting count per type |
| `init.synthetic.cells_default` | 1000 | count for types not listed in `cells` |
| `init.synthetic.packing` | 0.9 | fraction of the starting disc covered by cells (sets its radius) |
| `init.synthetic.domain_scale` | 3.0 | cells are confined to a disc this many times the starting radius |
| `init.synthetic.seed` | 0 | layout seed |
| `init.reference.crop_um` | 500 | side of the square crop (µm); the domain is that box |
| `init.reference.center_um` | `null` | `[x, y]` crop centre; `null` = the crop-sized window that best contains every type |
| `contact_factor` | 1.0 | signalling contacts: distance < contact_factor × (f_i R_i + f_j R_j) |
| `defaults.*` | see below | parameters of every type unless overridden |
| `types.<type>.*` | tcell: radius 5, speed 2, persistence 5, division 2e-4; macrophage: radius 10, speed 0.5, persistence 10, division 0 | per-type overrides of any `defaults` key |

**Parameters** (`defaults:`). Each can be overridden per type, and any rule can override it per cell:

| parameter | default | unit | meaning |
|---|---|---|---|
| `radius` | 8.41 | µm | cell radius (PhysiCell default volume) |
| `repulsion` | 10.0 | µm/min | repulsion strength when discs overlap |
| `adhesion` | 0.4 | µm/min | adhesion strength within the adhesion range |
| `max_adhesion_distance` | 1.25 | × radius | adhesion range (also sets who counts as a signalling contact) |
| `speed` | 0.0 | µm/min | motility speed (0 = not motile) |
| `persistence_min` | 1.0 | min | mean time between random direction changes |
| `division_rate` | 7.2e-4 | per min | ~16 h doubling time |
| `apoptosis_rate` | 5.31667e-5 | per min | death rate |
| `apoptosis_duration_min` | 516 | min | dying cells stop moving and signalling and are removed after this |

**Adding a parameter.** Add a key under `defaults:` and rules can set it per cell straight away. It
has no effect, though, until the code reads it with `tissue.par("<name>")`. Unknown types and
misspelt parameters under `types:` are rejected when the config loads.

### 6.4 `configs/rules.yaml` (behaviour rules)

| key | default | meaning |
|---|---|---|
| `chemotaxis_radius_um` | 40 | cells within this distance are used to estimate a ligand gradient |
| `plugins` | `[]` | Python files that add new signals and behaviours (§6.5) |
| `rules` | 9 example rules | list of rule rows (below) |

**Rule row keys:**

| key | required | default | meaning |
|---|---|---|---|
| `signal` | yes | — | what the rule reads (table below) |
| `behaviour` | yes | — | what it sets (table below) |
| `half_max` | yes | — | signal value giving half the response (> 0) |
| `max` / `min` | one of them | — | target value when the response saturates: `max` for increases, `min` for decreases |
| `direction` | no | inferred from `max` / `min` | `increases` or `decreases` |
| `hill_power` | no | 2 | steepness `n` of `H = s^n / (s^n + half_max^n)`, with `s = max(signal, 0)` |
| `cell_type` | no | `"*"` | `"*"` (all), one type, or a list `[tumor, tcell]` |
| `enabled` | no | `true` | `false` keeps the row but switches it off |
| `name` | no | — | free label |

**Signals:**

| signal | value | notes |
|---|---|---|
| `<module>` | module score in **baseline SDs**: (score − type mean) / type SD, measured at t = 0 | 2 = unusually high for that type |
| `<gene>` | the gene's log2 change vs baseline, in baseline SDs | any gene in `genes.tsv` |
| `pressure` | Σ over overlapping neighbours of (1 − d/(R_i+R_j))², relative to the starting tissue's median | 1 = typical at t = 0 |
| `neighbours` | number of signalling contacts | |
| `contact_<type>` | number of contacts of that type | e.g. `contact_macrophage` |
| plugin signals | anything | e.g. `distance_from_centre`, `dying_neighbours` in `plugins/examples.py` |

Only positive signal values trigger a response, because `s = max(signal, 0)`. To respond to *low*
expression, use an `anti` gene module, or a plugin signal that returns −z.

**Behaviours:**

| behaviour | base (no rule active) | value |
|---|---|---|
| any tissue parameter (`division_rate`, `apoptosis_rate`, `speed`, `adhesion`, `repulsion`, `persistence_min`, `radius`, `apoptosis_duration_min`, `max_adhesion_distance`, or one you add) | the type's value in `abm.yaml` | the parameter's own unit |
| `gene:<GENE>` | 0 | log change of the gene's transcription ON-rate (1.0 ≈ 2.7× more bursts) |
| `chemotaxis:<LIGAND>` | 0 | 0–1 bias of the motile direction up the gradient of that ligand's expression in nearby cells |
| plugin behaviours | 0 (or plugin-defined) | e.g. `module:<MODULE>` shifts every pro gene up and every anti gene down |

**How rules combine.** All rules on the same behaviour give one value per cell:

```
value = (1 − D)·[(1 − U)·base + U·max] + D·min
U = 1 − Π(1 − H_up),   D = 1 − Π(1 − H_down)
```

`max` and `min` are the H-weighted means of the rules' targets. Decreasing rules win over increasing
ones when both saturate.

**Validation.** Every row is checked when the rules load. One error message lists every problem:
unknown key, type, signal or behaviour, a missing `half_max`, and so on.

**Example rows:**

```yaml
# contact inhibition
- {signal: pressure, behaviour: division_rate, min: 0.0, half_max: 1.5, hill_power: 4}
# T cells with high CXCR4 climb the CXCL12 gradient
- {cell_type: tcell, signal: CXCR4, behaviour: "chemotaxis:CXCL12", max: 0.8, half_max: 1.0}
# crowded cells raise a gene's transcription
- {signal: pressure, behaviour: "gene:BCL2", max: 2.0, half_max: 1.5, hill_power: 4}
# tumour cells touching macrophages loosen adhesion
- {cell_type: tumor, signal: contact_macrophage, behaviour: adhesion, min: 0.1, half_max: 2.0}
```

### 6.5 Plugins: new kinds of signal or behaviour

A plugin is a Python file listed under `plugins:` in `rules.yaml`. Paths are relative to the project
folder. No other code changes.

```python
import numpy as np
from grn_abm.rules import register_signal, register_behaviour

@register_signal("distance_from_centre")                # rules can now use signal: distance_from_centre
def _(rules, sim, arg, cache):
    return np.linalg.norm(sim.tissue.pos, axis=1)       # one value per cell

@register_behaviour("module:", check=lambda rules, sim, arg: arg in rules.modules)
def _(rules, sim, arg, value, out):                     # value = rule-combined per-cell value
    ...                                                 # write effects into out["gene_mod"] / out["params"]
```

- **Prefixes.** A name ending in `_` or `:` is a prefix: `contact_tumor` calls the `contact_`
  function with `arg="tumor"`.
- **What a behaviour can write.** A behaviour function writes into `out`:
  - `out["params"][<tissue parameter>]` — per-cell values;
  - `out["gene_mod"]` — a cells × genes log on-rate change;
  - `out["bias"]` and `out["bias_dir"]` — chemotaxis.
- **Template.** `plugins/examples.py` is a working example.

**What you can read in a plugin:**

| object | contents |
|---|---|
| `sim.tissue` | `pos` (cells × 2, µm), `type_idx`, `dying`, `cell_id`, `parent_id`, `par(name)` (per-cell parameter), `adjacency()` (sparse contacts) |
| `sim.engine` | `S`, `U` (cells × genes, molecules), `genes`, `gi` (gene → column), `ref` (type × gene baseline), `type_idx` |
| `rules` | `modules`, `types`, `signal(name, sim)`, `pressure(sim)`, `gradient(gene, sim)` |

---

## 7. Outputs

### 7.1 A run folder: `results/runs/<name>/`

| file | contents |
|---|---|
| `snapshots.h5ad` | all snapshots stacked: one row per cell per snapshot (details below) |
| `tracks.npz` | every `--track-min` minutes: `time_h`, `cell_id`, `type` (index into `types`), `x`, `y`, `dying` |
| `events.tsv` | `kind` (`division`, `death_start`, `removed`), `cell_id`, `type`, `time_h` |
| `lr_flow.npz` | `flow[time, pair, sender type, receiver type]` = Σ over contacts of ligand(sender) × receptor gate(receiver); `pairs`, `types`, `time_h` |
| `grn_activity.npz` | `tf_activity[time, type, TF]` = mean Hill activation; `W` (TF × gene signed weights; edge activity = activation × W); `tfs`, `genes`, `types`, `time_h` |
| `run_info.yaml` | configs as run, hours, number of snapshots and rows |
| `config_engine.yaml`, `config_abm.yaml`, `config_rules.yaml` | copies of the configs used |
| `plots/`, `movement.gif`, `analysis/` | written by 09 and 10 |

### 7.2 `snapshots.h5ad`

| slot | contents |
|---|---|
| `X`, `layers['spliced']` | observed spliced counts, `Poisson(capture · s)` |
| `layers['unspliced']` | observed unspliced counts, `Poisson(capture · u)` |
| `layers['S_true']`, `layers['U_true']` | true spliced and unspliced molecules |
| `obs.time_h` | snapshot time (h) |
| `obs.cell_id`, `obs.parent_id` | stable cell ID; the cell it divided from (−1 = present at t = 0). After a division the mother keeps her ID and the daughter gets a new one |
| `obs.birth_h` | time the cell appeared (0 for starting cells) |
| `obs.cell_type`, `obs.dying` | type; whether apoptosis has started |
| `obs.x`, `obs.y`, `obsm['spatial']` | position (µm) |
| `obs.module_<m>` | module score (log2 vs the type baseline, pro − anti) |
| `obs.signal_<s>` | the signals the rules read (genes and modules in baseline SDs) |
| `obs.rule_<behaviour>` | value set by the rules: division / death rate (per min), speed, gene on-rate change, chemotaxis bias; for behaviours without rules, the type's value |
| `obs.neighbours` | signalling contacts |
| `obsm['lr_received']` | cells × L-R pairs: receptor-gated signal received (order = `uns['lr_pairs']`) |
| `var` | `role`, `modules`, `halflife_h`, `beta_per_min`, `gamma_per_min`, `is_tf` |
| `uns` | `cell_types`, `lr_pairs`, `capture`, `baseline_counts` (type × gene) |

**True velocity** (molecules/min) is `var.beta_per_min · U_true − var.gamma_per_min · S_true`.

**Size.** 25 snapshots of 3,000–6,000 cells come to about 160 MB, and tracks every 6 min to about
50 MB. For long or large runs, raise `--snapshot-min` and `--track-min`.

### 7.3 Plots (`plots/`)

| file | what it shows |
|---|---|
| `01_umap` | UMAP of the last snapshot and of all snapshots pooled, by type and by time |
| `02a_velocity_stream_{final,pooled}` | velocity streams on the UMAP: scvelo estimate from counts vs **true** velocity |
| `02b_velocity_accuracy_{final,pooled}` | per-cell cosine and per-gene correlation of estimated vs true velocity; accuracy vs half-life |
| `02c_phase_portraits_{final,pooled}` | unspliced vs spliced for the most dynamic genes, with the true steady-state line and the scvelo fit |
| `03_pseudotime` | diffusion pseudotime per type vs true time |
| `03b_velocity_pseudotime` | velocity pseudotime (estimated and true) vs true time |
| `04a_tf_expression_over_time` | the 30 most dynamic TFs per type (log2 vs baseline) |
| `04b_module_scores_over_time` | module scores per type (median and IQR) |
| `04c_grn_edge_activity_over_time` | the 30 most changing GRN edges per type |
| `04d_lr_communication_over_time` | top L-R pairs over time; sender → receiver type matrices at start and end |
| `05a_emergent_behaviour` | population, division / death rates, speed, segregation, infiltration, spread, MSD |
| `05b_neighbourhood_enrichment` | type–type contact enrichment (z vs shuffled labels), start vs end |
| `06_spatial_cell_types` | tissue snapshots coloured by type |
| `07a_cell_tracks` (+ `movement.gif`) | 40 tracks per type; animation of the whole tissue |
| `08_spatial_gene_expression` | rule genes, spatially variable genes and module scores in space at start / middle / end |
| `09a_emergent_states`, `09b_states_in_space` | Leiden states within each type: share over time, behaviour, niche, markers, location |
| `10_division_death_events` | where and when cells divide and die |
| `11_pressure` | pressure maps; pressure vs division rate |
| `12_spatially_variable_genes` | Moran's I of the top 25 genes over time |
| `13_behaviour_drivers` | genes most associated with each behaviour (red = genes the rules actually read) |

**Choosing and adding plots.** `--plots 01_umap,06_spatial_types` runs only some plots. To add one,
write a function `_NN_name()` in `09_analysis_plots.py` with the `@plot` decorator.

`analysis/final.h5ad` and `analysis/pooled.h5ad` hold the embeddings, and the pooled file also holds
the Leiden `state`, for your own analysis.

---

## 8. Recipes

**Your own GRN.**
1. Set `grn.path` and `grn.format` in `network.yaml`, or in a copy of it.
2. Run `04 --config <it> --out data/processed/network_mine`.
3. Set `network_dir: data/processed/network_mine` in a copy of `engine.yaml`.
4. Run `05` on that engine config.

**Your own reference and cell types.**
1. Edit `cell_types` so each entry is `short_name: label in obs[label_key]`.
2. Set `baseline.reference.*` to your files.
3. Run `03 --config <engine.yaml>`.
4. Run `04`, which uses the new table to rank L-R pairs, then `05`.

**No reference.** Set `baseline.source: synthetic`, list as many `cell_types` as you like (values
`null`), and keep `init.source: synthetic`. `configs/engine_synthetic5.yaml` is a 5-type example.

**More or fewer cells.** Change `init.synthetic.cells` and `cells_default`.

**Real tissue layout.** Set `init.source: reference` and choose `crop_um` and `center_um`. Only the
configured types are kept, so cells of other types leave gaps.

**Switch signalling off** (a control). Set `lr_gain: 0`.

**Switch all behaviour coupling off.** Use `--rules none`.

**A knockout-style perturbation.** Add a rule that forces a gene down everywhere, for example
`{signal: neighbours, behaviour: "gene:XYZ", min: -8, half_max: 0.01}`, or write a plugin. Compare
it with an unperturbed run using the same `--seed`.

**Replicates.** Run with `--seed 0,1,2,…` and different `--name` values.

**Different time resolution.**
- `dt_min` sets the step length. Rates are per minute, so dynamics stay the same; only the
  resolution changes.
- `--snapshot-min` sets how often data is recorded.

---

## 9. Checking a configuration

Run 05, 06 and 07 after changing configs. They compare the code with the model's own equations, not
with biology. Each writes a markdown report.

**Expected results with the defaults:**

| check | expected |
|---|---|
| engine calibration | Pearson r ≥ 0.97 between simulated and baseline means |
| contact cut | ≥ 95% of downstream TFs move against their link's sign when contacts are cut (TFs whose L-R input is all one sign; a TF fed by both signs is confounded by the all-contacts cut) |
| velocity | r ≈ 0.99 against the actual change |
| two-cell equilibrium | 0.00% error |
| MSD vs theory | within about 6% |
| event counts | Poisson z within ±2 |
| per-rule response table | every rule's behaviour moves with its signal in the stated direction |

---

## 10. What the model does not do (and caveats)

- **The biology is synthetic.** Where the inputs give no signs or weights, they are random but
  seeded. L-R → TF links are random unless you give a table. Validation is internal: the simulator
  reproduces its own baselines and rules, not measured dynamics.
- **2D only.** There are no diffusing fields: signalling is by contact, and chemotaxis uses the
  expression of cells within `chemotaxis_radius_um`.
- **Cell types are fixed.** There is no differentiation or type switching. Lineage is limited to
  mother/daughter IDs.
- **t = 0 is a steady state.** Counts at t = 0 are negative binomial (bursting steady state plus
  Poisson capture), after a 48 h burn-in with virtual divisions. Real tissue growth then changes
  neighbourhoods, so expression drifts as the tissue reorganises. That drift is model behaviour, not
  an artefact.
- **Rule strength depends on the signal scale.** Gene and module signals are in baseline SDs, so a
  `half_max` of 2–3 fires only in clearly unusual cells. Small `max` values on `gene:` rules give
  small effects.
- **Velocity from counts is hard.** Counts are sparse (a median of about 0.2 per gene and cell at
  capture 0.2), so scvelo's estimate agrees only weakly with the truth (median per-cell cosine about
  0.05–0.3). The true velocity is always in the output.
- **Pseudotime.** Without differentiation, pseudotime tracks gradual drift of state, not a
  trajectory.
- **Half-lives are human.** For another species, provide your own `halflives` table.
- **Cost.** A 72 h run of 3,000 → 6,000 cells takes about 5 min. Mechanics dominate, and cost grows
  roughly linearly with cell number.

---

## 11. Troubleshooting

| message / symptom | fix |
|---|---|
| `rules.yaml problems: ...` | read the listed problems for each rule (typo in a signal, behaviour, type or key) |
| `abm.yaml problems: unknown parameter` | add the parameter to `defaults:` before using it under `types:` |
| `no mean_<type> column ...` | run `03` with the engine config that lists your types |
| `KeyError` on a gene in `gene:` / `chemotaxis:` | the gene is not in `genes.tsv`; add it to a module or rebuild the network |
| `IORegistryError ... encoding_type='null'` when reading in `velo_py` | the file was written by newer anndata with `None` in `uns`; 09 strips these, so rerun 09 |
| `Unable to open file (file locking ...)` | `export HDF5_USE_FILE_LOCKING=FALSE` |
| cells run off or pile up at the edge | increase `init.synthetic.domain_scale` or lower speeds |
| run very slow | fewer cells, a larger `--snapshot-min` / `--track-min`, or a larger `dt_mech_min` (≤ 0.2) |

## 12. Comparing with scMultiSim

`comparison/` runs scMultiSim on the same inputs and compares the two simulators. The same inputs
are:

- the GRN with its signs and relative weights;
- the L-R pairs;
- the L-R → TF links (added as receptor → TF edges);
- three cell types;
- the same cell count;
- capture 0.2.

```bash
sbatch -J sms_spatial  comparison/run_scmultisim.sbatch spatial  3000 results/comparison/scmultisim_spatial
sbatch -J sms_velocity comparison/run_scmultisim.sbatch velocity 3000 results/comparison/scmultisim_velocity
sbatch -J compare scripts/run_py.sbatch comparison/compare_scmultisim.py --ours results/runs/<run> --out results/comparison/<name>
```

The report (`comparison_report.md`, plus `plots/c1`–`c5`) compares:

- realism against the reference at matched depth;
- UMAPs and layouts;
- GRN visibility (AUROC of true edges vs non-edges);
- L-R specificity (a neighbour's ligand vs the receiver's matching vs non-matching readouts, after
  removing neighbour composition);
- dynamics, velocity, runtime and features.

## 13. Code map

| file | role |
|---|---|
| `src/grn_abm/engine.py` | gene engine (`GeneEngine`) |
| `src/grn_abm/abm.py` | tissue (`Tissue`): mechanics, motility, events, per-cell parameters |
| `src/grn_abm/rules.py` | rules (`Rules`), the signal / behaviour registries, plugin loading |
| `src/grn_abm/sim.py` | `Simulation`: start-up, `step()`, `run()`; `load_simulation()` |
| `src/grn_abm/record.py` | `Recorder`: snapshots, tracks, events, signalling and GRN activity |
| `plugins/examples.py` | plugin template |
| `docs/PLAN.md` | design decisions and their reasons |
