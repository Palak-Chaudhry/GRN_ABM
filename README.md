# GRN_ABM — gene-regulated agent-based model of a tumour–immune niche

**New users: start with [docs/USER_GUIDE.md](docs/USER_GUIDE.md).** It covers inputs, every config
parameter, how to run, and outputs. This README is the developer log.

Agents (tumour cells, T cells, macrophages) move, divide and die according to PhysiCell-style
rules. Each agent carries its own gene state. A gene engine advances that state every step, driven by:

- a signed, cell-type-specific GRN derived from SpaceTravLR,
- ligand–receptor signalling from neighbouring agents.

Gene-set scores feed back into behaviour.

Everything is Python. scMultiSim is used only as an external baseline for comparison.

## Inputs

| what | where |
|---|---|
| Base GRN (TF → target, motif prior), CellChat L–R | `Collabs/Consensus_signaling/SpaceTravLR/data/{human_base_grn.parquet,cellchat_human.csv}` |
| Trained SpaceTravLR model (signed per-cell betas) | `Collabs/Consensus_signaling/atera_breast_subset_hvg3500_biomarkers_10xcelltype_2026-09-25/` |
| Reference spatial data (Xenium WTA breast DCIS) | `Collabs/Consensus_signaling/data/atera_breast/` (subset h5ad: 169,641 cells × 3,652 genes; full panel: 18,028 genes) |
| mRNA half-lives | `data/external/saluki_halflives/` (Agarwal & Kelley 2022) → `data/processed/mrna_halflife_human.tsv` |

The three cell types are `11q13 High Grade DCIS Cells` (tumour), `T Lymphocytes` and `Macrophages`.
There is no ATAC layer.

## Layout

- `scripts/NN_*.py` — numbered, rerunnable pipeline steps
- `src/grn_abm/` — the package: gene engine, ABM, coupling
- `data/external/` — downloaded files, untouched
- `data/processed/` — derived tables
- `docs/PLAN.md` — the plan and the decisions behind it

## Steps so far

1. `scripts/01_mrna_halflives.py` builds per-gene half-lives in hours. It takes the Agarwal & Kelley
   2022 consensus ranking and maps it onto the GM12878 4sU study. Result: median 5.0 h, 5–95% range
   2–33 h.
2. `scripts/02_module_gene_availability.py` is a diagnostic: which behaviour genes are in the Atera
   panel or in the trained SpaceTravLR run.
3. `scripts/03_reference_type_expression.py` computes per-type detection and expression for all 18k
   panel genes. It submits with `sbatch --job-name=03 scripts/run_py.sbatch scripts/03_...py`.
4. `scripts/04_build_network.py` builds the simulator's network from `configs/network.yaml`.
4b. `scripts/04b_build_lr_downstream.py` turns `data/external/pathway_tfs.tsv` (curated CellChat
   pathway -> TF effectors, signed) into the receptor->TF table `lr.downstream` expects. Without it
   each L-R pair drives `lr_downstream_tfs` RANDOM TFs, so what a signal does is arbitrary.
   `configs/network_curated.yaml` is the variant wired up to use it.

## Network format (what the simulator reads)

`data/processed/network/`:

| file | columns |
|---|---|
| `genes.tsv` | gene, role (TF / ligand / receptor / target), modules, halflife_h |
| `edges.tsv` | regulator, target, sign (+1/-1), weight |
| `lr_pairs.tsv` | ligand, receptor (subunits joined by `_`), score |
| `lr_downstream.tsv` | ligand, receptor, [pathway,] tf, sign, weight (`sign` -1 = engagement represses the TF) |
| `report.md` | what was selected, plus the reachability check |

### Using a different GRN

Edit only the `grn:` block in `configs/network.yaml`:

- `format: edge_list`, pointing at a table with columns `regulator, target` and, optionally, `sign`
  and `weight`. This accepts `.tsv`, `.csv` or `.parquet`.
- `format: celloracle_wide`, for CellOracle / SpaceTravLR base GRNs.

Then run `python scripts/04_build_network.py`. Signs, weights and receptor→TF links the new GRN
doesn't provide are filled from `defaults:`.

### Changing the behaviour genes

Edit `modules:`. A new module is a new name with `pro:` and `anti:` gene lists.

### Trying alternatives without overwriting

Run `python scripts/04_build_network.py --config my.yaml --out some/dir`.

## Gene engine (`src/grn_abm/engine.py`, config `configs/engine.yaml`)

Each cell carries a promoter state (on/off), unspliced RNA and spliced RNA for every gene. One step
of `dt_min` does the following:

1. Compute the promoter on-rate:
   `on-rate = baseline × exp(Hill-TF input + receptor-gated ligand input from neighbours + ABM rule input)`.
2. Switch promoters with the exact two-state probabilities.
3. Update unspliced and spliced RNA exactly, using the per-gene half-lives.

Velocity is β·u − γ·s.

| to change | edit |
|---|---|
| cell types (any number) | `cell_types:` in `configs/engine.yaml` |
| with reference | `baseline.source: reference`, then run `scripts/03` to build `mean_<type>` columns |
| without reference | `baseline.source: synthetic` (seeded gene levels; `type_sd_log` sets how different types are) |
| signalling strength / off | `lr_gain` |

There's a 5-type example with no reference: `configs/engine_synthetic5.yaml`.

**Calibration.** Each type's baseline on-rate is solved so that the starting tissue reproduces that
type's baseline expression, including the average ligand input in the starting layout. After that,
only changes move expression: in a cell's neighbourhood, or from ABM rule inputs.

**Self-checks.** `python scripts/05_engine_checks.py --config configs/<x>.yaml` writes to
`results/engine_checks/<x>/`. The checks are internal only:

- stability
- calibration against the engine's own baseline
- feedback cycles
- downstream TFs move against each link's sign when contacts are cut (activating links fall, repressive rise)
- velocity consistency
- speed

## Tissue / ABM (`src/grn_abm/abm.py`, `src/grn_abm/sim.py`, config `configs/abm.yaml`)

Cells are discs that follow these rules:

- **Mechanics.** PhysiCell-form repulsion and adhesion, updated with Adams–Bashforth 2 steps of
  `dt_mech_min`.
- **Motility.** Persistent random walk.
- **Division.** Happens at a rate; the daughter is placed half a radius away.
- **Death.** Happens at a rate. A dying cell stops moving and stops signalling, and is removed after
  `apoptosis_duration_min`.

**Parameters.** Every key under `defaults:` in `configs/abm.yaml` is a parameter. Each type can override
it under `types:`, and code reads it per cell with `tissue.par(name)`. Rules can override any parameter
per cell for one step, so a key added to `defaults:` is immediately rule-controllable. It changes
nothing, though, until some code reads it. Unknown types and misspelt parameters are rejected when
the config loads.

**Start-up** (`sim.py`):

1. The layout relaxes.
2. The gene engine is calibrated on that layout, including each type's mRNA dilution by division.
3. Each cell starts from the bursting model's steady state. Its mRNA level is gamma-distributed, so
   observed counts at t = 0 are negative binomial.
4. The genes then run a burn-in of `burn_in_h` (48 h by default) on the still tissue, split into
   `calibration_rounds`. During it, mRNA is halved at each type's division rate (virtual divisions,
   with the cell count unchanged). After each round, gene means are nudged back to baseline.

**Each phenotype step** (`engine.dt_min`, 6 min):

1. Rules read the current state and return per-cell parameter overrides, chemotaxis and gene inputs.
2. Mechanics and motility substeps run.
3. The signalling contacts and the rule gene inputs go into the gene engine.
4. Divisions and deaths are drawn.
5. Both layers are updated, so row *i* is the same cell in each.

```python
from grn_abm import load_simulation
sim = load_simulation("configs/engine.yaml", "configs/abm.yaml", "configs/rules.yaml", root=".")
sim.run(hours=72, every_min=60, callback=lambda s: print(s.t, s.counts()))   # or call sim.step() yourself
```

| to change | edit |
|---|---|
| per-type behaviour (radius, speed, division/death rates, …) | `types:` in `configs/abm.yaml`; unlisted types use `defaults:` |
| a new tissue parameter | add it to `defaults:`; read it with `tissue.par("name")` where it should act |
| starting layout | `init.source: synthetic` (counts per type, `cells_default` for others) or `reference` (crop of the Atera tissue; auto-centred on the window that best contains every type) |

**Self-checks.** `python scripts/06_abm_checks.py --engine configs/<e>.yaml --abm configs/<a>.yaml`
writes to `results/abm_checks/<e>__<a>/`. It covers:

- two-cell equilibrium vs the analytic distance;
- mean squared displacement vs persistent-random-walk theory;
- Poisson event counts;
- packed-tissue relaxation;
- the reference crop;
- a coupled 24 h run.

## Behaviour rules (`src/grn_abm/rules.py`, config `configs/rules.yaml`)

Rules use the PhysiCell grammar: one row per rule.

```yaml
- {cell_type: tcell, signal: CXCR4, behaviour: "chemotaxis:CXCL12", direction: increases, max: 0.8, half_max: 1.0}
```

- **Response.** `H = s^n / (s^n + half_max^n)`.
- **Combining.** Rules on the same behaviour combine as
  `(1 − D)[(1 − U)·base + U·max] + D·min`.
- **Row defaults.** `cell_type: "*"` and `hill_power: 2`. `direction` follows from whether the row
  gives `max` or `min`. `cell_type` may be a list. `enabled: false` switches a row off.
- **Validation.** Every row is checked when it loads, and all of its problems are reported at once.

| signals (built in) | |
|---|---|
| `<module>`, `<gene>` | expression vs the type baseline, in baseline SDs (measured after the burn-in) |
| `pressure` | compression; 1 = the starting tissue's median |
| `neighbours`, `contact_<type>` | number of signalling contacts (of that type) |

| behaviours (built in) | |
|---|---|
| any tissue parameter | `division_rate`, `apoptosis_rate`, `speed`, `adhesion`, `persistence_min`, … (base = the type's value) |
| `gene:<GENE>` | log change of the gene's transcription on-rate (base 0) |
| `chemotaxis:<LIGAND>` | 0–1 bias up the gradient of that ligand's expression around the cell |

**Adding things:**

- **A rule:** add a row.
- **A tissue parameter:** add it to `abm.yaml` `defaults:`.
- **A new kind of signal or behaviour:** write one function in a plugin file and list the file under
  `plugins:` in `rules.yaml`.

`plugins/examples.py` is the template. It adds the signals `distance_from_centre` and
`dying_neighbours`, and the behaviour `module:<MODULE>`, which shifts every pro gene of a module up and
every anti gene down. Neither the engine, the tissue nor the simulation loop needs to change.

**Self-checks.** `python scripts/07_rules_checks.py --engine ... --abm ... --rules ...` works for any
rule set. It checks:

- the formula;
- forced module → death (Poisson);
- gradient estimation and chemotactic drift;
- pressure rules vs none;
- gene rules on vs off;
- rejection of a bad config;
- a 48 h run with a per-rule response table (does each behaviour move with its signal, in the right
  direction?);
- extensibility (a rule on adhesion; plugin signals and behaviours).

## Running and analysing a simulation

```bash
sbatch -J run scripts/run_py.sbatch scripts/08_run_simulation.py --hours 72 --snapshot-min 180 --name my_run
sbatch -J analysis scripts/run_analysis.sbatch results/runs/my_run        # 09 (plots) then 10 (velocity)
```

**Recording.** `08` records the run via `src/grn_abm/record.py` into `results/runs/<name>/`:

- `snapshots.h5ad` holds, for every snapshot:
  - observed spliced and unspliced counts;
  - true mRNA;
  - positions;
  - module scores;
  - what the rules saw and set;
  - received L-R signal.
- `tracks.npz` holds positions at every step.
- `events.tsv` lists every division and death.
- `lr_flow.npz` holds signalling from each cell type to each other type, per L-R pair.
- `grn_activity.npz` holds TF activations; multiplying them by the network weights gives GRN edge
  activity.

**Plots.** `09` writes `plots/` in the `spatial` env. Each plot is one registered function; choose a
subset with `--plots`. They cover:

- UMAP;
- pseudotime;
- GRN, module and L-R dynamics;
- emergent behaviour and neighbourhood enrichment;
- spatial cell types;
- tracks and `movement.gif`;
- spatial gene maps;
- emergent states;
- event maps;
- pressure;
- spatially variable genes;
- behaviour drivers.

**Velocity.** `10` runs in the `velo_py` env. It plots scvelo's velocity alongside the simulator's
**true** velocity on the same UMAP, with accuracy, phase portraits and velocity pseudotime.
