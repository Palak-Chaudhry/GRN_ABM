# Plan

Status as of 2026-11-24.

## Architecture

Everything is in Python, split into two layers that alternate each phenotype step (Δt ≈ 6 min).

### 1. ABM (PhysiCell-style)

- Mechanics: adhesion and repulsion.
- Motility: persistent random walk, with optional chemotaxis.
- Division, death and removal.
- Hill-response rules that turn gene-set scores into behaviour rates.
- Rules may also change transcription of the migration L/R genes.

### 2. Gene engine

State per cell is unspliced and spliced mRNA for every gene, plus an on/off promoter state per
gene. Each step:

1. The promoter switches on or off at random (telegraph bursting). The ON rate follows
   SERGIO-style Hill regulation, driven by:
   - TFs, using signed weights per cell type;
   - ligand from neighbours, gated by the cell's own receptor.
2. Transcription happens only while the promoter is ON.
3. Unspliced RNA is spliced, then degraded with per-gene half-lives.
4. Velocity = β·u − γ·s.

All genes are stepped together, so GRN feedback loops are allowed.

### The loop, each step

1. The ABM moves cells and builds the neighbour graph.
2. Gene-set scores pass through the rules to give event rates and L/R transcription changes.
3. The gene engine advances by Δt. Daughter cells split mRNA binomially; dead cells stop
   signalling.

### Output

- AnnData with spliced and unspliced layers, spatial coordinates, and technical noise.
- Ground truth: the GRN per type, L–R contacts over time, events and lineage.
- Optional: rank-transfer onto the reference's marginal distributions for realism.

## Phases

| # | Phase | Done when |
|---|---|---|
| 0 | Network: about 50 TFs, 50 ligands, 50 receptors and 50 behaviour targets, with signed per-type weights from SpaceTravLR | Every readout gene has an incoming edge, and every module is reachable from at least one ligand |
| 1 | Gene engine (bursting, Hill regulation, unspliced/spliced, L–R gating) | It reaches a steady state, feedback loops stay stable, and 3k cells × 720 steps run in minutes |
| 2 | Baseline: scMultiSim, SERGIO and the engine compared against the reference | Mean–variance relationship, zero fraction, library size and gene–gene correlation are compared |
| 3 | ABM physics | Two-cell equilibrium, MSD matching theory, and stability at 3k+ cells |
| 4 | Coupling and rules | Demo: a ligand from A raises apoptosis in neighbouring B cells, which then die more often |
| 5 | Output and ground truth | Ground truth reloads, and scVelo agrees with the true velocity |
| 6 | Realism transfer (optional) | Per-type statistics match the reference, and edges are still recoverable |

## Decisions, and why

- **Python, not an R fork of scMultiSim.** We measured it at 2.5× faster. The existing
  scMultiSim spatial loop would take about 9 h for 3 days; see
  `Collabs/scMultiSim_abm/engine_sketch/`.
- **Gene state persists between steps.** scMultiSim's memoryless steps would erase ABM changes
  to L/R expression after a single step.
- **SERGIO's equations, not its code.** SERGIO's code hangs on GRN cycles, has no individual
  cells, and cannot be stepped one step at a time.
- **Bursting is added** to match the count-level realism of scMultiSim and SymSim.
- **No ATAC.**
- **mRNA half-lives** come from the Agarwal & Kelley 2022 consensus, scaled to hours with GM12878.
  This underestimates how fast the very shortest-lived transcripts decay.

## Open decisions

1. **Behaviour genes missing from the trained SpaceTravLR model.** Most canonical apoptosis genes
   and the EMT TFs are in the full 18k panel, but not in the 3.5k subset that SpaceTravLR was
   trained on.
2. **The tumour population.** Use `11q13 High Grade DCIS Cells` alone, or merge in its Mitotic and
   G1/S states.
