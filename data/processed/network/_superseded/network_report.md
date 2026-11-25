# Network report

## Behaviour genes
75 unique genes across 7 modules.
- not a target of any TF in the base GRN: none
- not measured in the full Atera panel (simulated anyway, half-life default 5.0 h): none

## TFs
50 selected: 7 behaviour genes that are TFs ['E2F1', 'MYC', 'SNAI1', 'SNAI2', 'TP53', 'TWIST1', 'ZEB1'], then round-robin over cell types, each time the TF expressed in that type with the highest log(1 + edges into behaviour genes) x detection.
- expressed per type (detected >= 8% and >= 20% of peak): tumor 41, tcell 38, macrophage 30
  - tumor: E2F1, MYC, TP53, FOXA1, BHLHE40, CHD2, FOS, GATA3, MBD2, MAFB, ELF1, STAT1, MEF2A, XBP1, EP300, ETS2, ESR1, MAX, REL, ZBTB7B, NR3C1, JUN, GTF2I, IRF1, KLF3, SETDB1, ELK4, MXI1, NFATC2, CEBPD, MYB, FOXD4L5, RUNX1, VEZF1, KLF13, STAT3, YY1, ZBTB7A, ID2, BACH1, SREBF1
  - tcell: MYC, TP53, ZEB1, ETS1, MAF, BHLHE40, CHD2, MBD2, MAFB, ELF1, STAT1, MEF2A, EP300, ETS2, MAX, REL, ZBTB7B, NR3C1, JUN, GTF2I, IRF1, KLF3, SETDB1, ELK4, IRF8, MXI1, NFATC2, CEBPD, FOXD4L5, RUNX1, KLF13, STAT3, YY1, PRDM1, FLI1, ZBTB7A, ID2, BACH1
  - macrophage: MAF, CHD2, FOS, MBD2, MAFB, ELF1, STAT1, MEF2A, EP300, ETS2, MAX, REL, NR3C1, JUN, GTF2I, IRF1, KLF3, SETDB1, ELK4, IRF8, NFATC2, CEBPD, FOXD4L5, RUNX1, STAT3, YY1, PRDM1, FLI1, ID2, BACH1

## Ligand-receptor pairs
1095 CellChat pairs pass expression (ligand and every receptor subunit detected in >= 2% of cells of some type).
- 50 selected (26 because they involve a behaviour gene, the rest by expression score); 34 ligand genes, 40 receptor subunit genes
- behaviour genes that are ligands/receptors but in no passing pair: ['HAVCR2', 'IL1B']
- pathways: {'THBS': 6, 'LAMININ': 5, 'CXCL': 3, 'CLEC': 2, 'CD86': 2, 'CDH': 2, 'FN1': 2, 'ICAM': 2, 'NECTIN': 2, 'TNF': 2, 'NOTCH': 2, 'SEMA4': 2, 'CCL': 1, 'CDH1': 1, 'MIF': 1, 'FASLG': 1, 'IFN-II': 1, 'IL10': 1, 'CHAD': 1, 'CD45': 1, 'PDL2': 1, 'PD-L1': 1, 'TIGIT': 1, 'CD46': 1, 'APP': 1, 'GRN': 1, 'IL16': 1, 'VISFATIN': 1, 'SEMA3': 1, 'CD226': 1}

## TF -> target edges
978 edges over 169 genes (165 genes have >= 1 regulator). Base GRN edges among selected TFs only; 158 genes capped at 6 regulators (up to 2 sampled per cell type from TFs it expresses, weighted by detection). TF->TF edges: 294 (feedback allowed).
- behaviour genes with no regulator: none

## Structure diagnostics
- TF out-degree: max 63 (EP300), median 16, TFs with none 2; top-6 TFs hold 28% of edges
- (L-R pair, receiver type) combinations: 78, with 78 distinct downstream TF sets; most-shared set used 1 times
- receivers per type: tumor 31 pairs, tcell 30 pairs, macrophage 17 pairs

## Reachability per receiver type (ligand -> receptor -> TF -> ... -> gene, only through TFs expressed in that type)
| module | genes | tumor | tcell | macrophage |
|---|---|---|---|---|
| apoptosis | 17 | 17 | 17 | 17 |
| proliferation | 17 | 17 | 17 | 17 |
| migration | 15 | 15 | 15 | 15 |
| chemotaxis | 3 | 3 | 3 | 3 |
| necrosis | 4 | 4 | 4 | 4 |
| state_tcell | 10 | 10 | 10 | 10 |
| state_macrophage | 10 | 10 | 10 | 10 |

Acceptance (each shared module reachable in every type; state modules in their own type): PASS
Note: the base GRN is dense, so reachability is easy to satisfy; how strongly a ligand moves a module is decided by the weights and is tested in the gene engine (Phase 1).

## Gene budget
169 genes: 50 TFs, 34 ligands, 40 receptor subunits, 51 other behaviour targets (genes can hold several roles).

## Defaults filled in (not from the base GRN)
- edge signs: activating with p = 0.7, else repressing (seed 0)
- edge |weight| ~ LogNormal(0.0, 0.5)
- 'expressed in a type' = detected in >= abs_min of its cells (TFs 8%, ligands/receptors 2%) and >= 20% of the gene's peak detection (removes spill-over)
- regulators per target capped at 6: up to 2 sampled per cell type from TFs it expresses
- each (L-R pair, receiver type) activates 3 TFs expressed in that type, sampled with weight = edges into its pathway's module (config pathway_modules; unmapped pathways use all modules), weight 1.0
- half-life default for genes without a consensus value: 5.05 h (19 genes)
