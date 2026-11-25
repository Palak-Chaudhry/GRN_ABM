# scMultiSim on the SAME inputs as GRN_ABM (data/processed/network): same GRN edges (signs, relative
# weights), same L-R pairs, same L-R -> TF links, 3 cell types, same cell count, same capture.
#
# Mapping (scMultiSim has no receptor-gated TF input, so the same information is expressed its way):
#   GRN               edges.tsv -> (target, regulator, effect = sign * weight * k), self-loops dropped
#                     (scMultiSim does not support them); k scales the mean |effect| to that of scMultiSim's
#                     own example GRN (GRN_params_100, mean 2.88) so the GRN is neither too weak nor too strong
#   L-R -> TF         lr_downstream.tsv -> extra GRN edges receptor subunit -> TF (effect lr_weight * k)
#   L-R pairs         lr_pairs.tsv -> cci params (receptor subunit, ligand subunit, effect k); complexes split into
#                     subunits; homophilic pairs (ligand = receptor, e.g. CDH1-CDH1) cannot be expressed and are dropped
#   cell types        Phyla3() tree, discrete populations, equal sizes (scMultiSim cannot take per-type reference
#                     expression, so its types are its own random profiles)
#   all L-R pairs act between all type pairs (cell.type.interaction off), random mixing (same.type.prob 1/3),
#   del.lr.pair = FALSE, dynamic GRN on (5 edges out + 5 in per epoch)
#   observed counts   add_expr_noise(protocol = "UMI", alpha_mean = 0.2) = capture 0.2 as in GRN_ABM
#
# Usage: Rscript comparison/run_scmultisim.R <spatial|velocity> <num_cells> <out_dir>
#   spatial  : with cell-cell interaction on a grid (scMultiSim gives no velocity in this mode)
#   velocity : no CCI, velocity on (the only way to get scMultiSim velocity)
suppressPackageStartupMessages(library(scMultiSim))
args <- commandArgs(TRUE)
MODE <- args[1]; NC <- as.integer(args[2]); OUT <- args[3]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
nd <- "/vast/djishnu/Palak/GRN_ABM/data/processed/network"
E <- read.delim(file.path(nd, "edges.tsv"), stringsAsFactors = FALSE)
LR <- read.delim(file.path(nd, "lr_pairs.tsv"), stringsAsFactors = FALSE)
D <- read.delim(file.path(nd, "lr_downstream.tsv"), stringsAsFactors = FALSE)
k <- 2.876 / mean(E$weight)
sub <- function(x) strsplit(x, "_", fixed = TRUE)[[1]]

grn <- data.frame(target = E$target, regulator = E$regulator, effect = E$sign * E$weight * k, source = "grn")
rt <- do.call(rbind, lapply(seq_len(nrow(D)), function(i)
  data.frame(target = D$tf[i], regulator = sub(D$receptor[i]), effect = D$weight[i] * k, source = "receptor_tf")))
grn <- rbind(grn, rt)
grn <- grn[grn$target != grn$regulator, ]
grn <- grn[!duplicated(grn[, 1:2]), ]
lr <- unique(do.call(rbind, lapply(seq_len(nrow(LR)), function(i)
  expand.grid(target = sub(LR$receptor[i]), regulator = sub(LR$ligand[i]), stringsAsFactors = FALSE))))
lr <- lr[lr$target != lr$regulator, ]
lr$effect <- k
write.csv(grn, file.path(OUT, "grn_used.csv"), row.names = FALSE)
write.csv(lr, file.path(OUT, "lr_used.csv"), row.names = FALSE)
cat(sprintf("GRN %d edges (%d receptor->TF), %d genes; L-R %d ligand->receptor subunit pairs; scale k = %.2f\n",
            nrow(grn), sum(grn$source == "receptor_tf"), length(unique(c(grn$target, grn$regulator))), nrow(lr), k))

opts <- list(
  GRN = grn[, 1:3], num.cells = NC, num.cifs = 50, tree = Phyla3(), discrete.cif = TRUE,
  discrete.pop.size = as.integer(rep(NC %/% 3L, 3) + c(NC %% 3L, 0L, 0L)),
  diff.cif.fraction = 0.8, do.velocity = MODE == "velocity", speed.up = TRUE,
  dynamic.GRN = list(seed = 0, num.steps = 200, cell.per.step = max(1L, NC %/% 1000L),
                     num.changing.edges = 5, weight.sd = 1, create.tf.edges = FALSE))
if (MODE == "spatial") {
  opts$cci <- list(params = lr, max.neighbors = 4, grid.size = ceiling(sqrt(NC) * 1.3), layout = "enhanced",
                   same.type.prob = 1 / 3, del.lr.pair = FALSE, start.layer = NC)
}
set.seed(0)                                  # rand.seed is ignored upstream
t0 <- Sys.time()
res <- sim_true_counts(opts)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
set.seed(1)
add_expr_noise(res, protocol = "UMI", alpha_mean = 0.2, alpha_sd = 0.02, randseed = 1)
cat(sprintf("simulated %d cells in %.0f s\n", NC, secs))

w <- function(x, f) write.csv(x, gzfile(file.path(OUT, f)), row.names = TRUE)
w(res$counts, "true_counts.csv.gz")                       # genes x cells
obs <- res$counts_obs; if (is.list(obs)) obs <- obs$counts
w(obs, "observed_counts.csv.gz")
meta <- res$cell_meta
if (!is.null(res$cci_locs)) meta <- cbind(meta, x = res$cci_locs[, 1], y = res$cci_locs[, 2])
write.csv(meta, file.path(OUT, "cell_meta.csv"), row.names = FALSE)
if (MODE == "velocity") {
  w(res$unspliced_counts, "unspliced_counts.csv.gz")
  w(res$velocity, "velocity.csv.gz")
}
# dynamic GRN: weights of every edge that changes, per cell (simulation order)
G <- res$cell_specific_grn
if (!is.null(G)) {
  W <- simplify2array(G)
  vary <- which(apply(W, c(1, 2), function(v) diff(range(v)) > 1e-8), arr.ind = TRUE)
  long <- t(sapply(seq_len(nrow(vary)), function(i) W[vary[i, 1], vary[i, 2], ]))
  rg <- colnames(G[[1]])
  if (all(grepl("^[0-9]+$", rg))) rg <- rownames(res$counts)[as.integer(rg)]   # regulator ids -> gene names
  rownames(long) <- paste0(rg[vary[, 2]], "->", rownames(G[[1]])[vary[, 1]])
  w(long, "dynamic_grn_changing_edges.csv.gz")
  cat(sprintf("dynamic GRN: %d of %d edges change across cells\n", nrow(vary), sum(W[, , 1] != 0)))
}
if (MODE == "spatial") {
  ok <- tryCatch({ w(res$cci_cell_type_param, "cci_cell_type_param.csv.gz"); TRUE }, error = function(e) FALSE)
}
writeLines(sprintf("mode: %s\ncells: %d\nseconds: %.1f\nscale_k: %.4f", MODE, NC, secs, k), file.path(OUT, "run_info.txt"))
cat("wrote", OUT, "\n")
