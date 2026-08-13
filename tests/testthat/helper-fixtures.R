# Synthetic multi-sample CycifStack fixture, used across the CellSelection /
# CellFeatures / computeCN test files. Kept small (two samples, ~100 cells) so
# tests run fast; large enough that per-sample subsampling and neighbor-graph
# computation are non-trivial.

make_test_cycif <- function(nm, n_cells, abs, cts_levels = c("TypeA", "TypeB", "TypeC")) {
  ct_vec <- sample(cts_levels, n_cells, replace = TRUE, prob = c(0.5, 0.3, 0.2))

  xy <- data.frame(X_centroid = runif(n_cells, 0, 1000),
                    Y_centroid = runif(n_cells, 0, 1000))

  expr <- as.data.frame(matrix(rnorm(n_cells * length(abs)), nrow = n_cells,
                                dimnames = list(NULL, abs)))

  abs_list_df <- data.frame(ab = abs, cycle = seq_along(abs))

  cell_state_def <- data.frame(cell_types = cts_levels)
  for (ab in abs) cell_state_def[[ab]] <- NA_real_

  ctc <- methods::new("CellTypes",
             cell_lineage_def = data.frame(),
             cell_state_def = cell_state_def,
             markers = data.frame(),
             n_cycles = 1,
             sample_names = rep(nm, n_cells),
             cell_types = factor(ct_vec, levels = cts_levels),
             is_strict = rep(TRUE, n_cells))

  methods::new("Cycif",
      name = nm,
      file_paths = methods::new("file_paths"),
      mask_type = "cellRing",
      n_cycles = 1,
      n_cells = n_cells,
      abs_list = abs_list_df,
      dna = data.frame(DNA1 = rnorm(n_cells)),
      raw = expr,
      xy_coords = xy,
      segment_property = data.frame(),
      used_cells = matrix(1, nrow = n_cells, ncol = 1),
      rois = list(),
      within_rois = rep(TRUE, n_cells),
      dna_thres = data.frame(),
      log_normalized = expr,
      logTh_normalized = expr,
      cell_types = list(default = ctc),
      ld_coords = list(),
      cell_neighborhood = methods::new("CellNeighborhood"),
      calls = list())
}

# `n_per_smpl` controls per-sample cell counts (named, one per sample name).
make_test_cycif_stack <- function(smpl_names = c("S1", "S2"),
                                   n_per_smpl = c(S1 = 60, S2 = 45),
                                   abs = c("CD3", "CD68", "PDL1"),
                                   cts_levels = c("TypeA", "TypeB", "TypeC")) {
  lst <- lapply(smpl_names, function(nm) {
    make_test_cycif(nm, n_per_smpl[[nm]], abs, cts_levels = cts_levels)
  })
  names(lst) <- smpl_names

  cs <- list2CycifStack(lst)
  cs@cell_types <- list(default = TRUE) # ct_names(CycifStack) reads this stack-level slot
  cs
}
