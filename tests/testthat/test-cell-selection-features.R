# End-to-end coverage for the v2 object-redesign pipeline:
#   computeCN() -> computeSelection() -> subsetCells() -> runUMAP() -> clusterCells()
# Emphasis is on row-alignment: every join in this pipeline is done via an
# explicit id/mask lookup, never row position, and these tests check that
# invariant directly against independently-recomputed ground truth rather
# than just checking shapes.

set.seed(1)

abs <- c("CD3", "CD68", "PDL1")
cs_raw <- make_test_cycif_stack(abs = abs)

test_that("computeCN() produces a full-population CellNeighborhood", {
  cs <- computeCN(cs_raw, type = "knn", k = 5,
                   used.cts = c("TypeA", "TypeB", "TypeC"), mc.cores = 1)
  cn <- cell_neighborhood(cs)

  expect_s4_class(cn, "CellNeighborhood")
  expect_length(cn@is_used, sum(nCells(cs)))
  expect_equal(sum(cn@is_used), nrow(cn@neighbor_counts))
  expect_true(sum(cn@is_used) > 0)
})

test_that("computeSelection() respects used.cts and n_per_sample caps", {
  cs <- computeCN(cs_raw, type = "knn", k = 5,
                   used.cts = c("TypeA", "TypeB", "TypeC"), mc.cores = 1)

  sel <- computeSelection(cs, used.cts = c("TypeA", "TypeB"), n_per_sample = 10, seed = 42)

  expect_s4_class(sel, "CellSelection")
  expect_length(sel@is_used, sum(nCells(cs)))
  expect_true(sum(sel@is_used) <= 20) # 2 samples x 10 cap

  cts_full <- cell_types(cs)$cell_types
  expect_true(all(cts_full[sel@is_used] %in% c("TypeA", "TypeB")))
})

test_that("computeSelection() errors on a badly-named n_per_sample vector", {
  cs <- computeCN(cs_raw, type = "knn", k = 5,
                   used.cts = c("TypeA", "TypeB", "TypeC"), mc.cores = 1)

  expect_error(
    computeSelection(cs, used.cts = "TypeA", n_per_sample = c(wrong_name = 5, another = 3)),
    "must be a single number"
  )
})

test_that("subsetCells() row-aligns sample/cell_types/xy/expression against ground truth", {
  cs <- computeCN(cs_raw, type = "knn", k = 5,
                   used.cts = c("TypeA", "TypeB", "TypeC"), mc.cores = 1)
  sel <- computeSelection(cs, used.cts = c("TypeA", "TypeB"), n_per_sample = 10, seed = 42)
  cf <- subsetCells(cs, sel)

  expect_s4_class(cf, "CellFeatures")
  dt <- cf@data
  sel_idx <- which(sel@is_used)

  expect_equal(nrow(dt), length(sel_idx))
  expect_identical(dt$global_id, sel_idx)

  full_cts <- cell_types(cs)
  expect_identical(as.character(dt$sample), as.character(full_cts$sample[sel_idx]))
  expect_identical(as.character(dt$cell_types), as.character(full_cts$cell_types[sel_idx]))

  full_xy <- data.frame(data.table::rbindlist(cyApply(cs, xys)))
  expect_identical(dt$X_centroid, full_xy$X_centroid[sel_idx])
  expect_identical(dt$Y_centroid, full_xy$Y_centroid[sel_idx])

  full_expr <- exprs(cs, type = "log")
  for (ab in abs) {
    expect_identical(dt[[ab]], full_expr[[ab]][sel_idx])
  }
})

test_that("subsetCells() neighbor_counts columns match CellNeighborhood via is_used lookup, not row position", {
  cs <- computeCN(cs_raw, type = "knn", k = 5,
                   used.cts = c("TypeA", "TypeB", "TypeC"), mc.cores = 1)
  sel <- computeSelection(cs, used.cts = c("TypeA", "TypeB"), n_per_sample = 10, seed = 42)
  cf <- subsetCells(cs, sel)
  dt <- cf@data
  sel_idx <- which(sel@is_used)

  cn <- cell_neighborhood(cs)
  cn_row_of <- rep(NA_integer_, length(sel@is_used))
  cn_row_of[cn@is_used] <- seq_len(sum(cn@is_used))

  has_cn <- !is.na(cn_row_of[sel_idx])
  expect_true(any(has_cn)) # sanity: some overlap expected given used.cts overlap

  expect_counts <- cn@neighbor_counts[cn_row_of[sel_idx[has_cn]], , drop = FALSE]
  got_counts <- as.matrix(dt[has_cn, paste0("nbhd_", colnames(cn@neighbor_counts)), with = FALSE])
  dimnames(got_counts) <- dimnames(expect_counts)
  expect_equal(got_counts, expect_counts, ignore_attr = TRUE)
})

test_that("runUMAP() adds reproducible coordinate columns to CellFeatures", {
  cs <- computeCN(cs_raw, type = "knn", k = 5,
                   used.cts = c("TypeA", "TypeB", "TypeC"), mc.cores = 1)
  sel <- computeSelection(cs, used.cts = c("TypeA", "TypeB"), n_per_sample = 10, seed = 42)
  cf <- subsetCells(cs, sel)

  cf2 <- runUMAP(cf, used.abs = abs, ld_name = "umap", n_neighbors = 5, umap.seed = 1)
  expect_true(all(c("umap_x", "umap_y") %in% names(cf2@data)))
  expect_equal(nrow(cf2@data), nrow(cf@data))
  expect_false(anyNA(cf2@data$umap_x))
  expect_false(anyNA(cf2@data$umap_y))

  cf2b <- runUMAP(cf, used.abs = abs, ld_name = "umap", n_neighbors = 5, umap.seed = 1)
  expect_equal(cf2@data$umap_x, cf2b@data$umap_x)
  expect_equal(cf2@data$umap_y, cf2b@data$umap_y)
})

test_that("clusterCells() adds a cluster column to CellFeatures", {
  cs <- computeCN(cs_raw, type = "knn", k = 5,
                   used.cts = c("TypeA", "TypeB", "TypeC"), mc.cores = 1)
  sel <- computeSelection(cs, used.cts = c("TypeA", "TypeB"), n_per_sample = 10, seed = 42)
  cf <- subsetCells(cs, sel)
  cf2 <- runUMAP(cf, used.abs = abs, ld_name = "umap", n_neighbors = 5, umap.seed = 1)

  cf3 <- clusterCells(cf2, used.abs = abs, cluster_name = "cluster", k.param = 5)
  expect_true("cluster" %in% names(cf3@data))
  expect_equal(nrow(cf3@data), nrow(cf@data))
  expect_false(anyNA(cf3@data$cluster))
})
