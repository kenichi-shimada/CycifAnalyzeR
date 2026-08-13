# Strategy / deferred work

Ideas explored during the v2 object-redesign work (branch `v2-object-redesign`)
that were intentionally **not** applied, to keep output consistent with
existing cached results and downstream code. Revisit these separately, on
their own branch/commit, once the impact on existing consumers (cached RDS
files, the EP notebook, any other analyses reading `cell_types`) has been
checked.

## 1. Replace `"outOfROI"` with real `NA` in `cell_types`

Currently, `defineCellTypes()` assigns the literal string `"outOfROI"` as a
cell type for every cell outside the ROI (`R/celltype-define.R`). At least 9
downstream sites across 7 files then defensively filter it out by name
(`celltype-graph.R`, `dimensionality_reduction.R` x2, `plot-cell-types.R`,
`plot-slide.R`, `plot-summary.R`, `spatial-area.R`, `spatial-cell-neighborhood.R`,
`spatial-helpers.R`, `spatial-tumor-stroma.R`).

Nothing structurally requires it to be a real stored factor level -- switching
to `NA` (with `within_rois(x)` as the authoritative source of truth) and
dropping the now-dead `!= "outOfROI"` filters would be a legitimate
simplification. Deferred because:
- it changes the shape/values of `cell_types()` output, which could silently
  affect any code (in this package, in the EP notebook, or elsewhere) that
  matches on the string `"outOfROI"` directly rather than via `is.na()`.
- need to check whether cached `CellTypes`/`Cycif`/`CycifStack` RDS objects on
  disk already have `"outOfROI"` baked in, which would need re-running
  `defineCellTypes()` to pick up the new convention (inconsistent old vs. new
  objects otherwise).

Also worth folding in while doing this: `cell_types` can already be `NA` for
reasons *other* than being outside the ROI -- e.g. a within-ROI cell whose
gating cascade dead-ends at a non-leaf parent partway through the hierarchy
(gets zeroed to `NA` by the leaf-only `factor(..., levels=uniq.cts)` step in
`defineCellTypes()`). Top-level unresolved cells already get the explicit
`"all_other"` leaf label, so this only affects mid-hierarchy dead-ends. Worth
deciding whether to distinguish "out of ROI" from "ambiguous/unresolved
gating" explicitly (e.g. two different reason codes) rather than collapsing
both to a single `NA`.

## 2. Expose `is_strict` as its own column in `cell_types()`

`cell_types(x, strict=TRUE)` currently masks `cell_types` to `NA` for
non-strict calls (cells matching more than one child type above `p_thres`)
and discards the underlying `is_strict` flag -- callers can't apply their own
strictness filter without re-deriving it.

Proposed: have `cell_types()` for `Cycif`/`CycifStack` always return
`is_strict` as its own column (in addition to whatever `cell_types` masking
the `strict=` argument already does), so callers can filter on it themselves
downstream. Purely additive (new column, no change to existing columns) --
low risk, but deferred alongside item 1 above to keep this round of changes
scoped to output-consistency only.
