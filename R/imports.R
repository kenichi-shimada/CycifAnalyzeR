#_ -------------------------------------------------------

# Package-wide imports and NSE column-name bindings.
#
# R CMD check's "no visible global function definition" / "no visible binding
# for global variable" NOTEs are largely pre-existing across this package:
# base-package functions called bare (without ::), and data.table/dplyr
# column names referenced inside NSE expressions (filter(), mutate(), `[`,
# .SD, etc.) that static analysis can't resolve. Consolidated here rather
# than scattered across every file that happens to trigger one.

#' @importFrom grDevices dev.flush dev.hold dev.off pdf rainbow
#' @importFrom graphics barplot frame hist layout legend lines mtext rect title
#' @importFrom methods getClassDef is new slotNames
#' @importFrom stats aov as.dendrogram chisq.test cycle dist hclust median na.omit order.dendrogram quantile reorder sd setNames
#' @importFrom dplyr any_of arrange group_by summarize summarize_at
#' @importFrom ggplot2 element_blank element_text scale_y_reverse sym syms theme theme_bw ylab
#' @importFrom ggrepel geom_text_repel
#' @importFrom sf st_coordinates st_polygon st_sfc
NULL

# NSE column-name symbols (data.table/dplyr expressions), not real global
# objects -- not needed at runtime, silences the corresponding NOTE only.
utils::globalVariables(c(
  ".", ".I", ".N", ".SD",
  "Child", "Parent", "Text", "X_centroid", "Y_centroid", "ab", "ancestor", "descendant",
  "cell_id", "cell_type", "ct.cols", "id", "idx",
  "is.used", "is.used1", "is.used2", "mask_types", "projected_crs", "rn",
  "roi_type", "sm", "smpl", "sub.pd", "this_ab", "type", "x", "xy2", "y"
))
