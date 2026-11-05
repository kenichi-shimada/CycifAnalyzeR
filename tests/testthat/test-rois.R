## importROIs()/setWithinROIs(): raw cell coordinates are compared with the raw ROI polygons (no Y flip), and
## boundary handling of negative ROIs.
make_roi_df <- function(text, pts) {
  data.frame(type = "Polygon", Text = text, all_points = pts, stringsAsFactors = FALSE)
}

roi_cycif <- function(xy) {
  cy <- make_test_cycif("S1", nrow(xy), abs = c("A", "B"))
  cy@xy_coords <- xy
  cy
}

test_that("a cell on the edge of a NEGATIVE ROI is excluded; on the edge of a POSITIVE ROI it is not included", {
  ## positive ROI = [0,1000]^2 box; negative ROI = box with x in [400,600], y in [400,600]
  xy <- data.frame(X_centroid = c(500, 400, 100, 1000),   # inside neg | on neg edge | outside neg | on pos edge
                   Y_centroid = c(500, 500, 100, 500))
  rois <- rbind(make_roi_df("pos1", "0,0 1000,0 1000,1000 0,1000"),
                make_roi_df("neg1", "400,400 600,400 600,600 400,600"))
  cy <- importROIs(roi_cycif(xy), rois)
  expect_equal(cy@within_rois, c(FALSE,   # strictly inside the negative ROI
                                 FALSE,   # on the negative ROI's edge -> excluded
                                 TRUE,    # in the positive ROI, outside the negative one
                                 FALSE))  # on the positive ROI's edge -> not strictly inside
})

test_that("cells and ROI polygons are compared in the same (unflipped) frame", {
  ## cells at Y = 100, 900, 500; polygon covers Y in [-10, 200]: only the Y = 100 cell is inside
  xy <- data.frame(X_centroid = c(500, 500, 100), Y_centroid = c(100, 900, 500))
  cy <- importROIs(roi_cycif(xy), make_roi_df("pos1", "0,-10 1000,-10 1000,200 0,200"))
  expect_equal(cy@within_rois, c(TRUE, FALSE, FALSE))
})
