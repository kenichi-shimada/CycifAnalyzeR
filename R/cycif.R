#_ -------------------------------------------------------
# fun: cumUsedCells Cycif ----
# fun: statUsedCells Cycif ----

#' @title Calculate Cumulative and Summary Statistics for Used Cells
#'
#' @description
#' The `cumUsedCells` function calculates the number of cells that are available (ie, not dropped by `dnaFilter`)
#' from 1st cycle to each cycle in a Cycif object.
#'
#' The `statUsedCells` function summarizes the available cells in each cycle,
#' with options to return cumulative counts and ratios.
#'
#' @param x A Cycif object.
#'
#' @param cumulative Logical. If TRUE, return the cumulative number of used cells (cumUsedCells) or
#' cumulative counts in each cycle (statUsedCells). If FALSE, return counts in each cycle (statUsedCells).
#' Default is TRUE.
#'
#' @param ratio Logical. If TRUE, return the ratio of available cells over all the cells (statUsedCells).
#' Default is TRUE.
#'
#' @param use_rois Logical. If TRUE, compute available cells within pre-specified ROIs.
#' Default is TRUE.
#'
#' @usage
#' cumUsedCells(x, use_rois = TRUE)
#' statUsedCells(x, cumulative = TRUE, ratio = TRUE, use_rois = TRUE)
#'
#' @return
#' For cumUsedCells, a matrix with the cumulative number of used cells at each cycle is returned.
#' For statUsedCells, a table summarizing the available cells is returned, with options for cumulative counts
#' and ratios.
#'
#' @export
#' @importFrom methods setMethod
setGeneric("cumUsedCells", function(x,...) standardGeneric("cumUsedCells"))

#' @rdname cumUsedCells
#' @export
setMethod("cumUsedCells", "Cycif",
          function(x,use_rois=TRUE){
            u <- x@used_cells

            u <- sapply(seq(ncol(u)),function(i){
              id <- rowSums(u[,seq(i),drop=F]==1)==i
              return(id)
            })

            if(use_rois){
              w <- x@within_rois
              if(length(w)>0){
                u <- u & w
              }
            }

            return(u)
          }
)

#' @rdname cumUsedCells
#' @export
setGeneric("statUsedCells", function(x,...) standardGeneric("statUsedCells"))

#' @rdname cumUsedCells
#' @export
setMethod("statUsedCells", "Cycif",
          function(x,cumulative=TRUE,ratio=TRUE,use_rois=TRUE){
            .Defunct(msg=paste(
              "statUsedCells() assumed dnaFilter()'s progressive per-cycle DNA-quality attrition,",
              "which is retired (dnaFilter() required interactive locator()/readline() input and",
              "is incompatible with non-interactive execution). Where @used_cells is now set",
              "uniformly across cycles (e.g. via ROI membership -- see set_used_cells_from_rois()",
              "in the hr_apm/TALAVE projects), a per-cycle retention ratio is meaningless: every",
              "cycle returns the same value. Summarize x@used_cells directly instead, e.g.",
              "table(x@used_cells) or mean(x@used_cells[,1])."
            ))
          }
)

#_ -------------------------------------------------------
