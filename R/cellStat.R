#' Show available cells at each cycle
#'
#' @param x A Cycif object
#' @param cumulative logical. If TRUE, available cells through cycles 1 to n are computed.
#'   If FALSE, available cells at each cycle are computed.
#' @param ratio Logical. If TRUE, show ratio of available and unavailable cells
#'  at each cycle. If FALSE, actual actual number is provided.
#'
#' @return cumUsedCells() returns a matrix where rows correspond to cells
#'   and columns correspond to cycles. statUsedCells() returns a table sumarizing
#'   the numbers of available and unavailable cells at each cycle.
#'
#' @export
setGeneric("cumUsedCells", function(x,...) standardGeneric("cumUsedCells"))
setMethod("cumUsedCells", "Cycif",
          function(x){
            u <- x@used_cells
            u <- sapply(seq(ncol(u)),function(i){
              id <- rowSums(u[,seq(i),drop=F]==1)==i
              return(id)
            })
            return(u)
          }
)

#' Summarize available cells at each cycle
#' @rdname cumUsedCells
#' @order 2
#' @export
setGeneric("statUsedCells", function(x,...) standardGeneric("statUsedCells"))
setMethod("statUsedCells", "Cycif",
          function(x,cumulative=TRUE,ratio=TRUE,...){
            stopifnot(nrow(x@used_cells)>0)
            if(cumulative){
              mat <- cumUsedCells(x)
            } else{
              mat <- x@used_cells == 1 # 0 - dropped, 1 - alive, 2 - bunched
            }
            tab <- apply(mat,2,function(l){
              f <- factor(l,levels=c("TRUE","FALSE"))
              table(f)
            })
            if(ratio){
              nc.true <- tab["TRUE",]
              tab <- nc.true/nc.true[1]
            }
            return(tab)
          }
)
