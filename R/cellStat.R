## used_cells, cumulative to the cycle, output is the table
#' @export
setGeneric("cumUsedCells", function(x,...) standardGeneric("cumUsedCells"))

#' @export
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

## used_cells, cumulative to the cycle, output is the table
#' @export
setGeneric("statUsedCells", function(x,...) standardGeneric("statUsedCells"))

#' @export
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
