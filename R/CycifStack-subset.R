#' Subsetting samples from a CycifStack object
#'
#' @param x A CycifStack object
#' @param i A numeric or character vector. If numeric, the values are integers
#'   representing indices for the Cycif objects. If character, names of the Cycif
#'   objects.
#' @param value A Cycif object to be inserted.
#'
#' @rdname CycifStack-subset
#' @export
setMethod("[",
          "CycifStack",
          function(x, i = NULL){
            if(is.null(i)){
              return(x)
            } else {
              x@samples <- x@samples[i]

              n_samples <- length(x@samples)
              nms <- 	names(x@samples) <- sapply(x@samples,names)
              n_cells <- 	sapply(x@samples,nCells)
              n_cycles <- sapply(x@samples,nCycles)
              max_cycles <- max(n_cycles)

              idx <- which(n_cycles == max_cycles)[1]
              abs_list <- abs_list(x@samples[[idx]])

              x@names <- nms
              x@abs_list <- abs_list
              x@n_samples <- n_samples
              x@n_cycles <- n_cycles
              x@max_cycles <- max_cycles
              x@n_cells <- n_cells
              return(x)
            }
          }
)

#' @rdname CycifStack-subset
#' @export
setMethod("[[",
          "CycifStack",
          function(x, i = NULL){
            stopifnot(length(i)==1)
            y <- x@samples[[i]]
            return(y)
          }
)

#' @rdname CycifStack-subset
#' @export
setMethod("[<-",
          "CycifStack",
          function(x, i, value){
            stopifnot(class(value)=="Cycif")
            x@samples[[i]] <- value

            n_samples <- length(x@samples)
            nms <- 	names(x@samples) <- sapply(x@samples,names)
            n_cells <- 	sapply(x@samples,nCells)
            n_cycles <- sapply(x@samples,nCycles)
            max_cycles <- max(n_cycles)

            idx <- which(n_cycles == max_cycles)[1]
            abs_list <- abs_list(x@samples[[idx]])

            x@names <- nms
            x@abs_list <- abs_list
            x@n_samples <- x@n_samples
            x@n_cycles <- n_cycles
            x@max_cycles <- max_cycles
            x@n_cells <- n_cells

            # validObject(x)
            return(x)
          }
)

#' @rdname CycifStack-subset
#' @export
setMethod("[[<-",
          "CycifStack",
          function(x, i, value){
            stopifnot(class(value)=="Cycif")
            x@samples[[i]] <- value

            n_samples <- length(x@samples)
            nms <- 	names(x@samples) <- sapply(x@samples,names)
            n_cells <- 	sapply(x@samples,nCells)
            n_cycles <- sapply(x@samples,nCycles)
            max_cycles <- max(n_cycles)

            idx <- which(n_cycles == max_cycles)[1]
            abs_list <- abs_list(x@samples[[idx]])

            x@names <- nms
            x@abs_list <- abs_list
            x@n_samples <- x@n_samples
            x@n_cycles <- n_cycles
            x@max_cycles <- max_cycles
            x@n_cells <- n_cells

            # validObject(x)
            return(x)
          }
)
