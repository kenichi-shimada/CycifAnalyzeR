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
              uniq_abs <- abs_list(x@samples[[idx]])

              x@names <- nms
              x@uniq_abs <- uniq_abs
              x@n_samples <- n_samples
              x@n_cycles <- n_cycles
              x@max_cycles <- max_cycles
              x@n_cells <- n_cells
              return(x)
            }
          }
)

#' @export
setMethod("[[",
          "CycifStack",
          function(x, i = NULL){
            stopifnot(length(i)==1)
            y <- x@samples[[i]]
            return(y)
          }
)

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
            uniq_abs <- abs_list(x@samples[[idx]])

            x@names <- nms
            x@uniq_abs <- uniq_abs
            x@n_samples <- x@n_samples
            x@n_cycles <- n_cycles
            x@max_cycles <- max_cycles
            x@n_cells <- n_cells

            # validObject(x)
            return(x)
          }
)

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
            uniq_abs <- abs_list(x@samples[[idx]])

            x@names <- nms
            x@uniq_abs <- uniq_abs
            x@n_samples <- x@n_samples
            x@n_cycles <- n_cycles
            x@max_cycles <- max_cycles
            x@n_cells <- n_cells

            # validObject(x)
            return(x)
          }
)
