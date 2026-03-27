#_ -------------------------------------------------------

# fun: subsetting [, [[, [<-, [[<- ----

#' Subsetting samples from a CycifStack object
#'
#' @param x A CycifStack object
#' @param i A numeric or character vector. If numeric, the values are integers
#'   representing indices for the Cycif objects. If character, names of the Cycif
#'   objects.
#' @param value A Cycif object to be inserted or replaced.
#'
#' @rdname CycifStack_subset
#' @export
setMethod("[",
          "CycifStack",
          function(x, i = NULL){
            if(is.null(i)){
              return(x)
            } else {
              x@samples <- x@samples[i]

              n_samples <- length(x@samples)
              nms <- names(x@samples) <- sapply(x@samples,names)
              n_cells <- sapply(x@samples,nCells)
              n_cycles <- x@n_cycles
              max_cycles <- max(n_cycles)

              idx <- 1
              abs_list <- abs_list(x@samples[[idx]])

              x@names <- nms
              x@abs_list <- abs_list
              x@n_samples <- n_samples
              x@n_cycles <- n_cycles
              x@max_cycles <- max_cycles
              x@n_cells <- n_cells

              ph <- x@phenoData
              if(is.data.frame(ph) & nrow(ph)>0){
                ph <- ph %>% dplyr::filter(id %in% nms)
                x@phenoData <- ph
              }

              if(0){
                ct_names <- ct_names(x)
                for(ct_name in ct_names){
                  ct <- x@cell_types[[ct_name]]
                  ct@name <- nms
                }
              }

              return(x)
            }
          }
)

#' @rdname CycifStack_subset
#' @export
setMethod("[[",
          "CycifStack",
          function(x, i = NULL){
            stopifnot(length(i)==1)
            y <- x@samples[[i]]
            return(y)
          }
)

#' @rdname CycifStack_subset
#' @export
setMethod("[<-",
          "CycifStack",
          function(x, i, value){
            stop("Replacing multiple samples in a CycifStack object is not supported.",
                 "Replace individual Cycif objects one by one using `[[<-`")
          }
)

#' @rdname CycifStack_subset
#' @export
setMethod("[[<-",
          "CycifStack",
          function(x, i, value){
            stopifnot(is(value,"Cycif"))
            n_cycles <- x@n_cycles
            if(n_cycles != value@n_cycles){
              stop("Number of cycles in the new sample does not match the existing samples.")
            }
            x@samples[[i]] <- value

            n_samples <- length(x@samples)
            nms <- names(x@samples) <- sapply(x@samples,names)
            n_cells <- sapply(x@samples,nCells)

            idx <- 1
            abs_list <- abs_list(x@samples[[idx]])

            x@names <- nms
            x@abs_list <- abs_list
            x@n_cells <- n_cells

            return(x)
          }
)
