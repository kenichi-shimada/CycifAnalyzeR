#_ -------------------------------------------------------

# fun: constructor CycifStack ----
#' @param ft_filenames quantification file (.csv)
#' @param path path to the dir containign filename
#' @param mask_type a common mask_type of each channel. If applicable, the mask_type is
#'   removed from the channel names.
#' @param object a CycifStack object

#' @title Load Cycif or CycifStack samples into a CycifStack object
#'
#' @description
#' The `CycifStack` function loads Cycif or CycifStack samples from provided feature table filenames
#' and creates a CycifStack object containing information from these samples.
#'
#' @param ft_filenames Character vector of feature table filenames.
#' @param path The directory path where feature tables are located. Defaults to the current directory.
#' @param mask_type Character vector specifying the mask types to use for sample loading.
#'   Default is c("cellRing", "cell").
#' @param mcmicro Logical. Whether the input files are from the MCMicro platform. Defaults to FALSE.
#' @param use_scimap Logical. Whether to use scimap for additional processing. Defaults to FALSE.
#'
#' @return
#' A CycifStack object containing information from the loaded Cycif or CycifStack samples.
#'
#' @importFrom utils read.csv
#' @importFrom dplyr left_join
#' @importFrom magrittr %>%
#'
#' @examples
#' # Example usage:
#' ft_filenames <- c("unmicst-sample1_cellRing.csv", "unmicst-sample2_cell.csv")
#' stack <- CycifStack(ft_filenames)
#'
#' @export
CycifStack <- function(ft_filenames,
                       path=".",
                       mask_type=c("cellRing","cell"),
                       mcmicro=FALSE,
                       use_scimap=FALSE){
  stopifnot(all(file.exists(file.path(path,ft_filenames))))
  cat("Trying to load ",length(ft_filenames)," samples.\n",sep="")

  n.samples <- length(ft_filenames)
  samples <- lapply(ft_filenames,function(ft_filename){
    has.mask_types <- sapply(mask_type,function(mt)grepl(paste0("_",mt,"\\."),ft_filename))
    if(any(has.mask_types)){
      mk <- mask_type[min(which(has.mask_types))]
      smpl <- sub(paste0("unmicst-(.+)_",mk,"\\.csv"),"\\1",ft_filename)
    }else{
      mk <- mask_types
      smpl <- sub(paste0("unmicst-(.+)\\.csv"),"\\1",ft_filename)
    }
    cat("Loading ",smpl,"...\n",sep="")
    smpl <- Cycif(ft_filename=ft_filename,
                  path=path,
                  mask_type=mk,
                  mcmicro=mcmicro,
                  use_scimap=use_scimap)
    return(smpl)
  })

  nms <- names(samples) <- sapply(samples,function(x){
    filename <- names(x)
    if(grepl(mask_type,filename)){
      name <- sub(paste0("unmicst-(.+)",mask_type,"\\.csv"),"\\1",filename)
    }else{
      name <- sub(paste0("unmicst-(.+)\\.csv"),"\\1",filename)
    }
  })
  names(nms) <- NULL

  n_cells <- 	sapply(samples,nCells)
  n_cycles <- sapply(samples,nCycles)
  max_cycles <- max(n_cycles)

  if(0){
    ###### "Validate if samples with max_cycles use abs_list."
  }else{
    idx <- which(n_cycles == max_cycles)[1]
  }

  abs_list <- abs_list(samples[[idx]])
  new("CycifStack",
      names = nms,
      mask_type = mask_type,
      n_samples = n.samples,
      n_cycles = n_cycles,
      max_cycles = max_cycles,
      n_cells = n_cells,
      abs_list = abs_list,
      cell_types = list(),
      ld_coords = list(),
      samples = samples
  )
}

# fun: setValidity CycifStack ----
setValidity("CycifStack", function(object) {
  if (!is(object, "CycifStack")) {
    stop("Invalid object: 'object' is not of class 'CycifStack'.")
  }

  # if (!all(sapply(object@rois, inherits, "roi"))) {
  #   "All elements of roi_list must be of class 'roi'"
  # }
})



# fun: show CycifStack ----
#' @rdname CycifStack
#' @export
setMethod("show", "CycifStack", function(object) {
  df <- data.frame(
    SampleName = object@names,
    nCycles = object@n_cycles,
    nCells = object@n_cells)
  rownames(df) <- c()

  n.ch <- max(table(abs_list(object)$cycle))
  m <- do.call(cbind,tapply(abs_list(object)$ab,abs_list(object)$cycle,
                            function(x){
                              tmp <- rep(NA,n.ch)
                              tmp[seq(x)] <- x
                              return(tmp)
                            }))

  max.cycles <- max(nCycles(object))
  rownames(m) <- paste0("Ch",seq(n.ch))
  colnames(m) <- paste0("Cycle",seq(max.cycles))
  m <- data.frame(m)

  cat(is(object)[[1]], "\n")
  print(df)
  cat("Abs:\n")
  print(m)
})



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

              ph <- x@phenoData
              if(is.data.frame(ph) & nrow(ph)>0){
                ph <- ph %>% dplyr::filter(id %in% nms)
                x@phenoData <- ph
              }

              if(0){ # to be done
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
            stopifnot(is(value,"Cycif"))
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

#' @rdname CycifStack_subset
#' @export
setMethod("[[<-",
          "CycifStack",
          function(x, i, value){
            stopifnot(is(value,"Cycif"))
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


#_ -------------------------------------------------------

# fun: cyApply CycifStack ----
#' Apply-like loop function for a CycifStack object
#' @param x A CycifStack object
#' @param fun function to apply to each Cycif object
#' @param simplify logical (FALSE by default). If TRUE, the function returns a vector or a matrix rather
#'   than a list, if applicable. It uses 'sapply' instead of 'lapply' internally.
#' @param as.CycifStack logical (TRUE by default). If TRUE, the function attempts to convert the
#'   output into a CycifStack object.
#' @param ... Additional arguments passed to sapply/lapply functions.
#'
#' @usage
#' cyApply(x,fun,simplify=FALSE,as.CycifStack=TRUE,...)
#'
#' @export
setGeneric("cyApply", function(x,...) standardGeneric("cyApply"))
setMethod("cyApply", "CycifStack", function(x,fun,simplify=FALSE,as.CycifStack=TRUE,...){
  if(simplify){
    out <- sapply(x@samples,fun,...)
  }else{
    out <- lapply(x@samples,fun,...)
  }
  if(all(sapply(out,is,"Cycif")) && as.CycifStack){
    out <- list2CycifStack(out)

    if(length(x@cell_types)>0){
      out@cell_types <- x@cell_types
    }
    if(length(x@ld_coords)>0){
      out@ld_coords <- x@ld_coords
    }
    if(nrow(x@phenoData)>0){
      pData(out) <- x@phenoData
    }
    if(length(x@calls)>0){
      out@calls <- x@calls
    }
  }

  return(out)
})

#_ -------------------------------------------------------

