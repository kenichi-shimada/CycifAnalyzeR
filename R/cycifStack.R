#_ -------------------------------------------------------
# utils CycifStack ----
#' Accessing slots in a CellTypeCycifStack object
#'
#' @param x A CycifStack object
#'
#' @export
setGeneric("nSamples", function(x)standardGeneric("nSamples"))
setMethod("nSamples", "CycifStack", function(x) x@n_samples)

#' @export
setMethod("length","CycifStack",function(x)x@n_samples)

#' @export
setMethod("names", "CycifStack", function(x) x@names)
setMethod("names<-", "CycifStack", function(x,value){
  ori.names <- x@names
  for(i in seq(ori.names)){
    x@samples[[i]]@name <- as.character(value[i])
  }
  x@names <- value
  # validObject(x)
  return(x)
})

#' @export
setMethod("abs_list", "CycifStack", function(x) x@abs_list)

#' @export
setGeneric("maxCycles", function(x)standardGeneric("maxCycles"))
setMethod("maxCycles", "CycifStack", function(x) x@max_cycles)

#' @export
setMethod("nCells", "CycifStack", function(x){
  cyApply(x,nCells,simplify=T)
})

#' @rdname nCycles
#' @export
setMethod("nCycles", "CycifStack", function(x){
  if(!is.null(x@n_cycles)){
    return(x@n_cycles)
  }else{
    return(cyApply(x,nCycles,simplify=T))
  }
})

#' @rdname nCycles<-
#' @export
setMethod("nCycles<-", "CycifStack", function(x,value){
  if(!is.numeric(value) | length(value) != 1){
    warning("CycifStack@n_cycles will not accept a vector but a scalar in the future update")
  }

  ## subsetting samples - trimming ones not reaching the cycle specified
  ncys <- cyApply(x,nCycles,simplify=T)
  is.used.smpls <- ncys >= value
  x <- x[is.used.smpls]

  x <- cyApply(x,function(this){
    smpl <- names(this)
    nCycles(this) <- value
    return(this)
  })

  x@n_cycles <- value
  x@max_cycles <- value # will be discontinued

  x@abs_list <- abs_list(x[[1]]) ## all the samples should have same # of cycles
  this.abs <- as.character(abs_list(x)$ab)

  x@names <- smpls <- names(x)
  x@n_samples <- length(x)
  x@samples <- x@samples[smpls]

  if(nrow(x@phenoData)>0){
    x@phenoData <- x@phenoData %>%
      filter(id %in% smpls) %>%
      arrange(match(smpls, id))
  }

  return(x)
})

#' @export
setMethod("within_rois", "CycifStack", function(x){
  unlist(cyApply(x,within_rois))
})
#_ -------------------------------------------------------

# fun: constructor CycifStack ----
#' Instantiate and show a CycifStack object
#'
#' @param ft_filenames quantification file (.csv)
#' @param path path to the dir containign filename
#' @param mask_type a common mask_type of each channel. If applicable, the mask_type is
#'   removed from the channel names.
#' @param object a CycifStack object
#'
#' @rdname CycifStack
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
      adata = list(),
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

# fun: check_abs ----
setGeneric("check_abs", function(x)standardGeneric("check_abs"))
setMethod("check_abs", "CycifStack", function(x){
  abs <- cyApply(x,function(x){
    ab <- abs_list(x)$ab
  })
  n.abs <- sapply(abs,length)
  if(length(unique(n.abs))>1 | n.abs[1]==0){
    warning("not all the samples have the same # of abs;\nrun 'barplot(cyApply(x,function(cy)nrow(abs_list(cy)),simplify=T),las=2)'\n")
    warning("Run nCycles() to set a fixed # of cycles and/or drop samples that don't meet the criteria.")
  }else if(length(unique(n.abs))>1 | n.abs[1]==0){
    warning("not all the samples have the same Ab names;\nrun 'apply(do.call(cbind,cyApply(x,function(cy)abs_list(cy)$ab)),1,unique)'\n")
    warning("Run set_abs() to set the identical names in the Abs.")
  }
})

# fun: set_abs ----
setMethod("set_abs<-", "CycifStack", function(x,value){
    new.abs <- value
    abs <- abs_list(x)$ab
    if(length(abs)!=length(new.abs)){
      stop("The number of 'new.abs' should be the same as the number of abs in abs_list(x)")
    }
    if(!all(names(new.abs)==abs)){
      ("The names of the 'new.abs' should be the identical to abs_list(x)$ab")
    }
    x@abs_list$ab <- new.abs
    x <- cyApply(x,function(cy){
      set_abs(cy) <- new.abs
      return(cy)
    },)
    return(x)
  }
)

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

# fun: list2CycifStack list ----

#' Converting a list to a CycifStack object
#'
#' @param x A list or a Cycif object (if Cycif, the length should be one)
#'
#' @return A CycifStack object
#' @export
#'
#' @export
setGeneric("list2CycifStack", function(x) standardGeneric("list2CycifStack"))

#'@export
setMethod("list2CycifStack", "list",function(x){
  stopifnot(all(sapply(x,is,"Cycif")) | all(sapply(x,is,"CycifStack")))
  n.samples <- sum(sapply(x,length))
  xs <- do.call(c,lapply(x,function(cs){
    if(is(cs,"CycifStack")){
      out <- cs@samples
    }else if(is(cs,"Cycif")){
      out <- cs
    }
  }))

  mask_type <- xs[[1]]@mask_type
  nms <- 	names(xs) <- sapply(xs,names)
  n_cells <- 	sapply(xs,nCells)
  n_cycles <- sapply(xs,nCycles)
  max_cycles <- max(n_cycles)

  idx <- which(n_cycles == max_cycles)[1]
  al <- abs_list(xs[[idx]])[1:2]
  for(i in seq(xs)){
    al <- al %>% dplyr::left_join(abs_list(xs[[i]]),by=c("ab","cycle"))
  }


  new("CycifStack",
      n_samples = n.samples,
      names = nms,
      abs_list = al,
      n_cycles = n_cycles,
      max_cycles = max_cycles,
      n_cells = n_cells,
      samples = xs,
      mask_type = mask_type
      # threshold = thres
  )
})

#_ -------------------------------------------------------

# fun: Cycif2CycifStack Cycif ----

#' @export
setGeneric("Cycif2CycifStack", function(x) standardGeneric("Cycif2CycifStack"))

#' @export
setMethod("Cycif2CycifStack", "Cycif",function(x){
  Cycif2CycifStack(list(x))
})

#_ -------------------------------------------------------

# fun: subsetting [, [[, [<-, [[<- ----

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
# setMethod("[<-",
#           "CycifStack",
#           function(x, i, value){
#             stopifnot(is(value,"Cycif"))
#             x@samples[[i]] <- value
#
#             n_samples <- length(x@samples)
#             nms <- 	names(x@samples) <- sapply(x@samples,names)
#             n_cells <- 	sapply(x@samples,nCells)
#             n_cycles <- sapply(x@samples,nCycles)
#             max_cycles <- max(n_cycles)
#
#             idx <- which(n_cycles == max_cycles)[1]
#             abs_list <- abs_list(x@samples[[idx]])
#
#             x@names <- nms
#             x@abs_list <- abs_list
#             x@n_samples <- x@n_samples
#             x@n_cycles <- n_cycles
#             x@max_cycles <- max_cycles
#             x@n_cells <- n_cells
#
#             # validObject(x)
#             return(x)
#           }
# )

#' @rdname CycifStack-subset
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
#' @param fun function to apply
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
    if(length(x@adata)>0){
      out@adata <- x@adata
    }
    if(length(x@interaction)>0){
      out@interaction <- x@interaction
    }
    if(length(x@neighborhoods)>0){
      out@neighborhoods <- x@neighborhoods
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

