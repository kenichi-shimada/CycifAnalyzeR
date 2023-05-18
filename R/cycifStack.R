#_ -------------------------------------------------------
# utils CycifStack-slots ----
#' Accessing slots in a CellTypeCycifStack object
#'
#' @param x A CycifStack object
#'
#' @export
setGeneric("nSamples", function(x)standardGeneric("nSamples"))
setMethod("nSamples", "CycifStack", function(x) x@n_samples)

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
setMethod("nCells", "CycifStack", function(x) x@n_cells)

#' @export
setMethod("length","CycifStack",function(x)length(x@samples))

#' @rdname nCycles
#' @export
setMethod("nCycles", "CycifStack", function(x){
  if(!is.null(x@n_cycles)){
    return(x@n_cycles)
  }else{
    return(cyApply(x,nCycles))
  }
})

#' @rdname nCycles<-
#' @export
setMethod("nCycles<-", "CycifStack", function(x,value){
  if(!is.numeric(value) | length(value) != 1){
    stop("CycifStack@nCycles should be a single numeric value now (not accepting a vector with length of x)")
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

  if(all(this.abs %in% names(x@log_normalized))){
    x@log_normalized <- x@log_normalized[this.abs]
  }

  if(all(this.abs %in% names(x@logTh_normalized))){
    x@logTh_normalized <- x@logTh_normalized[this.abs]
  }

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

#' @rdname gates
#' @export
setMethod("gates", "CycifStack", function(x)x@cell_type@gates)

#_ -------------------------------------------------------

# fun: constructor CycifStack ----
#' Instantiate and show a CycifStack object
#'
#' @param filenames quantification file (.csv)
#' @param path path to the dir containign filename
#' @param mask_type a common mask_type of each channel. If applicable, the mask_type is
#'   removed from the channel names.
#' @param object a CycifStack object
#'
#' @rdname CycifStack
#' @export
CycifStack <- function(filenames,path,mask_type="_cellRing"){
  stopifnot(all(file.exists(file.path(path,filenames))))
  cat("Trying to load ",length(filenames)," samples.\n",sep="")

  n.samples <- length(filenames)
  samples <- lapply(filenames,function(filename){
    if(grepl(mask_type,filename)){
      name <- sub(paste0("unmicst-(.+)",mask_type,"\\.csv"),"\\1",filename)
    }else{
      name <- sub(paste0("unmicst-(.+)\\.csv"),"\\1",filename)
    }
    cat("Loading ",name,"...\n",sep="")
    smpl <- Cycif(filename,path=path,mask_type=mask_type)
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
      abs_list = abs_list,
      n_samples = n.samples,
      n_cycles = n_cycles,
      max_cycles = max_cycles,
      n_cells = n_cells,
      samples = samples
  )
}

# fun: show CycifStack ----
#' @rdname CycifStack
#' @export
setMethod("show", "CycifStack", function(object) {
  df <- data.frame(
    SampleName = object@names,
    nCycles = object@n_cycles,
    nCells = object@n_cells)
  rownames(df) <- c()

  m <- matrix(object@abs_list$ab,nrow=3)
  rownames(m) <- c("R","G","B")
  colnames(m) <- paste0("Cycle",seq(object@max_cycles))
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
  abs_list <- abs_list(xs[[idx]])
  # thres <- as.data.frame(sapply(xs,function(y){
  #   if(.hasSlot(y,"threshold")){
  #     thres <- y@threshold
  #   }
  # }))
  new("CycifStack",
      n_samples = n.samples,
      names = nms,
      abs_list = abs_list,
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
    if(all(sapply(out,is,"Cycif")) && as.CycifStack){
      out <- list2CycifStack(out)
    }
  }
  return(out)
})

#_ -------------------------------------------------------

