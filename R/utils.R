#_ -------------------------------------------------------

# names Cycif,CycifStack ----

#' @export
setMethod("names", "CycifStack", function(x) x@names)
setMethod("names", "Cycif", function(x) x@name)

#_ -------------
# abs_list Cycif,CycifStack ----

#' @title Get the list of antibodies used in a CycifStack object
#'
#' @description
#' The `abs_list` function returns the absolute list of markers in a CycifStack object.
#'
#' @param x A CycifStack object.
#'
#' @return
#' A data frame containing the absolute list of markers with columns "ab" and "cycle."
#'
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
#' @export
setGeneric("abs_list", function(x) standardGeneric("abs_list"))
setMethod("abs_list", "Cycif", function(x) x@abs_list)
setMethod("abs_list", "CycifStack", function(x) x@abs_list)

#_ -------------

# dna Cycif ----

# create helpfile for the dna function below
#' @title Get the DNA channels in a Cycif object
#'
#' @description
#' The `dna` function returns the DNA channels in a Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return
#' A character vector containing the names of DNA channels in the Cycif object.
#'
#' @export
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
#' @examples
#' \dontrun{
#' dna(cycif_object)
#' }
#'
setGeneric("dna", function(x) standardGeneric("dna"))
setMethod("dna", "Cycif", function(x) x@dna)

#_ -------------

# xys Cycif ----

#' @title Get the XY coordinates in a Cycif object
#'
#' @description
#' The `xys` function returns the XY coordinates in a Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return
#' A data frame containing the XY coordinates in the Cycif object.
#'
#' @export
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
setGeneric("xys", function(x) standardGeneric("xys"))
setMethod("xys", "Cycif", function(x) x@xy_coords)

#_ -------------

## roxygen2 help file for segProp()
#'
#' @title Get the segmentation property in a Cycif object
#'
#' @description
#' The `segProp` function returns the segmentation property in a Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return
#' A character vector containing the segmentation property in the Cycif object.
#'
#' @export
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
setGeneric("segProp", function(x) standardGeneric("segProp"))
setMethod("segProp", "Cycif", function(x) x@segment_property)

## roxygen2 help file for dna_thres()
#' @title Get the DNA threshold in a Cycif object
#'
#' @description
#' The `dna_thres` function returns the DNA threshold in a Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return
#' A data frame containing the DNA threshold in the Cycif object.
#'
#' @export
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
setGeneric("dna_thres", function(x) standardGeneric("dna_thres"))
setMethod("dna_thres", "Cycif", function(x) x@dna_thres)

## roxygen2 help file for used_cells()
#' @title Get the used cells in a Cycif object
#'
#' @description
#' The `used_cells` function returns the used cells in a Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return
#' A data frame containing the used cells in the Cycif object.
#'
#' @export
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
setGeneric("used_cells", function(x) standardGeneric("used_cells"))
setMethod("used_cells", "Cycif", function(x) rowSums(x@used_cells==1)==ncol(x@used_cells))

#' @title Check if samples in a CycifStack are within ROIs
#'
#' @description
#' The `within_rois` function checks if samples in a CycifStack object are within regions of interest (ROIs).
#'
#' @param x A CycifStack object.
#'
#' @return
#' A logical vector indicating whether each sample is within an ROI.
#'
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#'
#' @export
setGeneric("within_rois", function(x) standardGeneric("within_rois"))
setMethod("within_rois", "Cycif", function(x) x@within_rois)

#_ -------------

#' @title Get or Set the number of cycles in a CycifStack object
#'
#' @description
#' The `nCycles` function returns the number of cells in each sample of a CycifStack object.
#' The `nCycles<-` function sets the number of cycles in a CycifStack object. It trims samples that do not reach the specified cycle value.
#'
#' @param x A CycifStack object.
#' @param value A numeric scalar indicating the desired number of cycles.
#'
#' @return
#' The `nCycles` function returns a numeric vector containing the number of cells in each sample of the CycifStack object.
#' The `nCycles<-` function returns a modified CycifStack object with the updated number of cycles.
#'
#' @details
#' This `nCycles` function calculates the number of cells in each sample of a CycifStack object and returns the results as a numeric vector.
#' The `nCycles<-` function sets the number of cycles in a CycifStack object to the specified value. It trims samples that do not reach the specified cycle value, ensuring that all samples have the same number of cycles.
#'
#' @note
#' This function is designed to accept a numeric scalar (`value`) as the desired number of cycles. In the future, it will not accept a vector, but only a scalar value.
#'
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#' @export
setGeneric("nCycles", function(x) standardGeneric("nCycles"))

#' @rdname nCycles
#' @export
setMethod("nCycles", "Cycif", function(x) x@n_cycles)

#' @rdname nCycles
#' @export
setMethod("nCycles", "CycifStack", function(x){
  if(!is.null(x@n_cycles)){
    return(x@n_cycles)
  }else{
    return(cyApply(x,nCycles,simplify=T))
  }
})

#' @rdname nCycles
#' @export
setGeneric("nCycles<-", function(x,...,value) standardGeneric("nCycles<-"))

#' @rdname nCycles
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export
setMethod("nCycles<-", "Cycif", function(x,value){

  if(!is.numeric(value)){
    stop("value should be numeric")
  }
  value <- as.integer(value)

  current.nc <- nCycles(x)
  if(any(value > current.nc)){
    stop("New max n_cycles should be less than or equal to the current max n_cycles.")
  }

  ## n_cycles, abs_list, raw, log_normalized, logTh_normalized, dna,
  ## dna_thres
  x@n_cycles <- value

  x@abs_list <- abs_list(x) %>% filter(cycle <= value)
  this.abs <- as.character(abs_list(x)$ab)

  x@raw <- x@raw[this.abs]
  if(all(this.abs %in% names(x@log_normalized))){
    x@log_normalized <- x@log_normalized[this.abs]
  }
  if(all(this.abs %in% names(x@logTh_normalized))){
    x@logTh_normalized <- x@logTh_normalized[this.abs]
  }

  x@dna <- x@dna[paste0("DNA",seq(value))]

  if(nrow(x@used_cells)>0){
    x@used_cells <- x@used_cells[,seq(value)]
  }

  if(nrow(x@dna_thres)>0){
    x@dna_thres <- x@dna_thres[paste0("DNA",seq(value)),]
  }

  if(length(x@rois)>0){
    ncys <- sapply(x@rois,function(y)y$cycle)
    x@rois <- x@rois[ncys <= value]
  }
  return(x)
})

#' @rdname nCycles
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
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


#_ -------------------------------------------------------
# nSamples CycifStack ----
# Accessing slots in a CellTypeCycifStack object
#' @title Get the number of samples in a CycifStack object
#'
#' @description
#' The `nSamples` function returns the number of samples in a CycifStack object.
#'
#' @param x A CycifStack object.
#'
#' @return
#' The number of samples in the CycifStack object.
#'
#' @seealso
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#'
#' @export
setGeneric("nSamples", function(x)standardGeneric("nSamples"))

#' @rdname nSamples
#' @export
setMethod("nSamples", "CycifStack", function(x) x@n_samples)

#' @rdname nSamples
#' @export
setMethod("length","CycifStack",function(x)x@n_samples)

#_ -------
# names CycifStack ----

#' @export
setMethod("names<-", "CycifStack", function(x,value){
  ori.names <- x@names
  for(i in seq(ori.names)){
    x@samples[[i]]@name <- as.character(value[i])
  }
  x@names <- value
  # validObject(x)
  return(x)
})

#_ -------
# maxCycles CycifStack ----

#' @title Get the maximum number of cycles used in a CycifStack object
#'
#' @description
#' The `maxCycles` function returns the maximum number of cycles in a CycifStack object.
#'
#' @param x A CycifStack object.
#'
#' @return
#' The maximum number of cycles in the CycifStack object.
#'
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#'
#' @export
setGeneric("maxCycles", function(x)standardGeneric("maxCycles"))
setMethod("maxCycles", "CycifStack", function(x) x@max_cycles)

#_ -------
# nCells Cycif, CycifStack ----

#' @title Get the number of cells in each sample of a CycifStack object
#'
#' @description
#' The `nCells` function returns the number of cells in each sample of a CycifStack object.
#' Note that the cells are not assessed for their QC status. For the number of cells dropped
#' through DNAFiter or ROI selection, use other methods.
#'
#' @param x A CycifStack object.
#'
#' @return
#' A numeric vector containing the number of cells in each sample of the CycifStack object.
#'
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
#' @export
setGeneric("nCells", function(x) standardGeneric("nCells"))

#' @rdname nCells
#' @export
setMethod("nCells", "Cycif", function(x) x@n_cells)

#' @rdname nCells
#' @export
setMethod("nCells", "CycifStack", function(x){
  cyApply(x,nCells,simplify=T)
})

#_ -------
# set_abs CycifStack ----

#' @title Check and set antibodies for a Cycif or CycifStack object
#'
#' @description
#' The `check_abs` function checks whether all samples in a CycifStack object have the same number
#' of antibodies and whether they have identical antibody names. If not, it provides
#' suggestions on how to address the issue.
#'
#' The `set_abs` function allows you to set the
#' antibodies for a CycifStack object to a specified set of antibodies.
#'
#' @param x A CycifStack object.
#'
#' @param value A character vector containing the new set of antibodies (for `set_abs<-`).
#'
#' @return
#' For `check_abs`, a warning message is displayed if the samples do not have the same
#' number of antibodies or if the antibody names are not identical.
#'
#' For `set_abs<-`, a modified CycifStack object with the specified antibodies is returned.
#'
#' @export
setGeneric("check_abs", function(x)standardGeneric("check_abs"))

#' @rdname check_abs
#' @export
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

#' @rdname check_abs
#' @export
setGeneric("set_abs<-", function(x,...,value)standardGeneric("set_abs<-"))

#' @rdname check_abs
#' @export
setMethod("set_abs<-", "Cycif", function(x,value){
  new.abs <- value
  abs <- abs_list(x)$ab

  if(length(abs)!=length(new.abs)){
    stop("The number of 'new.abs' should be the same as the number of abs in abs_list(x)")
  }else if(!all(names(new.abs)==abs)){
    warning(names(new.abs))
    warning(paste(abs,collapse="\n"))
    stop("The names of the 'new.abs' should be the identical to abs_list(x)$ab")
  }

  x@abs_list$ab <- new.abs
  names(x@raw) <- new.abs

  if(nrow(x@log_normalized)>0){
    names(x@log_normalized) <- new.abs
  }

  if(nrow(x@logTh_normalized)>0){
    names(x@logTh_normalized) <- new.abs
  }

  return(x)
})

#' @rdname check_abs
#' @export
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
})

#_ -------------------------------------------------------
# fun: list2Cycif list ----
# fun: list2CycifStack list ----
# fun: Cycif2CycifStack Cycif ----

#' Converting a list to a Cycif object
#'
#' @param x A list object
#' @return A Cycif object
#'
#' @seealso
#' \code{\link{list2CycifStack}}, \code{\link{Cycif2CycifStack}},
#' \code{\link{Cycif-class}}, \code{\link{CycifStack-class}}
#'
#' @export
setGeneric("list2Cycif", function(x) standardGeneric("list2Cycif"))

#' @rdname list2Cycif
#' @export
setMethod("list2Cycif", "list",function(x){
  slots.in.Cycif <- slotNames(getClassDef("Cycif"))
  stopifnot(all(names(x) %in% slots.in.Cycif))
  n.samples <- sum(sapply(x,length))
  xs <- do.call(c,lapply(x,function(cs){
    if(is(cs,"Cycif")){
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
  new("Cycif",
      n_samples = n.samples,
      names = nms,
      abs_list = abs_list,
      n_cycles = n_cycles,
      max_cycles = max_cycles,
      n_cells = n_cells,
      samples = xs,
      mask_type = mask_type
  )
})

#' @title Convert a list of Cycif or CycifStack objects to a CycifStack object
#'
#' @description
#' The `list2CycifStack` and `Cycif2CycifStack` functions convert a list of Cycif or CycifStack objects into a single CycifStack object.
#'
#' @param x A list containing Cycif or CycifStack objects.
#'
#' @return
#' A CycifStack object containing information from the input list of Cycif or CycifStack objects.
#'
#' @details
#' The `list2CycifStack` function takes a list of Cycif or CycifStack objects and combines them into a single CycifStack.
#' It extracts information such as the number of samples, names, antibodies, cycles, and cell counts from each input object,
#' and creates a new CycifStack object in which each element is a Cycif object.
#' The `Cycif2CycifStack` function handles a special case where a single Cycif object is converted to a CycifStack object.
#'
#' @seealso
#' \code{\link{list2Cycif}}, \code{\link{list2CycifStack}}, \code{\link{Cycif2CycifStack}}, \code{\link{Cycif}}, \code{\link{CycifStack}}, \code{\link{CellTypes}}
setGeneric("list2CycifStack", function(x) standardGeneric("list2CycifStack"))

#' @rdname list2CycifStack
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

#' @rdname list2CycifStack
#' @export
setGeneric("Cycif2CycifStack", function(x) standardGeneric("Cycif2CycifStack"))

#' @rdname list2CycifStack
#' @export
setMethod("Cycif2CycifStack", "Cycif",function(x){
  list2CycifStack(list(x))
})
