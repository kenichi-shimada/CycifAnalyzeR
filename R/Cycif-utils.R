#' @include Cycif-class.R
#' @export
setMethod("names", "Cycif", function(x) x@name)

#' @export
setGeneric("abs_list", function(x) standardGeneric("abs_list"))
setMethod("abs_list", "Cycif", function(x) x@abs_list)

#' @export
setGeneric("abs_params", function(x) standardGeneric("abs_params"))
setMethod("abs_params", "Cycif", function(x) x@abs_params)

#' @export
setGeneric("abs_params<-", function(x,...) standardGeneric("abs_params<-"))
setMethod("abs_params<-", "Cycif", function(x,value){
  stopifnot(!missing(value))
  x@abs_params <- value
  validObject(x)
  return(x)
})

#' @export
setGeneric("dna", function(x) standardGeneric("dna"))
setMethod("dna", "Cycif", function(x) x@dna)

#' @export
setGeneric("nCycles", function(x) standardGeneric("nCycles"))
setMethod("nCycles", "Cycif", function(x) x@n_cycles)

#' @export
setGeneric("nCells", function(x) standardGeneric("nCells"))
setMethod("nCells", "Cycif", function(x) x@n_cells)

#' @export
setGeneric("xys", function(x) standardGeneric("xys"))
setMethod("xys", "Cycif", function(x) x@xy_coords)

#' @export
setGeneric("segProp", function(x) standardGeneric("segProp"))
setMethod("segProp", "Cycif", function(x) x@segment_property)

#' @export
setGeneric("dna_thres", function(x) standardGeneric("dna_thres"))
setMethod("dna_thres", "Cycif", function(x) x@dna_thres)

#' @export
setGeneric("cell_type", function(x) standardGeneric("cell_type"))
setMethod("cell_type", "Cycif", function(x)x@cell_type)

#' @export
setGeneric("cell_type<-", function(x,...) standardGeneric("cell_type<-"))
setMethod("cell_type<-", "Cycif", function(x,value,...){
  if(!class(value)=="CellType"){
    stop("cell_type should be CellType")
  }

  x@cell_type <- value
  ## check n_cycle, used_abs, etc

  return(x)
})



