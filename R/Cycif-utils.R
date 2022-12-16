#' @include Cycif-class.R
#' @rdname Cycif-slots
#' @export
setMethod("names", "Cycif", function(x) x@name)

#' @rdname Cycif-slots
#' @export
setGeneric("abs_list", function(x) standardGeneric("abs_list"))
setMethod("abs_list", "Cycif", function(x) x@abs_list)

#' @rdname Cycif-slots
#' @export
setGeneric("abs_params", function(x) standardGeneric("abs_params"))
setMethod("abs_params", "Cycif", function(x) x@abs_params)

#' @rdname Cycif-slots
#' @export
setGeneric("abs_params<-", function(x,...) standardGeneric("abs_params<-"))
setMethod("abs_params<-", "Cycif", function(x,value){
  stopifnot(!missing(value))
  x@abs_params <- value
  validObject(x)
  return(x)
})

#' @rdname Cycif-slots
#' @export
setGeneric("dna", function(x) standardGeneric("dna"))
setMethod("dna", "Cycif", function(x) x@dna)

#' @rdname Cycif-slots
#' @export
setGeneric("nCells", function(x) standardGeneric("nCells"))
setMethod("nCells", "Cycif", function(x) x@n_cells)

#' @rdname Cycif-slots
#' @export
setGeneric("xys", function(x) standardGeneric("xys"))
setMethod("xys", "Cycif", function(x) x@xy_coords)

#' @rdname Cycif-slots
#' @export
setGeneric("segProp", function(x) standardGeneric("segProp"))
setMethod("segProp", "Cycif", function(x) x@segment_property)

#' @rdname Cycif-slots
#' @export
setGeneric("dna_thres", function(x) standardGeneric("dna_thres"))
setMethod("dna_thres", "Cycif", function(x) x@dna_thres)

#' @rdname Cycif-slots
#' @export
setGeneric("used_cells", function(x) standardGeneric("used_cells"))
setMethod("used_cells", "Cycif", function(x) x@used_cells)

#' @rdname Cycif-slots
#' @export
setGeneric("within_rois", function(x) standardGeneric("within_rois"))
setMethod("within_rois", "Cycif", function(x) x@within_rois)

#' @rdname Cycif-slots
#' @export
setGeneric("within_rois<-", function(x,...) standardGeneric("within_rois<-"))
setMethod("within_rois<-", "Cycif", function(x,value){
  stopifnot(!missing(value))
  x@within_rois <- value
  validObject(x)
  return(x)
})

#' #' @rdname Cycif-slots
#' #' @export
#' setGeneric("cell_type", function(x) standardGeneric("cell_type"))
#' setMethod("cell_type", "Cycif", function(x)x@cell_type)
#' setMethod("cell_type", "CycifStack", function(x)x@cell_type)
#'
#' #' @rdname Cycif-slots
#' #' @export
#' setGeneric("cell_type<-", function(x,...) standardGeneric("cell_type<-"))
#' setMethod("cell_type<-", "Cycif", function(x,value,...){
#'   if(!class(value)=="CellType"){
#'     stop("cell_type should be CellType")
#'   }
#'
#'   x@cell_type <- value
#'
#'   return(x)
#' })



