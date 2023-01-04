#'@export
#' @export
setGeneric("ld_names", function(x)standardGeneric("ld_names"))
setMethod("ld_names", "Cycif",function(x)names(x@ld_coords))

#' @export
setMethod("ld_names", "CycifStack",function(x)names(x@ld_coords))

#' @export
setGeneric("ld_coords", function(x,...)standardGeneric("ld_coords"))
setMethod("ld_coords", "Cycif",function(x,ld_name){
  if(missing(ld_name)){
    stop("'ld_name' should be specified to retrieve ld_coords")
  }else if(!ld_name %in% ld_names(x)){
    stop("Specified 'ld_name' doesn't exist")
  }
  return(x@ld_coords[[ld_name]])
})

#'@export
setMethod("ld_coords", "CycifStack",function(x,ld_name){
  if(missing(ld_name)){
    stop("'ld_name' should be specified to retrieve ld_coords")
  }else if(!ld_name %in% ld_names(x)){
    stop("Specified 'ld_name' doesn't exist")
  }
  return(x@ld_coords[[ld_name]])
})
