#_ -------------------------------------------------------

# fun: ct_names Cycif, CycifStack ----
#
#' @title Get the names of cell types in a Cycif or CycifStack object
#'
#' @description This function returns the names of cell types stored in a Cycif or CycifStack object.
#'
#' @param x A Cycif or CycifStack object.
#'
#' @return A character vector containing the names of cell types.
#'
#' @export
setGeneric("ct_names", function(x) standardGeneric("ct_names"))

#' @rdname ct_names
#' @export
setMethod("ct_names", "Cycif", function(x) names(x@cell_types))

#' @rdname ct_names
#' @export
setMethod("ct_names", "CycifStack", function(x) names(x@cell_types))
#_ -------------------------------------------------------

# fun: show CellTypes ----
#' @export
setMethod("show", "CellTypes", function(object){
  nmk <- length(object@markers$ab)
  cts <- object@cell_lineage_def$Child
  cts <- cts[!cts == "unknown"]
  nct <- length(cts)
  mty <- apply(object@cell_lineage_def[-c(1:2)],2,function(x){
    if(all(x %in% c("AND","OR","NOT","NOR",""))){
      return("lin")
    }else if(all(x %in% c("CAN",""))){
      return("str")
    }else{
      return("others")
    }
  })
  nty <- tapply(names(mty),mty,identity)
  mst <- colnames(object@cell_state_def)

  nlin <- length(nty$lin)
  nst <- length(mst)

  # cat(m)
  cat("[",is(object)[[1]],"]\n",
      # "Sample name:\t",object@name,"\n",
      "# cycles:\t",object@n_cycles,"\n\n",
      "# cell types:\t",nct,"\n",
      paste(cts,collapse=", "),"\n\n",

      "# markers in total:\t", nmk,"\n",
      "# cell lineage markers:",nlin,"\n",
      "# cell state markers:\t",nst,"\n\n",

      "lineage markers:",
      paste(nty$lin,collapse=", "),"\n",
      "state markers:",
      paste(mst,collapse=", ")
  )
})
#_ -------------------------------------------------------

# fun: cell_types Cycif, CycifStack, CellType ----

#' @title Get cell type information from a CellTypes, Cycif, or CycifStack object.
#'
#' @description This function retrieves cell type information from a CellTypes, Cycif, or CycifStack object.
#' You can specify the cell type calling method to use (if applicable) and whether
#' to include only strict cell type assignments.
#'
#' @param x A CellTypes, Cycif, or CycifStack object.
#' @param ct_name Name of the cell type calling method (default is "default").
#' @param strict Logical, whether to include only strict cell type assignments.
#' Strict assignment returns NA when more than one cell type can be assigned to a cell, e.g., CD8T cell and Tumor cell  (default is FALSE).
#'
#' @return A data frame with cell type information, including sample names and cell types.
#'
#' @seealso
#' \code{\link{Cycif}}, \code{\link{CycifStack}}, \code{\link{CellTypes}}
#'
#' @export
setGeneric("cell_types", function(x,...) standardGeneric("cell_types"))

#' @rdname cell_types
#' @export
setMethod("cell_types", "CellTypes",function(x,strict=FALSE){
  cts <- x@cell_types
  ctype <- x@cell_lineage_def

  if(strict){
    is.strict <- x@is_strict
    cts[!is.strict] <- NA
  }
  return(cts)
}) # fast

#' @rdname cell_types
#' @export
setMethod("cell_types","Cycif",
          function(x,
                   ct_name="default",
                   strict=FALSE){
  if(!ct_name %in% ct_names(x)){
    stop("ct_name doesn't exist: ",ct_name)
  }
  cts <- cell_types(x=x@cell_types[[ct_name]],strict=strict)

  smpls <- x@cell_types[[ct_name]]@sample_names
  df <- data.frame(sample=smpls,cell_types=cts)

  return(df)
})

#' @rdname cell_types
#' @export
setMethod("cell_types", "CycifStack",
          function(x,
                   ct_name="default",
                   strict=FALSE){
  if(!ct_name %in% ct_names(x)){
    stop("ct_name doesn't exist: ",ct_name)
  }
  cts <- cyApply(x,function(cy){
    cell_types(x=cy,
               ct_name=ct_name,
               strict=strict)
  })
  df <- data.frame(data.table::rbindlist(cts)) %>%
    mutate(sample=factor(sample,levels=names(x)))

  return(df)
})
