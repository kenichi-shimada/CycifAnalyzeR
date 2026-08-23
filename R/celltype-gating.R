#_ -------------------------------------------------------

# getGates Cycif, CycifStack

#' @title Get or set gates information in a Cycif or CycifStack object
#'
#' @description This function allows you to get or set gates information in a Cycif or CycifStack object.
#'
#' @param x A Cycif or CycifStack object.
#' @param gates.df (For `setGates`) a data.frame of gate definitions.
#' @param run_normalize (For `setGates`) logical, whether to re-run \code{normalize()} after setting gates. Default TRUE.
#' @param p_thres (For `setGates`) numeric threshold passed to \code{normalize()}. Default 0.5.
#' @param trim (For `setGates`) numeric trim value passed to \code{normalize()}. Default 1e-3.
#' @param ... Additional arguments specific to the get or set operation.
#'
#' @return If used as `getGates(x)`, it returns the gates information.
#' If used as `setGates(x, ...)`, it sets the gates information and returns the modified object.
#'
#' @examples
#' \dontrun{
#' # Illustrative only -- requires a real Cycif/CycifStack object and gates_data.
#' # Get gates information from a Cycif object
#' getGates(ct_obj)
#'
#' # Set gates information in a Cycif object
#' setGates(ct_obj, gates_data)
#'
#' # Get gates information from a CycifStack object
#' getGates(ct_stack_obj)
#'
#' # Set gates information in a CycifStack object
#' setGates(ct_stack_obj, gates_data)
#' }
#'
#' @rdname getGates
#' @importFrom dplyr left_join select
#' @importFrom tibble rownames_to_column
#'
#' @export
setGeneric("getGates", function(x)standardGeneric("getGates"))

#' @rdname getGates
#' @export
setMethod("getGates", "Cycif", function(x){
  gate.names <- paste0("gates_",names(x))
  if(!all(gate.names %in% names(abs_list(x)))){
    stop("Run 'setGates()' to insert gates in a CycifStack object")
  }else{
    out <- abs_list(x) %>% dplyr::select(any_of(c("ab",gate.names)))
  }
  return(out)
})

#' @rdname getGates
#' @export
setMethod("getGates", "CycifStack", function(x){
  gate.names <- paste0("gates_",names(x))
  if(!all(gate.names %in% names(abs_list(x)))){
    stop("Run 'setGates()' to insert gates in a CycifStack object")
  }else{
    out <- abs_list(x) %>% dplyr::select(any_of(c("ab",gate.names)))
  }
  return(out)
})

# fun: setGates Cycif, CycifStack ----

#' @rdname getGates
#' @export
setGeneric("setGates", function(x,...)standardGeneric("setGates"))

#' @rdname getGates
#' @export
setMethod("setGates", "Cycif", function(x,gates.df,run_normalize=TRUE,p_thres=0.5,trim=1e-3){
  if(!is(gates.df,"data.frame")){
    "input should be in the data.frame"
  }

  if(all(rownames(gates.df) == as.character(seq(nrow(gates.df))))){
    if(!"ab" %in% names(gates.df)){
      stop("input obj should be a data.frame which contains two columns named 'ab', and 'gates_{smpl}'")
    }
  }else if(any(rownames(gates.df) %in% abs_list(x)$ab)){
    gates.df <- gates.df %>% tibble::rownames_to_column("ab")
  }

  this.col <- paste0("gates_",names(x))
  if(length(this.col)!=1 || !this.col %in% names(gates.df)){
    stop("the name of the column should be 'gates_*', where * are sample names")
  }else{
    gates.df <- gates.df %>% dplyr::select(any_of(c("ab",this.col)))
  }

  ## If gates_{smpl} already exists on x (e.g. setGates() run again on the same object,
  ## a rerun within a session), overwrite it rather than erroring -- drop the stale column
  ## before joining the fresh one in.
  existing.abs_list <- abs_list(x) %>% dplyr::select(-any_of(this.col))
  x@abs_list <- existing.abs_list %>% dplyr::left_join(gates.df,by="ab")

  if(run_normalize){
    x <- normalize(x,method="logTh",p_thres=p_thres,trim=trim)
  }

  return(x)
})

#' @rdname getGates
#' @export
setMethod("setGates", "CycifStack", function(x,gates.df,run_normalize=TRUE,p_thres=0.5,trim=1e-3){
  if(!is(gates.df,"data.frame")){
    "input should be in the data.frame"
  }

  if(all(rownames(gates.df) == as.character(seq(nrow(gates.df))))){
    if(!"ab" %in% names(gates.df)){
      stop("either rownames should be abs or one of the column's names should be 'ab'")
    }
  }else if(any(rownames(gates.df) %in% abs_list(x)$ab)){
    gates.df <- gates.df %>% tibble::rownames_to_column("ab")
  }

  this.col <- paste0("gates_",names(x))
  gates.df <- gates.df %>% dplyr::select(any_of(c("ab",this.col)))

  if(any(which(!this.col %in% names(gates.df)))){
    non.existing <- sub("gates_","",this.col[!this.col %in% names(gates.df)])
    stop("There are non-existing samples in the gates: ",non.existing)
  }

  x <- cyApply(x,function(cy){
    cat(names(cy),"..\n")
    nm <- names(cy)
    cy <- setGates(cy,gates.df,run_normalize=run_normalize,p_thres=p_thres,trim=trim)
    return(cy)
  },as.CycifStack=TRUE)

  return(x)
})
