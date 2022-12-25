#' A class that contains data of multiple samples from a CyCIF project.
#'
#' @include CellType-class.R
#' @export
setClass("CycifStack",
     slots = c(
       names = "character",
       suffix = "character",

       abs_list = "data.frame",

       n_samples = "numeric",
       n_cycles = "numeric",
       max_cycles = "numeric",
       n_cells = "numeric",

       cell_type = "CellTypeCycifStack",
       cell_type_full = "CellTypeCycifStack",

       ld_coords = "data.frame",
       samples = "list",

       phenoData = "data.frame",

       functions = "list"

     )
)

