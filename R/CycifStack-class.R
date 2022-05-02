#' A class that contains data of multiple samples from a CyCIF project.
#'
#' @include CellType-class.R
#' @export
setClass("CycifStack",
         slots = c(
           names = "character",
           suffix = "character",

           uniq_abs = "data.frame",
           used_abs = "character",

           n_samples = "numeric",
           n_cycles = "numeric",
           max_cycles = "numeric",
           n_cells = "numeric",

           cell_type = "CellTypeCycifStack",
           threshold = "data.frame",

           raw = "data.frame",
           normalized = "data.frame",
           normalize.method = "character",
           ld_coords = "data.frame",
           samples = "list"
         )
)

