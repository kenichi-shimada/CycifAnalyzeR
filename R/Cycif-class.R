#' A class that contains data of one sample in a CyCIF project.
#'
#' @include CellType-class.R
#'
#' @export
setClass("Cycif",
   slots = c(
     ## using qts
     name = "character",
     suffix = "character",

     abs_list = "data.frame",

     raw = "data.frame",
     log_normalized = "data.frame",
     logTh_normalized = "data.frame",
     dna = "data.frame",

     n_cycles = "numeric",
     n_cells = "numeric",

     xy_coords = "data.frame",
     segment_property = "data.frame",

     cell_type = "CellTypeCycif",
     cell_type_full = "CellTypeCycif",

     used_cells = "matrix", # matrix
     within_rois = "logical", # matrix
     dna_thres = "data.frame",
     positive_rois= "list",

     ld_coords = "data.frame",
     clusters = "numeric", # numeric

     tif_files = "data.frame", # full path
     dim = "data.frame", # size of the tif_file
     rect = "list",

     segmentation_masks = "character",

     commands = "list",

     functions = "list"
   )
)

