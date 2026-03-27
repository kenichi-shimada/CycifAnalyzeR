#_ -------------------------------------------------------

# fun: constructor CycifStack ----
#' @param ft_filenames quantification file (.csv)
#' @param path path to the dir containign filename
#' @param mask_type a common mask_type of each channel. If applicable, the mask_type is
#'   removed from the channel names.
#' @param object a CycifStack object

#' @title Load Cycif or CycifStack samples into a CycifStack object
#'
#' @description
#' The `CycifStack` function loads Cycif or CycifStack samples from provided feature table filenames
#' and creates a CycifStack object containing information from these samples.
#'
#' @param ft_filenames Character vector of feature table filenames.
#' @param path The directory path where feature tables are located. Defaults to the current directory.
#' @param mask_type Character vector specifying the mask types to use for sample loading.
#'   Default is c("cellRing", "cell").
#' @param mcmicro Logical. Whether the input files are from the MCMicro platform. Defaults to FALSE.
#' @param use_scimap Logical. Whether to use scimap for additional processing. Defaults to FALSE.
#'
#' @return
#' A CycifStack object containing information from the loaded Cycif or CycifStack samples.
#'
#' @importFrom utils read.csv
#' @importFrom dplyr left_join
#' @importFrom magrittr %>%
#'
#' @examples
#' # Example usage:
#' ft_filenames <- c("unmicst-sample1_cellRing.csv", "unmicst-sample2_cell.csv")
#' stack <- CycifStack(ft_filenames)
#'
#' @export
CycifStack <- function(ft_filenames,
                       path=".",
                       mask_type=c("cellRing","cell"),
                       mcmicro=FALSE,
                       n_cycles,
                       use_scimap=FALSE){
  stopifnot(all(file.exists(file.path(path,ft_filenames))))
  cat("Trying to load ",length(ft_filenames)," samples.\n",sep="")

  n.samples <- length(ft_filenames)
  samples <- lapply(ft_filenames,function(ft_filename){
    has.mask_types <- sapply(mask_type,function(mt)grepl(paste0("_",mt,"\\."),ft_filename))
    if(any(has.mask_types)){
      mk <- mask_type[min(which(has.mask_types))]
      smpl <- sub(paste0("unmicst-(.+)_",mk,"\\.csv"),"\\1",ft_filename)
    }else{
      mk <- mask_types
      smpl <- sub(paste0("unmicst-(.+)\\.csv"),"\\1",ft_filename)
    }
    cat("Loading ",smpl,"...\n",sep="")
    smpl <- Cycif(ft_filename=ft_filename,
                  path=path,
                  mask_type=mk,
                  mcmicro=mcmicro,
                  use_scimap=use_scimap)
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

  n_cells <- sapply(samples,nCells)
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
      n_samples = n.samples,
      n_cycles = n_cycles,
      max_cycles = max_cycles,
      n_cells = n_cells,
      abs_list = abs_list,
      cell_types = list(),
      ld_coords = list(),
      samples = samples
  )
}

# fun: setValidity CycifStack ----
setValidity("CycifStack", function(object) {
  if (!is(object, "CycifStack")) {
    stop("Invalid object: 'object' is not of class 'CycifStack'.")
  }

  # if (!all(sapply(object@rois, inherits, "roi"))) {
  #   "All elements of roi_list must be of class 'roi'"
  # }
})
