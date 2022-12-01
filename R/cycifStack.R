#' Instantiate and show a CycifStack object
#'
#' @include CycifStack-class.R
#'
#' @param filenames quantification file (.csv)
#' @param suffix a common suffix of each channel. If applicable, the suffix is
#'   removed from the channel names.
#' @param object a CycifStack object
#'
#' @rdname CycifStack
#' @export
CycifStack <- function(filenames,path,suffix="_cellRing"){
  stopifnot(all(file.exists(file.path(path,filenames))))
  cat("Trying to load ",length(filenames)," samples.\n",sep="")

  n.samples <- length(filenames)
  samples <- lapply(filenames,function(filename){
    if(grepl(suffix,filename)){
      name <- sub(paste0("unmicst-(.+)",suffix,"\\.csv"),"\\1",filename)
    }else{
      name <- sub(paste0("unmicst-(.+)\\.csv"),"\\1",filename)
    }
    cat("Loading ",name,"...\n",sep="")
    smpl <- Cycif(filename,path=path,suffix=suffix)
    return(smpl)
  })

  nms <- names(samples) <- sapply(samples,function(x){
      filename <- names(x)
      if(grepl(suffix,filename)){
        name <- sub(paste0("unmicst-(.+)",suffix,"\\.csv"),"\\1",filename)
      }else{
        name <- sub(paste0("unmicst-(.+)\\.csv"),"\\1",filename)
      }
  })
  names(nms) <- c()

  n_cells <- 	sapply(samples,nCells)
  n_cycles <- sapply(samples,nCycles)
  max_cycles <- max(n_cycles)

  if(0){
    ###### "Validate if samples with max_cycles use uniq_abs."
  }else{
    idx <- which(n_cycles == max_cycles)[1]
  }

  uniq_abs <- abs_list(samples[[idx]])
  new("CycifStack",
      names = nms,
      suffix = suffix,
      uniq_abs = uniq_abs,
      n_samples = n.samples,
      n_cycles = n_cycles,
      max_cycles = max_cycles,
      n_cells = n_cells,
      samples = samples
  )
}

#' @rdname CycifStack
#' @export
setMethod("show", "CycifStack", function(object) {
  df <- data.frame(
    SampleName = object@names,
    nCycles = object@n_cycles,
    nCells = object@n_cells)
  rownames(df) <- c()

  m <- matrix(object@uniq_abs$ab,nrow=3)
  rownames(m) <- c("R","G","B")
  colnames(m) <- paste0("Cycle",seq(object@max_cycles)-1)
  m <- data.frame(m)

  cat(is(object)[[1]], "\n")
  print(df)
  cat("Abs:\n")
  print(m)
})
