#' @include CellType-class.R

#' @export
setClass("CycifStack",
         slots = c(
           names = "character",
           suffix = "character",

           uniq_abs = "data.frame",

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

#' @export
CycifStack <- function(filenames){
  stopifnot(all(file.exists(filenames)))
  cat("Trying to load ",length(filenames)," samples.\n",sep="")

  n.samples <- length(filenames)
  samples <- lapply(filenames,function(filename){
    name <- sub("unmicst-(.+)\\.csv","\\1",filename)
    cat("Loading ",name,"...\n",sep="")
    smpl <- Cycif(filename)
    return(smpl)
  })

  suffix <- samples[[1]]@suffix
  nms <- 	names(samples) <- sapply(samples,names)
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
