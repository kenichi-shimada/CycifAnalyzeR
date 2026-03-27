#_ -------------------------------------------------------

# fun: show Cycif ----
#' @rdname Cycif
#' @export
setMethod("show", "Cycif", function(object) {
  n.ch <- max(table(abs_list(object)$cycle))
  m <- do.call(cbind,tapply(abs_list(object)$ab,abs_list(object)$cycle,
              function(x){
                tmp <- rep(NA,n.ch)
                tmp[seq(x)] <- x
                return(tmp)
  }))

  rownames(m) <- paste0("Ch",seq(n.ch))
  colnames(m) <- paste0("Cycle",seq(nCycles(object)))
  m <- data.frame(m)

  nd <- nrow(dna_thres(object)) > 0
  nln <- nrow(object@log_normalized) > 0
  nltn <- nrow(object@logTh_normalized) > 0
  nct <- length(object@cell_types) > 0
  cat("[",is(object)[[1]], " object]\n\n",
      "Name: ", object@name, "\n",
      "nCells:  ", object@n_cells, "\n",
      "nCycles: ", object@n_cycles, "\n\n",
      "DNA filtered: ", nd, "\n",
      "Log-normalized: ", nln, "\n",
      "LogTh-normalized: ", nltn, "\n",
      "CellType-idenitified: ", nct, "\n\n",
      "Abs: \n",sep="")
  print(m)
})


#_ -------------------------------------------------------

# fun: show CycifStack ----
#' @rdname CycifStack
#' @export
setMethod("show", "CycifStack", function(object) {
  df <- data.frame(
    SampleName = object@names,
    nCycles = object@n_cycles,
    nCells = object@n_cells)
  rownames(df) <- c()

  n.ch <- max(table(abs_list(object)$cycle))
  m <- do.call(cbind,tapply(abs_list(object)$ab,abs_list(object)$cycle,
                            function(x){
                              tmp <- rep(NA,n.ch)
                              tmp[seq(x)] <- x
                              return(tmp)
                            }))

  max.cycles <- max(df$nCycles)
  rownames(m) <- paste0("Ch",seq(n.ch))
  colnames(m) <- paste0("Cycle",seq(max.cycles))
  m <- data.frame(m)

  cat(is(object)[[1]], "\n")
  print(df)
  cat("Abs:\n")
  print(m)
})
