#'@export
CellTypeCycifStack <- function(x,ctype,cstate,gates.df){
  require(dplyr)
  if(class(x)=="CycifStack"){
    abs <- abs_list(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be a CycifStack object")
  }

  gates.smpls <- names(gates.df)
  smpls <- names(x)
  if(!all(smpls %in% gates.smpls)){
    stop("All the samples should be gated before cell type calling")
  }

  nc <- nCycles(x)
  min.i <- min(which(nc==max(nc)))
  smpl <- names(x)[min.i]
  x1 <- x[[smpl]]

  ctc <- x1@cell_type

  lmks <- colnames(ctype)[-c(1:2)]
  smks <- colnames(cstate)

  mks <- c(lmks,smks)
  gates.list <- as.data.frame(cyApply(x,function(cy){
    ctc <- CellTypeCycif(cy,ctype,cstate,gates.df)
    ctc@gates[mks]
  },simplify=TRUE))

  is.ungated <- sapply(gates.list,function(g)all(is.na(g)))

  if(any(is.ungated)){
    warning("some samples are not gated:",paste(smpls[is.ungated],collapse=","))
  }

  new("CellTypeCycifStack",
      n_samples = x@n_samples,
      max_cycles = x@max_cycles,
      cell_lineage_def = ctc@cell_lineage_def,
      cell_state_def = ctc@cell_state_def,
      markers = ctc@markers,
      gates = gates.list
  )
}

setMethod("show", "CellTypeCycifStack", function(object){
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
  nstr <- length(nty$str)
  nst <- length(mst)

  is.gated <- sapply(object@gates,function(g){
    any(!is.na(g))
  })

  # cat(m)
  cat("[",is(object)[[1]],"]\n",
      "# samples:\t",object@n_samples,
      paste0("(",sum(is.gated)," gated)"),"\n",
      "# max cycles:\t",object@max_cycles,"\n\n",
      "# cell types:\t",nct,"\n",
      paste(cts,collapse=", "),"\n\n",

      "# markers in total:\t", nmk,"\n",
      "# cell lineage markers:",nlin,"\n",
      "# stratifying markers:\t",nstr,"\n",
      "# cell state markers:\t",nst,"\n\n",

      "lineage markers:",
      paste(nty$lin,collapse=", "),"\n",
      "stratifying markers:",
      paste(nty$str,collapse=", "),"\n",
      "state markers:",
      paste(mst,collapse=", ")
  )
})
