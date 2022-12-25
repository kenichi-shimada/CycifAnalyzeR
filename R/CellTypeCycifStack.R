#'@export
CellTypeCycifStack <- function(x,ctype.full=FALSE){
  require(dplyr)
  if(class(x)=="CycifStack"){
    abs <- abs_list(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be a CycifStack object")
  }

  nc <- nCycles(x)
  min.i <- min(which(nc==max(nc)))
  smpl <- names(x)[min.i]
  x1 <- x[[smpl]]

  if(ctype.full){
    ctc <- x1@cell_type_full
  }else{
    ctc <- x1@cell_type
  }

  ## load ctype and cstate
  ctype <- ctc@cell_lineage_def
  cstate <- ctc@cell_state_def

  ## marker abs
  lmks <- colnames(ctype)[-c(1:2)]
  smks <- colnames(cstate)
  mks <- c(lmks,smks)

  ## compile gates
  gates.df <- as.data.frame(cyApply(x,function(cy){
    if(ctype.full){
      ctc <- cy@cell_type_full
    }else{
      ctc <- cy@cell_type
    }
    ctc@gates[mks]
  },simplify=TRUE))

  gates.smpls <- names(gates.df)
  smpls <- names(x)
  if(!all(smpls %in% gates.smpls)){
    stop("Check gates; All the samples should be gated before cell type calling")
  }

  is.ungated <- sapply(gates.df,function(g)all(is.na(g)))

  if(any(is.ungated)){
    warning("some samples are not gated:",paste(smpls[is.ungated],collapse=","))
  }

  abs <- abs %>% filter(ab %in% mks)

  new("CellTypeCycifStack",
      n_samples = x@n_samples,
      max_cycles = x@max_cycles,
      cell_lineage_def = ctype,
      cell_state_def = cstate,
      markers = abs,
      gates = gates.df
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
