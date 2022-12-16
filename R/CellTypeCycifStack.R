CellTypeCycifStack <- function(x,lineage_df,state_df,gates.df){
  require(dplyr)
  if(class(x)=="CycifStack"){
    abs <- uniq_abs(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be a CycifStack object")
  }
  
  smpls <- names(x)
  
  n <- length(x)
  nc <- nCycles(x)
  min.i <- min(which(nc==max(nc)))
  smpl <- names(x)[min.i]
  x1 <- x[[smpl]]
  
  ctc <- CellTypeCycif(x1,lineage_df,state_df,gates.df)
  
  lmks <- colnames(lineage_df)[-c(1:2)]
  smks <- colnames(state_df)
  
  mks <- c(lmks,smks)
  gs <- data.frame(array(NA,dim=c(length(mks),nSamples(x)),dimnames=list(mks,names(x))))
  names(gs) <- sub("^X","",names(gs))
  for(ncy in smpls){
    cy <- x[[ncy]]
    used.mks <- mks[mks %in% abs_list(cy)$ab]
    gs[used.mks,ncy] <- gates.df[used.mks,ncy]
  }
  
  ungated <- sapply(gs,function(g){
    all(is.na(g))
  })
  if(any(ungated)){
    warning("some samples are not gated:",paste(names(gs)[ungated],collapse=","))
  }
  
  new("CellTypeCycifStack",
      n_samples = x@n_samples,
      max_cycles = x@max_cycles,
      cell_lineage_df = ctc@cell_lineage_df,
      cell_state_df = ctc@cell_state_df,
      markers = ctc@markers,
      gates = gs
  )
}

setMethod("show", "CellTypeCycifStack", function(object){
  nmk <- length(object@markers$ab)
  cts <- object@cell_lineage_df$Child
  cts <- cts[!cts == "unknown"]
  nct <- length(cts)
  mty <- apply(object@cell_lineage_df[-c(1:2)],2,function(x){
    if(all(x %in% c("AND","OR","NOT","NOR",""))){
      return("lin")
    }else if(all(x %in% c("CAN",""))){
      return("str")
    }else{
      return("others")
    }
  })
  nty <- tapply(names(mty),mty,identity)
  mst <- colnames(object@cell_state_df)
  
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