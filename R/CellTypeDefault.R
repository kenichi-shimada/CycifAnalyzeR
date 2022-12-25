CellTypeDefault <- function(x,ctype,cstate,ctype.full=TRUE){
  require(dplyr)
  if(class(x) %in% c("Cycif","CycifStack")){
    abs <- abs_list(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be either a Cycif or CycifStack object")
  }
  if(missing(ctype) || missing(cstate)){
    stop("both lineage and state definitions should be provided")
  }
  if(class(ctype) != "data.frame"){
    stop("cell lineage definition should be a data.frame")
  }
  if(class(cstate) != "data.frame"){
    stop("cell state definition should be a data.frame")
  }
  if(!all(as.character(ctype$Child)==rownames(cstate))){
    stop("ctype$Child and cell types in cstate should be identical")
  }
  mks <- unique(c(colnames(ctype)[-(1:2)],colnames(cstate)))
  if(!all(mks %in% used.abs)){
    unknown <- mks[!mks %in% used.abs]
    stop(paste0("some markers in the defs not found: ",paste(unknown,collapse=",")))
  }else{
    mks.info <- abs %>% filter(ab %in% mks)
  }

  elin <- expandLineageDef(ctype,ctype.full=ctype.full)

  est <- cstate[elin$names$idx,]
    rownames(est) <- elin$names$expanded

  new("CellTypeDefault",
      cell_lineage_def = ctype,
      cell_state_def = cstate,
      markers = mks.info
  )
}

setMethod("show", "CellTypeDefault", function(object){
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

  # cat(m)
  cat("[",is(object)[[1]],"]\n",
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
