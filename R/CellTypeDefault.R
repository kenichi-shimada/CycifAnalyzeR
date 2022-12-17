CellTypeDefault <- function(x,lineage_df,state_df){
  require(dplyr)
  if(class(x) %in% c("Cycif","CycifStack")){
    abs <- abs_list(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be either a Cycif or CycifStack object")
  }
  if(missing(lineage_df) || missing(state_df)){
    stop("both lineage and state definitions should be provided")
  }
  if(class(lineage_df) != "data.frame"){
    stop("cell lineage definition should be a data.frame")
  }
  if(class(state_df) != "data.frame"){
    stop("cell state definition should be a data.frame")
  }
  if(!all(as.character(lineage_df$Child)==rownames(state_df))){
    stop("lineage_df$Child and cell types in state_df should be identical")
  }
  mks <- unique(c(colnames(lineage_df)[-(1:2)],colnames(state_df)))
  if(!all(mks %in% used.abs)){
    unknown <- mks[!mks %in% used.abs]
    stop(paste0("some markers in the defs not found: ",paste(unknown,collapse=",")))
  }else{
    mks.info <- abs %>% filter(ab %in% mks)
  }

  elin <- expandLineageDef(lineage_df)

  est <- state_df[elin$names$idx,]
  rownames(est) <- elin$names$expanded

  new("CellTypeDefault",
      cell_lineage_df = lineage_df,
      cell_state_df = state_df,
      expanded_lineage_df = elin$lineage_df,
      expanded_state_df = est,
      markers = mks.info
  )
}

setMethod("show", "CellTypeDefault", function(object){
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
