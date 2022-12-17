CellTypeCycif <- function(x,lineage_df,state_df,gates.df){
  require(dplyr)
  if(class(x)=="Cycif"){
    abs <- abs_list(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be a Cycif object")
  }

  lmks <- colnames(lineage_df)[-c(1:2)]
  used.abs1 <- lmks[lmks %in% abs$ab]
  unused.abs1 <- lmks[!lmks %in% abs$ab]
  used1 <- lineage_df[,used.abs1]
  unused1 <- lineage_df[,unused.abs1]
  used.cts <- !apply(unused1=="AND",1,any)

  smks <- colnames(state_df)
  used.abs2 <- smks[smks %in% abs$ab]

  lineage_df.sub <- cbind(lineage_df[1:2],lineage_df[used.abs1])[used.cts,]
  state_df.sub <- state_df[used.cts,used.abs2]

  ctd <- CellTypeDefault(x,lineage_df.sub,state_df.sub)

  mks <- unique(c(used.abs1,used.abs2))
  g <- rep(NA,length(mks))
  names(g) <- mks

  smpl <- names(x)
  g[] <- gates.df[mks,smpl]

  new("CellTypeCycif",
      name = x@name,
      n_cycles = x@n_cycles,
      cell_lineage_df = ctd@cell_lineage_df,
      cell_state_df = ctd@cell_state_df,
      expanded_lineage_df = ctd@expanded_lineage_df,
      expanded_state_df = ctd@expanded_state_df,
      markers = ctd@markers,
      gates = g
  )
}

setMethod("show", "CellTypeCycif", function(object){
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
      "Sample name:\t",object@name,"\n",
      "# cycles:\t",object@n_cycles,"\n\n",
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
