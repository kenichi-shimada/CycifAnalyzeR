#'@export
CellTypeCycif <- function(x,ctype,cstate,gates.df,ctype.full=FALSE){
  require(dplyr)
  if(class(x)=="Cycif"){
    abs <- abs_list(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be a Cycif object")
  }

  if(is.matrix(gates.df)){
    gates.df <- as.data.frame(gates.df)
  }
  ## redefine ctype and cstate
  ctd <- CellTypeDefault(x,ctype,cstate,ctype.full=ctype.full)
  ctype <- ctd@cell_lineage_def
  cstate <- ctd@cell_state_def
  mks <- ctd@markers$ab

  ##
  g <- rep(NA,length(mks))
  names(g) <- mks

  smpl <- names(x)
  g[] <- gates.df[mks,smpl]

  abs <- abs %>% filter(ab %in% mks)

  new("CellTypeCycif",
      name = x@name,
      n_cycles = x@n_cycles,
      cell_lineage_def = ctype,
      cell_state_def = cstate,
      markers = abs,
      gates = g
  )
}

setMethod("show", "CellTypeCycif", function(object){
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
