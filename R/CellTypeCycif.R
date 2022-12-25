#'@export
CellTypeCycif <- function(x,ctype,cstate,gates.df,ctype.full=FALSE){
  require(dplyr)
  if(class(x)=="Cycif"){
    abs <- abs_list(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be a Cycif object")
  }

  ## redefine ctype and cstate
  ctd <- CellTypeDefault(x,ctype,cstate,ctype.full=ctype.full)
  ctype <- ctd@cell_lineage_def
  cstate <- ctd@cell_state_def

  ## Subsetting ctype and cstate so only used antibodies exist in the experiment
  lmks <- names(ctype)[-c(1:2)]
  used.abs1 <- lmks[lmks %in% abs$ab]
  unused.abs1 <- lmks[!lmks %in% abs$ab]
  used.ctype1 <- ctype[,used.abs1,drop=F]
  unused.ctype1 <- ctype[,unused.abs1,drop=F]
  is.used.ct <- !apply(unused.ctype1=="AND",1,any) #

  smks <- colnames(cstate)
  used.abs2 <- smks[smks %in% abs$ab]

  ctype.sub <- ctype[is.used.ct,c("Parent","Child",used.abs1)]
  cstate.sub <- cstate[is.used.ct,used.abs2]

  ##
  mks <- unique(c(used.abs1,used.abs2))
  g <- rep(NA,length(mks))
  names(g) <- mks

  smpl <- names(x)
  g[] <- gates.df[mks,smpl]

  abs <- abs %>% filter(ab %in% mks)

  new("CellTypeCycif",
      name = x@name,
      n_cycles = x@n_cycles,
      cell_lineage_def = ctype.sub,
      cell_state_def = cstate.sub,
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
