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

  lmks <- names(ctype)[-c(1:2)]
  smks <- colnames(cstate)

  mks <- unique(c(lmks,smks))
  if(!all(mks %in% used.abs)){
    ## here ctype and cstate should be subsetted based on available mks
    ## Subsetting ctype and cstate so only used antibodies exist in the experiment
    used.abs1 <- lmks[lmks %in% abs$ab]
    unused.abs1 <- lmks[!lmks %in% abs$ab]
    used.ctype1 <- ctype[,used.abs1,drop=F]
    unused.ctype1 <- ctype[,unused.abs1,drop=F]
    is.used.ct <- !apply(unused.ctype1=="AND",1,any) #

    used.abs2 <- smks[smks %in% abs$ab]

    ctype <- ctype[is.used.ct,c("Parent","Child",used.abs1)]
    cstate <- cstate[is.used.ct,used.abs2]

    ## up to here
    unknown <- mks[!mks %in% used.abs]
  }
  mks.info <- abs %>% filter(ab %in% mks)

  elin <- expandLineageDef(ctype=ctype,cstate=cstate,ctype.full=ctype.full)

  new("CellTypeDefault",
      cell_lineage_def = elin$ctype,
      cell_state_def = elin$cstate,
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
