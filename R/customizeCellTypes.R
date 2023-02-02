#'@export
customizeCellTypes <- function(ctype,uniq.cts,check=TRUE){
  pas <- ctype$Parent
  chs <- ctype$Child
  cts <- unique(c(pas,chs))
  
  if(!all(uniq.cts %in% cts)){
    non.exist <- uniq.cts[!uniq.cts %in% cts]
    stop("Some cell types are not defined in `ctype': ",paste0(non.exist,collapse=","))
  }
  
  ch.cts <- lapply(cts,function(x)x)
  names(ch.cts) <- cts
  ctlevs <- rev(CellTypeGraph(ctype,plot=F))
  
  for(lev in seq(ctlevs)){
    this.lev <- ctlevs[[lev]]
    ctype1 <- ctype %>% filter(Child %in% this.lev) %>% select(Parent,Child)
    lst1 <- tapply(ctype1$Child,ctype1$Parent,identity)
    for(pa in names(lst1)){
      lst1.ch <- unlist(ch.cts[c(lst1[[pa]],pa)])
      names(lst1.ch) <- c()
      ch.cts[[pa]] <- lst1.ch
    }
  }
  
  ##
  idx <- sapply(uniq.cts,function(nu){
    sapply(uniq.cts,function(de){
      is.under <- all(ch.cts[[nu]] %in% ch.cts[[de]])
    })
  })
  diag(idx) <- FALSE
  tab.idx <- which(idx,arr.ind=T)
  colnames(tab.idx) <- c("Parent","Child")
  rownames(tab.idx) <- c()
  hierarchy.cts <- as.data.frame(apply(tab.idx,2,function(i)uniq.cts[i]))
  if(any(idx)){
    cat("Provided uniq.cts include redundant cell types:\n")
    print(hierarchy.cts)
    stop("Provide non-redundant cell types.")
  }
  ch.cts1 <- ch.cts[uniq.cts]
  return(ch.cts1)
}

#'@export
switchCellTypes(x,ch.cts){
  return(ch.cts)
}
