#' @export
CellTypeCalling <- function(cy,p_thres=0.5,strict=FALSE){
  # return a character vector containing 'cell_types'
  # cy <- x[[1]]
  lth <- exprs(x[[1]],type="logTh_normalized")
  if(nrow(lth)==0){
    stop("run normalize(method=\"logTh\") before CellTypeCalling()")
  }
  ctc <- cy@cell_type
  ctype <- ctc@cell_lineage_df
  ctlevs <- CellTypeGraph(ctype,plot=F)

  for(l in seq(length(ctlevs)-1)){
    pa <- ctlevs[[l]]
    ch <- ctlevs[[l+1]]
    ctype %>% filter(Parent %in% pa & Child %in% ch)
  }

  if(names(which.max(ctlevs))!="all"){
    stop("lineage_df: top level cell type should be 'all'\n")
  }

  for(l in uniq.levs[-1]){


  }


  uniq.ctypes <- ctype$Child
  ab.used <- colnames(ctype)[-(1:2)] %in% abs_list(cy)$ab

  ctype1 <- ctype[,ab.used]
  ctype.used <- apply(!is.na(ctype1),1,any)
  ctype1 <- ctype1[ctype.used,]
  norm <- exprs(cy,type="normalized")[colnames(ctype1)]

  used.cells <- apply(norm,1,function(x)all(!is.na(x)))
  n.cells <- sum(used.cells)

  ctdef <- apply(ctype1,1,function(x){
    thisdef <- x[!is.na(x)] =="AND"
  })

  def1 <- sapply(ctdef,function(def){
    cat("*")
    apply(norm >= p_thres,1,function(x){
      all(x[names(def)] == def)
    })
  })

  cts <- rep(NA,nrow(def1))
  i0 <- which(rowSums(def1)==0)
  i1 <- which(rowSums(def1)==1)
  ct1 <- apply(def1[i1,],1,function(x)names(which(x)))
  cts[i0] <- "unknown"
  cts[i1] <- ct1

  inc <- which(rowSums(def1)>1)
  if(strict){
    cts[inc] <- "inconc"
    uniq.cts <- c(rownames(ctype1),"unknown","inconc")
  }else{
    ct.abs <- unique(unlist(sapply(ctdef,function(x)names(which(x)))))
    inc.ct <- sapply(inc,function(i){
      tc <- this.ct <- ctdef[names(which(def1[i,]))]
      names(tc) <- c()
      this.mkrs <- names(which(unlist(tc)))
      this.mkrs <- this.mkrs[this.mkrs %in% ct.abs]
      max.mkr <- names(which.max(norm[i,this.mkrs]))
      ct1 <- names(this.ct)[which(sapply(this.ct,function(x)any(names(x)==max.mkr)))]
      return(ct1)
    })
    cts[inc] <- inc.ct
    uniq.cts <- c(rownames(ctype1),"unknown")
  }

  cts <- factor(cts,levels=uniq.cts)
  # is.def <- !is.na(rowSums(def1))
  names(cts) <- rownames(norm)
  return(cts)
}

