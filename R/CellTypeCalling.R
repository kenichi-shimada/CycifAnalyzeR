#' @export
CellTypeCalling <- function(cy,ctype,p_thres=0.5,strict=FALSE){
  # cy <- x[[1]]

  ctlevs <- CellTypeGraph(cy,ctype,plot=F)
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
  if(0){
    used.abs <- used_abs(cy)
    is.def <- sapply(ctdef,function(x){
      all(names(x) %in% used.abs)
    })
  }

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
# setGeneric("defCellType", function(x,...) standardGeneric("defCellType"))
# setMethod("defCellType", "Cycif",
#
#   # this is not done yet.
#
# 	function(x,ctype.def,trim=1e-5,method="logTh",...){
# 		## default method is logTh
# 		if(missing(method)){
# 			method <- "logTh"
# 		}
#
# 		## 'used abs' set already?
# 		if("used_abs" %in% slotNames(x)){
# 			used.abs <- used_abs(x)
# 		}else{
# 			used.abs <- as.character(abs_list(x)$ab)
# 		}
#
# 		smpl <- names(x)
# 		raw <- x@raw
# 		is.used <- cumUsedCells(x)
# 		x@normalize.method <- method
#
# 		## treatment is different between methods
# 		if(method=="log"){
# 			norm <- as.data.frame(
# 				sapply(used.abs,function(ab){
# 					cycle <- abs_list(x)$cycle[abs_list(x)$ab==ab]
# 					cycle <- cycle + 1 # 0-origin
# 					is.used.1 <- is.used[,cycle]
#
# 					r <- raw[[ab]]
# 					n <- rep(NA,length(r))
#
# 					n[is.used.1] <- transform(r[is.used.1],method="log",trim=trim)
# 					return(n)
# 				})
# 			)
# 		}else if(method=="logTh"){
# 			if(!.hasSlot(x,"threshold")){
# 				stop(smpl,": set threshold first.\n")
# 		}
#
# 			thres <- threshold(x)
# 			used.abs <- used.abs[!is.na(thres[used.abs])]
#
# 			norm <- as.data.frame(
# 				sapply(used.abs,function(ab){
# 					cycle <- abs_list(x)$cycle[abs_list(x)$ab==ab]
# 					cycle <- cycle + 1 # 0-origin
# 					is.used.1 <- is.used[,cycle]
# 					r <- raw[[ab]]
# 					n <- rep(NA,length(r))
# 					th <- thres[ab]
# 					n[is.used.1] <- transform(r[is.used.1],method="logTh",th=th,trim=trim)
# 					return(n)
# 				})
# 			)
# 		}
# 		x@normalized <- norm
# 		return(x)
# 	}
# )
