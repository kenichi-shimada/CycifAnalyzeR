#' @export
CellTypeCalling <- function(cy,p_thres=0.5,ctype.full=TRUE){
  # return a character vector containing 'cell_types'
  # cy <- x[[1]]
  lth <- exprs(cy,type="logTh_normalized")
  if (nrow(lth)==0) {
    stop("run normalize(method=\"logTh\") before CellTypeCalling()")
  }

  if(ctype.full){
    ctc <- cy@cell_type_full
  }else{
    ctc <- cy@cell_type
  }
  ctype <- ctc@cell_lineage_def

  ctlevs <- CellTypeGraph(ctype,plot=F)

  cell.types <- rep("all",nrow(lth))
  is.strict <- rep(TRUE,nrow(lth))

  for(l in seq(length(ctlevs)-1)){
    pas <- ctlevs[[l]]
    chs <- ctlevs[[l+1]]
    i.others <- chs %in% "unknown" | grepl("_other$",chs)
    chs1 <- chs[!i.others]

    ct <- ctype %>% filter(Parent %in% pas & Child %in% chs1)
    uniq.pas <- unique(ct$Parent)

    ## compute probability of being each child node in the following block
    ## (not considering the parent node)
    prs <- sapply(chs1,function(ch){
      tmp <- unlist(ct %>% filter(Child == ch))
      pa <- tmp[1]
      tmp <- tmp[-(1:2)]
      abs.and <- names(which(tmp=="AND"))
      abs.or <- names(which(tmp=="OR"))
      abs.not <- names(which(tmp=="NOT"))
      suppressWarnings({
        if(length(abs.and)>0 & length(abs.or)>0){
          a <- apply(lth[abs.and],1,min,na.rm=T)
          a[a==-Inf] <- NA
          b <- apply(lth[abs.or],1,max,na.rm=T)
          b[b==-Inf] <- NA
          pos.out <- pmin(a,b)
        }else if(length(abs.and)>0){
          pos.out <- apply(lth[abs.and],1,min,na.rm=T)
          pos.out[pos.out==Inf] <- NA
        }else if(length(abs.or)>0){
          pos.out <- apply(lth[abs.or],1,max,na.rm=T)
          pos.out[pos.out==-Inf] <- NA
        }else{
          pos.out <- rep(1,nrow(lth))
        }

        this.ct <- pos.out
        if(length(abs.not)>0){
          neg.out <- apply(lth[abs.not],1,max,na.rm=T)
          neg.out[neg.out==-Inf] <- NA
          this.ct[which(neg.out > p_thres)] <- 0
          this.ct[is.na(neg.out)] <- NA
        }
      })
      return(this.ct)
    })

    ## for each parent, show which child is more likely to be the cell type each cell is
    for(pa in uniq.pas){
      this.ind <- cell.types==pa
      this.chs <- (ct %>% filter(Parent==pa))$Child
      this.chs <- colnames(prs)[colnames(prs) %in% this.chs]
      prs1 <- prs[this.ind,this.chs,drop=F]
      cell.types[this.ind] <- apply(prs1,1,function(pr){
          if(all(is.na(pr))){
            return(pa)
          }
          ind <- which(pr==max(pr,na.rm=T))
          if(length(ind)>1){
            return(pa)
          }
          if(pr[ind] > p_thres){
            return(this.chs[ind])
          }else{
            ch1 <- paste0(pa,"_other")
            if(ch1 == "all_other"){
              ch1 <- "unknown"
            }
            return(ch1)
          }
      })
      this.strict <- rowSums(prs1 > p_thres,na.rm=F) < 2
      this.strict[is.na(this.strict)] <- FALSE
      is.strict[this.ind] <- this.strict & is.strict[this.ind]
    }
  }

  ## convert to factor: Q: how to order cell types from ctype df?
  uniq.cts <- c("all",ctype$Child)
  leaves <- uniq.cts[!uniq.cts %in% ctype$Parent]

  uniq.cts <- c(uniq.cts[uniq.cts!="unknown"],"unknown")
  cts <- factor(cell.types,levels=uniq.cts)

  return(list(cell_type=cts,is_strict=is.strict))
}

