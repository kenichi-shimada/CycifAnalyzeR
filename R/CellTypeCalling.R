#' @export
CellTypeCalling <- function(cy,p_thres=0.5,strict=FALSE,expanded_df=TRUE){
  # return a character vector containing 'cell_types'
  # cy <- x[[1]]
  lth <- exprs(cy,type="logTh_normalized")
  if (nrow(lth)==0) {
    stop("run normalize(method=\"logTh\") before CellTypeCalling()")
  }
  if(strict){
    stop("strict=TRUE is not implemented yet")
  }

  ctc <- cy@cell_type
  if(expanded_df){
    ctype <- ctc@expanded_lineage_df
  }else{
    ctype <- ctc@cell_lineage_df
  }

  ctlevs <- CellTypeGraph(ctype,plot=F)

  cell.types <- rep("all",nrow(lth))
  for(l in seq(length(ctlevs)-1)[1]){
    #l=1
    pas <- ctlevs[[l]]
    chs <- ctlevs[[l+1]]
    i.others <- chs %in% "unknown" | grepl("_other",chs)
    chs1 <- chs[!i.others]

    ct <- ctype %>% filter(Parent %in% pas & Child %in% chs1)
    uniq.pas <- unique(ct$Parent)
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
          stop(ch," should have either AND or OR in the definition")
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

    # i0 <- which(cell.types=="dropped")
    for(pa in uniq.pas){
      this.ind <- cell.types==pa
      this.chs <- (ct %>% filter(Parent==pa))$Child
      this.chs <- colnames(prs)[colnames(prs) %in% this.chs]
      cell.types[this.ind] <- apply(prs[this.ind,this.chs,drop=F],1,function(pr){
        if(all(is.na(pr))){
          return(pa)
        }
        ind <- which(pr==max(pr,na.rm=T))
        if(length(ind)>1){
          return(pa)
        }
        if(pr[ind] > p_thres){
          return(colnames(prs)[ind])
        }else{
          ch1 <- paste0(pa,"_other")
          if(ch1 == "all_other"){
            ch1 <- "unknown"
          }
          return(ch1)
        }
      })
    }
  }
  ## convert to factor: Q: how to order cell types from ctype df?
  return(cell.types)
}

