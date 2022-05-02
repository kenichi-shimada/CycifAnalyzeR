#'@export
defCellType_v2 <- function(x,ctype,defs,depths,p_thres=0.4){
  # when l = 2, I didn't check if the tumors' are exclusively
  # classified into subclasses. That check should be done (To be implemented)
  # x = lt
  # p_thres=0.4
  mat <- exprs(x,type="normalized")
  mat <- mat[!grepl("smpl",names(mat))]
  # dim(mat)
  is.used.cells <- rowSums(is.na(mat))==0
  mat.1 <- mat[is.used.cells,]
  ctypes.all <- rep(NA,nrow(mat))
  names(ctypes.all) <- rownames(mat)

  ctypes <- rep("na",nrow(mat.1))

  for(l in as.numeric(levels(depths))){
    cat("l = ",l,"\n")
    if(l==1){
      ## depth.levs ==1
      ## check if each cell can be one or more of the phenotype
      cts1 <- names(depths)[which(depths==l)]

      idx1 <- rep(cts1,sapply(defs[cts1],length))
      abs1 <- unlist(lapply(defs[cts1],names))

      is.pos <- apply(mat.1[abs1],1,max) > p_thres
      max.ind <- idx1[apply(mat.1[abs1],1,which.max)]
      ctypes[!is.pos] <- "unknown"
      ctypes[is.pos] <- max.ind[is.pos]
      ctype1 <- ctypes
      tab <- table(ctype1)
      for(k in seq(tab)){
        cat(" ",names(tab)[k]," : ",tab[k],"\n")
      }
    }else if(l>1){
      cts1 <- depth.levs[[l-1]]
      for(ct1 in cts1){
        #ct1 <- cts1[1]
        cts2 <- ctype[ctype[,"parent"]==ct1,"child"]
        #cat(length(cts2),"\n")
        for(ct2 in cts2){
          #ct2 <- cts2[1]
          cat(" ",ct1," => ", ct2," : ")
          def <- defs[[ct2]]
          if(any(def=="AND")){
            this.prof <- def[def %in% c("AND","NOT")] # ignore "OR"
            this.prof <- this.prof=="AND"
            is.this.ct <- apply(mat.1[names(this.prof)] > p_thres,1,function(x){
              all(x==this.prof) # check both AND and NOT
            })
          }else if(any(def=="OR")){
            this.prof <- def[def %in% c("OR","NOT")] # shouldn't remove anything
            cells.or <- names(which(this.prof=="OR"))
            ## at least one pos for OR and all negative for NOT
            has.pos <- rowSums(mat.1[cells.or] > p_thres) > 0 # check OR
            if(any(def=="NOT")){
              cells.not <- names(which(this.prof=="NOT"))
              has.neg <- rowSums(mat.1[cells.not] > p_thres)==0
              is.this.ct <- has.pos & has.neg
            }else{
              is.this.ct <- has.pos
            }
          }else{
            stop("def should contain at least one \"AND\" or \"OR\"")
          }
          cat(sum(ctype1==ct1 & is.this.ct),"\n")
          ctypes[ctype1==ct1 & is.this.ct] <- ct2
        }
      }
      ctype1 <- ctypes
    }
  }
  ctypes.all[is.used.cells] <- ctypes
  return(ctypes.all)
}
