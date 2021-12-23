def2num <- function(ctype){
  is.and <- ctype == "AND"
  is.can <- ctype == "CAN"
  num.can <- num.and <- array(0,dim(is.and))

  num.and[is.and] <- 1
  num.can[is.can] <- 1

  uniq.bins <- 2^(seq(ncol(ctype))-1)
  idx <- lapply(seq(nrow(ctype)),function(i){
    and <- (num.and %*% uniq.bins)[i]
    idx <- which(num.can[i,]==1)
    can.mat <- cbind(0,uniq.bins[idx])
    cans <- rowSums(expand.grid(lapply(seq(nrow(can.mat)),function(j){
      can.mat[j,]
    })))
    all <- and + cans
    return(all)
  })
  names(idx) <- rownames(ctype)
  vec.idx <- unlist(idx)
  cell.names <- rep(names(idx),sapply(idx,length))
  df <- data.frame(idx=vec.idx,cell.name=cell.names,stringsAsFactors=F)
  df <- rbind(c(0,"others"),df)
  rownames(df) <- c()
  return(df)
}



