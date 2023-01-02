#' Get stats for positive cells after gating each abs
#' @export
setGeneric("barplotPosCells", function(x,...) standardGeneric("barplotPosCells"))

#' @export
setMethod("barplotPosCells", "Cycif",function(x,type=c("one","two"),mar,...){
  smpl <- names(x)
  lth <- exprs(x,type="logTh_normalized")            
  ct <-  x@cell_type
  used.abs <- c(names(ct@cell_lineage_def)[-(1:2)],names(ct@cell_state_def))
  lin.abs <- colnames(ct@cell_lineage_def)[-(1:2)]

  pos <- as.data.frame(lapply(used.abs,function(ab){
    tmp <- lth[[ab]] > 0.5
    tmp <- factor(tmp,levels=c("TRUE","FALSE"))
    return(tmp)
  }))
  names(pos) <- used.abs
  
  is.str <- ct@is_strict
  
  opar <- par()
  if(type=="one"){
    diag <- sapply(used.abs,function(ab){
      pos1 <- pos[[ab]]
      if(ab %in% lin.abs){
        no <- sum(pos1=="FALSE",na.rm=T) # false
        yes.relax <- sum(pos1=="TRUE" & !is.str, na.rm=T) # false
        yes.str <- sum(pos1=="TRUE" & is.str, na.rm=T) # false      
        ntab <- c(0,yes.str,yes.relax,no)
        n <- ntab/sum(ntab)
      }else{
        ntab <- table(pos1)                
        ntab <- c(ntab[1],0,0,ntab[2])
        n <- ntab/sum(ntab)
      }
      return(n)
    })

    if(missing(mar)){
    par(mar=c(7,5,4,7))      
    }else{
      par(mar=mar)
    }

    barplot(diag,beside=F,las=2,col=c(3,2,4,"grey80"),main=smpl,border=FALSE,ylab="Ratio of positive cells",...)
    par(xpd=T)
    legend(par()$usr[2],par()$usr[4],c("strict","non-strict","cell-state"),fill=c(2,4,3))
    suppressWarnings(par(opar))
    invisible(diag)    
  }else if(type=="two"){
    ## compute overlap of positive cells
    rs <- sapply(used.abs,function(ab1){
      sapply(used.abs,function(ab2){
        a <- pos[[ab1]]
        b <- pos[[ab2]]
        idx <- !is.na(a) & !is.na(b)
        n <- sum(a[idx]=="TRUE" & b[idx]=="TRUE",na.rm=T)
        d <- sum(a[idx]=="TRUE",na.rm=T)
        r <- n/d
        return(r)
      })
    })
    
    dimnames(rs) <- list(
      numerator=rownames(rs),
      denominator=colnames(rs))
  
    rs <- round(rs,2)
    
    ##  heatmap 
    if(missing(mar)){
      par(mar=c(1,10,10,1))
    }else{
      par(mar=mar)
    }

    uniq.cols <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"))(50)

    
    par(fg="white")
    image(seq(nrow(rs)),seq(ncol(rs)),t(rs[rev(seq(ncol(rs))),]),col=NA,xlab="",ylab="",axes=F,zlim=c(0,1))
    axis(2,at=seq(nrow(rs)),labels=rev(rownames(rs)),las=1)
    axis(3,at=seq(nrow(rs)),labels=rownames(rs),las=2)
    par(fg="black")

    image(seq(nrow(rs)),seq(ncol(rs)),rs[,rev(seq(ncol(rs)))],col=uniq.cols,xlab="",ylab="",axes=F,add=T,
          main="",zlim=c(0,1),...)
    title(main=smpl,line=7)
    mtext("Numerator",side=3,line=5)
    mtext("Denominator",side=2,line=5)  
    suppressWarnings(par(opar))
    invisible(rs)
  }
})

#' @export
setMethod("barplotPosCells", "CycifStack",function(x,ab,mar,...){
  ct <-  x@cell_type
  used.abs <- c(names(ct@cell_lineage_def)[-(1:2)],names(ct@cell_state_def))
  lin.abs <- names(ct@cell_lineage_def)[-(1:2)]
  
  if(missing(ab)){
    stop("ab should be specified")
  }else if(!ab %in% used.abs){
    stop("ab should be one of the used abs")
  }
  
  lth <- exprs(x,type="logTh_normalized")
  sms <- sub("\\..+","",rownames(lth))

  tmp <- lth[[ab]] > 0.5
  lab <- rep(NA,length(tmp))
  if(ab %in% lin.abs){
    is.str <- ct@is_strict 
    
    lab[!tmp] <- "no"
    lab[tmp & !is.str] <- "pos.relax"
    lab[tmp & is.str] <- "pos.str"

    lab <- factor(lab, levels=c("pos.str","pos.relax","no"))
    
    leg <- c("strict","non-strict","negative")
    cols <- c(2,4,"grey80")
    mar <- c(7,5,7,7)
  }else{
    lab[!tmp] <- "neg"
    lab[tmp] <- "pos"
    lab <- factor(lab,levels=c("pos","neg"))

    leg <- c("positive","negative")
    cols <- c(3,"grey80")
    mar <- c(7,5,7,7)    
  }
  ntab <- table(lab,sms)
  diag <- apply(ntab,2,function(tmp)tmp/sum(tmp))
  
  opar <- par()

  par(mar=mar)
  
  xs <- barplot(diag,beside=F,las=2,main="",ylab="Ratio of positive cells",col=cols,...)
  title(main=ab,line=4)
  
  par(xpd=T)
  legend(par()$usr[2],par()$usr[4],leg,fill=cols)
  par(xpd=F)    

  # suppressWarnings(par(opar))
  invisible(xs)    
  
})
  