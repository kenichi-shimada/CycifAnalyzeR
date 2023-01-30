#' Filtering out unavaialble cells based on nuclei staining
#'
#' @param x A Cycif object

#' @param manual logical. If TRUE, lower and upper limits of DNA intensity should be chosen
#'   interactively by clicking the limit values in the heatmap. FALSE by default.
#' @param ratio logical. If TRUE, the ratio of DNA intensity between each cycle and cycle 0 will be used.
#'   If FALSE, raw values of DNA intensity will be used.
#' @param n numeric. The number of breaks in the histogram.
#' @param n1 numeric. indices from 1 through n1 will be plugged into loess().
#' @param show.only logical. If TRUE, only show the summary of the plots but don't attempt to filter cells.
#'
#' @export
setGeneric("dnaFilter", function(x,...) standardGeneric("dnaFilter"))
setMethod("dnaFilter", "Cycif",
  function(x,show.only=FALSE){
    mat <- x@dna
    smpl <- x@name
    n <- n1 <- 1000

    dna.list <- names(mat)
    mat <- cbind(log1p(mat[[1]]),as.data.frame(lapply(mat,function(x)log1p(x/mat[[1]])))[-1])
    names(mat) <- dna.list

    dna.thres <- x@dna_thres

    if(nrow(dna.thres)>0){
      dna.ths1 <- x@dna_thres$low
      dna.ths2 <- x@dna_thres$high
      names(dna.ths1) <- names(dna.ths2) <- dna.list
    }else{
      dna.ths1 <- sapply(mat,min)
      dna.ths2 <- sapply(mat,max)
      names(dna.ths1) <- names(dna.ths2) <- dna.list
    }

    for(i in seq(length(dna.list))){
        channel <- dna.list[i]
        m <- mat[[channel]]

        xmax <- max(m)

        if(i==1){
          ttl <- channel
          brks <- seq(min(m),max(m),length=(n+1))
        }else{
          ttl <- paste0(channel," / ",dna.list[1])
          brks <- seq(0,xmax,length=(n+1))
        }

        l <- layout(matrix(c(2,1),nrow=2),heights=c(2,3))
        slidePlot(x,plot_type="dna",ab=channel,mar=c(3,3,0,3),ttl="",use.roi=FALSE)
        hst <- hist_fun(x=m,ths=c(dna.ths1[i],dna.ths2[i]),brks1=brks,ttl1=ttl)

        smoothened <- hst$smoothened
        a <- hst$a

        ##
        if(!show.only){
            cat(paste0("Filtering ",smpl,"...\n"))
            # cat("Do you want to see the auto_filters?")
            # ans <- readline(prompt="Y/N [N]:")
            auto_filter <- FALSE

            if(auto_filter){
                ## auto_filter thresholding
                th.up <- quantile(m,.9)

                ind.v <- which(diff(sign(diff(smoothened)))>0  & a$mids[c(-1,-n1)] < th.up) # idx of valleys
                ind.p <- which(diff(sign(diff(smoothened)))<0  & a$mids[c(-1,-n1)] < th.up) # idx of peaks
                idx <- sort(c(ind.v,ind.p))

                ## if idx is left to idx.lo, remove them (loess gets bumpy at the edge)
                idx.diff.th <- 13

                while(any(diff(idx)<= idx.diff.th)){
                  k <- which.min(diff(idx))
                  idx <- idx[-c(k,k+1)]
                  ind.v <- ind.v[ind.v %in% idx]
                  ind.p <- ind.p[ind.p %in% idx]
                }

                ## identify removed cells
                if(length(ind.v)>0){
                  ind.v <- ind.v
                  tmp <- a$mids[max(ind.v)]
                  if(tmp < 0.25){
                    x.drop.th <- tmp
                  }else{
                    x.drop.th <- 0
                  }
                }else{
                  x.drop.th <- 0
                }

                dna.ths1[i] <- max(min(m),x.drop.th)

                ## identify bunching cells
                idx.hi <- 40

                if(length(ind.p)>0){
                  if(max(ind.p) > idx.hi){
                    p.max <- max(ind.p)
                    x.max <- a$mids[p.max]
                    y.max <- a$density[p.max]
                    dist.half.max <- which.min(abs(a$density - y.max/2)[-seq(p.max)])
                    x.half.max <- dist.half.max
                    ind.bu <- p.max + x.half.max*3
                    if(ind.bu < n){
                      x.bunch.th <- a$mids[ind.bu]
                    }else{
                      x.bunch.th <- Inf
                    }
                    abline(v=x.max,col=1)
                    abline(v=x.half.max,col=1)
                    abline(v=x.bunch.th,col=2)
                  }else{
                    x.bunch.th <- Inf
                  }
                }

                dna.ths2[i] <- min(max(m),x.bunch.th)

            }
            if(0){
              ## show current dna.ths1 and dna.ths2
              in.rng <- factor((m > dna.ths1[i]) + (m > dna.ths2[i]) + 1,levels=c(2,1,3))
              cexs <- c(1,20)[(in.rng %in% as.character(c(1,3)))+1]

              slidePlot(x,plot_type="filter",within_filter_rng=in.rng,
                        uniq.cols=c("grey80",4,2),
                        cex=2,
                        mar=c(3,3,0,3),ttl="")
              hist_fun(x=m,ths=c(dna.ths1[i],dna.ths2[i]),brks1=brks,ttl1=ttl)
            }

            # if(manual){
            if(TRUE){
                cat("How do you want to modify the 'dna.ths'?\n")
                ans <- "init"
                lo <- dna.ths1[i]
                hi <- dna.ths2[i]
                while(ans!="0"){
                  ans <- readline(prompt="0:none, 1:lower th, 2:higher th, 3:both ths, 4:check fl dist, 5: check current ths [0-5]")
                  ans <- sub("^(.).*","\\1",ans)
                  if(ans %in% as.character(c(1:3,5))){
                    if(ans=="1"){
                      lo <- locator(1)$x
                      abline(v=lo,col=4)
                    }else if(ans=="2"){
                      hi <- locator(1)$x
                      abline(v=hi,col=2)
                    }else if(ans=="3"){
                      bdr <- locator(2)$x
                      bdr <- sort(bdr,decreasing=F)
                      lo <- bdr[1]
                      hi <- bdr[2]
                      abline(v=lo,col=4)
                      abline(v=hi,col=2)
                    }

                    in.rng <- factor((m > lo) + (m > hi) + 1,levels=c(2,1,3))
                    cexs <- c(1,20)[(in.rng %in% as.character(c(1,3)))+1]

                    ## update
                    slidePlot(x,plot_type="filter",within_filter_rng=in.rng,
                              uniq.cols=c("grey80",4,2),
                              cex=2,
                              mar=c(3,3,0,3),ttl="")
                    hist_fun(x=m,ths=c(lo,hi),brks1=brks,ttl1=ttl)
                  }else if(ans=="4"){
                    slidePlot(x,plot_type="dna",ab=channel,mar=c(3,3,0,3),ttl="")
                    hist_fun(x=m,ths=c(lo,hi),brks1=brks,ttl1=ttl)
                  }
                }
                dna.ths1[i] <- lo
                dna.ths2[i] <- hi
            }
        }
    }

    used.cells <- sapply(names(mat),function(channel){
      ind <- (mat[[channel]] > dna.ths1[channel]) + (mat[[channel]] > dna.ths2[channel])
      return(ind)
    })

    ## summarise used.cells
    uc <- used.cells==1
    ucs <- sapply(seq(ncol(used.cells)),function(i){
      uc1 <- rep(1,nrow(uc))
      for(j in seq(i)){
        uc1 <- uc1 * uc[,j]
      }
      return(uc1)
    })
    ret <- rowSums(ucs)

    x@dna_thres <- data.frame(low=dna.ths1,high=dna.ths2)
    x@used_cells <- used.cells

    ## final after dnaFitlter
    nc <- length(unique(ret))

    uniq.cols <- colorRampPalette(brewer.pal(9,"YlGnBu"))(nc+2)[-(nc+(0:1))]

    cat("Cell retention through each cycle:\n")
    l <- layout(matrix(c(1,2),nrow=2),heights=c(2,3))
    plotUsedCellRatio(x)
    slidePlot(x,plot_type="filter",within_filter_rng=ret,
              uniq.cols=uniq.cols,
              cex=2,ncells=1e4,
              mar=c(3,3,0,3),ttl="")

    # ## choose ROI
    # cat("Do you want to select ROIs?\n")
    # ans <- readline(prompt="Y/N [Y]")
    # pos.rois <- list()
    # if(!grepl("^[nN]",ans)){
    #   cat("Positive ROIs, negative ROIs, or Cancel?")
    #   ans <- readline(prompt="P/N/C []")
    #   while(!grepl("[Cc]",ans)){
    #     if(grepl("^[Pp]",ans)){
    #       x <- roiFilter(x,roi_type="positive")
    #     }else if(grepl("^[nN]",ans)){
    #       x <- roiFilter(x,roi_type="negative")
    #     }
    #     cat("More positive ROIs, negative ROIs, or Cancel?")
    #     ans <- readline(prompt="P/N/C []")
    #   }
    #   x <- isPassedROIs(x)
    # }
    return(x)
  }
)

#'@export
#'
hist_fun <- function(x,n=1000,ths,mar=c(3,4,4,2)+.1,brks1,ttl1){
  omar <- par()$mar
  par(mar=mar)

  n.ab <- trim_fun(x,trim_th=1e-2)
  min.i <- which.min(abs(brks1-min(n.ab)))
  max.i <- which.min(abs(brks1-max(n.ab)))
  # stop(list(min.i,max.i))
  uniq.cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"Spectral")))(max.i-min.i+1)
  # stop(length(uniq.cols))
  cols <- c(rep(uniq.cols[1],min.i-1),uniq.cols,rep(rev(uniq.cols)[1],n-max.i+1))
  a <- hist(x,breaks=brks1,main=ttl1,freq=FALSE,xlab="",col=cols,border=NA)

  ## smoothening the trail of histogram
  loessMod <- loess(a$density[seq(n)] ~ brks1[seq(n)], span=0.02)
  smoothened <- predict(loessMod)
  lines(smoothened, x=brks1[seq(n)], col=1,lwd=2)

  # cat("Showing current dna_thres\n")
  abline(v=ths[1],col=4,lty=2,lwd=2)
  abline(v=ths[2],col=2,lty=2,lwd=2)
  par(mar=omar)
  invisible(list(a=a,smoothened=smoothened))
}
