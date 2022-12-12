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
    n=1000
    n1=1000
    manual=TRUE

    mat <- x@dna
    smpl <- x@name

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
        slidePlot(x,type="dna",ab="DNA3",mar=c(3,3,0,3),ttl="")
        hist_fun(x=m,n=n1,ths=c(dna.ths1[i],dna.ths2[i]),brks1=brks,ttl1=ttl)

        ##
        if(!show.only){
            cat(paste0("Filtering ",smpl,"...\n"))

            ans <- readline(prompt="Do you want to run auto_filter instead of using existing thresholds?[(y)/n]:")
            auto_filter <- !grepl("^[nN]",ans)

            if(auto_filter){
                ## auto_filter thresholding
                th.up <- quantile(m,.9)

                ind.v <- which(diff(sign(diff(smoothed)))>0  & a$mids[c(-1,-n1)] < th.up) # idx of valleys
                ind.p <- which(diff(sign(diff(smoothed)))<0  & a$mids[c(-1,-n1)] < th.up) # idx of peaks
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

                in.rng <- factor((m > dna.ths1[i]) + (m > dna.ths2[i]) + 1,levels=c(1:3))
                cexs <- c(1,20)[(in.rng %in% as.character(c(1,3)))+1]

                ## show current dna.ths1 and dna.ths2
                slidePlot(x,type="filter",cell_type=in.rng,
                          uniq.cols=c("blue","grey80","red"),
                          cex=2,
                          mar=c(3,3,0,3),ttl="")
                hist_fun(x=m,n=n1,ths=c(dna.ths1[i],dna.ths2[i]),brks1=brks,ttl1=ttl)

            }

            if(manual){
                cat("How do you want to modify the 'dna.ths'?\n")
                ans <- "init"
                lo <- dna.ths1[i]
                hi <- dna.ths2[i]
                while(ans!="0"){
                  ans <- readline(prompt="0:none, 1:lower th, 2:higher th, 3:both ths, 4:check fl dist [0-4]")
                  ans <- sub("^(.).*","\\1",ans)
                  if(ans %in% as.character(1:3)){
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

                    in.rng <- factor((m > lo) + (m > hi) + 1,levels=c(1:3))
                    cexs <- c(1,20)[(in.rng %in% as.character(c(1,3)))+1]

                    ## update
                    slidePlot(x,type="filter",cell_type=in.rng,
                              uniq.cols=c("blue","grey80","red"),
                              cex=2,
                              mar=c(3,3,0,3),ttl="")
                    hist_fun(x=m,n=n1,ths=c(lo,hi))
                  }else if(ans=="4"){
                    slidePlot(x,type="dna",ab=channel,mar=c(3,3,0,3),ttl="")
                    hist_fun(x=m,n=n1,ths=c(lo,hi),brks1=brks,ttl1=ttl)
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

    x@dna_thres <- data.frame(low=dna.ths1,high=dna.ths2)
    x@used_cells <- used.cells
    # validObject(x)
    return(x)
  }
)

#'@export
#'
hist_fun <- function(x,n,ths,mar=c(3,4,4,2)+.1,brks1,ttl1,...){
  omar <- par()$mar
  par(mar=mar)
  a <- hist(x,breaks=brks1,main=ttl1,freq=FALSE,xlab="",col="grey60",border=NA)

  ## smoothening the trail of histogram
  loessMod <- loess(a$density[seq(n)] ~ brks1[seq(n)], span=0.02)
  smoothed <- predict(loessMod)
  lines(smoothed, x=brks1[seq(n)], col=1,lwd=2)

  # cat("Showing current dna_thres\n")
  abline(v=ths[1],col=4,lty=2)
  abline(v=ths[2],col=2,lty=2)
  par(mar=omar)
}
