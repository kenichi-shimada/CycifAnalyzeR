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
  function(x,manual=TRUE,n=1000,n1=1000,show.only=FALSE){
    mat <- x@dna
    smp <- x@name

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
        a <- hist(m,breaks=brks,main=ttl,freq=FALSE,xlab="",col="grey60",border=NA)

        ## smoothening the trail of histogram
        loessMod <- loess(a$density[seq(n1)] ~ brks[seq(n1)], span=0.02)
        smoothed <- predict(loessMod)
        lines(smoothed, x=brks[seq(n1)], col=1,lwd=2)

        if(nrow(dna.thres)>0 & show.only){
          cat("Showing current dna_thres\n")
          abline(v=dna.ths1[i],col=4,lty=2)
          abline(v=dna.ths2[i],col=2,lty=2)
        }

        ##
        if(!show.only){
            cat(paste0("Filtering ",smp,"...\n"))

            cat("Showing current dna_thres\n")
            abline(v=dna.ths1[i],col=4,lty=2)
            abline(v=dna.ths2[i],col=2,lty=2)

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

              ## show current dna.ths1 and dna.ths2
              cat("Showing updated threholds for",channel,".\n")
              abline(v=dna.ths1[i],col=4)
              abline(v=dna.ths2[i],col=2)
            }


            if(manual){
              cat("Do you want to modify the 'dna.ths'?\n")
              ans <- "init"
              lo <- dna.ths1[i]
              hi <- dna.ths2[i]
              while(ans!="0"){
                ans <- readline(prompt="0:no, 1:low, 2:high, 3:both [0]")
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
              }

              dna.ths1[i] <- lo
              dna.ths2[i] <- hi
              #              inds.v[[i]] <- ind.v <- ind.v[ind.v %in% idx]
              #              inds.p[[i]] <- ind.p <- ind.p[ind.p %in% idx]
            }
          }

        # ## all DNA channels
        # inds.v <- inds.p <- list()
        # par(mar=c(3,3,1,1))
        # par(bg="white",fg="black")
        #
        # xmax <- max(mat)
        # brks <- seq(0,xmax,length=(n+1))
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
