#' @export
setGeneric("dimPlot", function(x,...) standardGeneric("dimPlot"))

#' @export
setMethod("dimPlot", "CycifStack",
          function(x,pch=".",type=c("cell_type","smpl","exp"),celltype=cts,ab,uniq.cols,main,...){
            require(RColorBrewer)
            ld <- x@ld_coords
            rn <- rownames(ld)
            if(missing(uniq.cols)){
              uniq.cols <- brewer.pal(11,"Spectral")
            }
            if(type=="cell_type"){
              dc <- do.call(c,cts)
              uniq.cols <- c(RColorBrewer::brewer.pal(11,"Spectral"),"grey80","black")
              # smpls <- factor(sub("\\..+","",rn))
              # cols <- uniq.cols[dc[rn]]
              cols[smpls != "5633"] <- NA
              ttl <- "Cell type"
            }else if(type=="smpl"){
              smpls <- factor(sub("\\..+","",rn))
              # stop(table(smpls))
              cols <- as.numeric(smpls)
              ttl <- "Samples"
            }else if(type=="exp"){
              stopifnot(!missing(ab))
              n <- x@normalized
              n <- n[rn,ab]

              rng <- range(n,na.rm=T)
              nn <- (n - rng[1])/diff(rng)

              nc <- 100
              uniq.cols <- rev(RColorBrewer::brewer.pal(11,"Spectral"))
              cols <- colorRampPalette(uniq.cols)(nc+1)[round(nn*nc)+1]
              ttl <- paste0(ab, "(expression)")
            }
            if(missing(main)){
              main <- ttl
            }
          }
)

#' @export
setMethod("dimPlot", "CycifStack",
          function(x,pch=".",type=c("cell_type","smpl","exp"),celltype=cts,ab,uniq.cols,main,...){
            require(RColorBrewer)
            ld <- x@ld_coords
            rn <- rownames(ld)
            if(missing(uniq.cols)){
              uniq.cols <- brewer.pal(11,"Spectral")
            }
            if(type=="cell_type"){
              dc <- do.call(c,cts)
              uniq.cols <- c(RColorBrewer::brewer.pal(11,"Spectral"),"grey80","black")
              cols <- uniq.cols[dc[rn]]
              ttl <- "Cell type"
            }else if(type=="smpl"){
              smpls <- factor(sub("\\..+","",rn))
              # stop(table(smpls))
              cols <- as.numeric(smpls)+1
              ttl <- "Samples"
              pts <- levels(smpls)
            }else if(type=="exp"){
              stopifnot(!missing(ab))
              n <- x@normalized
              n <- n[rn,ab]

              rng <- range(n,na.rm=T)
              nn <- (n - rng[1])/diff(rng)

              nc <- 100
              uniq.cols <- rev(RColorBrewer::brewer.pal(11,"Spectral"))
              cols <- colorRampPalette(uniq.cols)(nc+1)[round(nn*nc)+1]
              ttl <- paste0(ab, "(expression)")
            }
            if(missing(main)){
              main <- ttl
            }

            par(mar=c(4,4,4,10))
            plot(ld,col=cols,pch=".",main=main)
            par(xpd=F)
            legend(par()$usr[2],par()$usr[4],)
            par(xpd=T)
          }
)
