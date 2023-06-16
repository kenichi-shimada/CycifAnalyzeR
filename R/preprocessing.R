#_ -------------------------------------------------------

# fun: importROIs Cycif ----

#' Import Polygon coordinates of ROIs from OMERO
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate select
#' @export
setGeneric("importROIs", function(x,...) standardGeneric("importROIs"))

#' @rdname importROIs
#' @export
setMethod("importROIs", "Cycif",
          function(x,exported.rois){
            if(missing(exported.rois) || nrow(exported.rois) == 0){
              stop("'exported.rois' with one or more rows should be provided.")
            }
            # exported.rois <- exported.rois[[smpl]]
            exported.rois <- exported.rois %>%
              dplyr::filter(type %in% c("Rectangle","Polygon","Ellipse")) %>%
              dplyr::filter(grepl("^(pos|neg)([0-9]+)$",Text)) %>%
              dplyr::mutate(Dir = sub("^(pos|neg)([0-9]+)$","\\1",Text)) %>%
              dplyr::mutate(Cycle = as.numeric(sub("^(pos|neg)([0-9]+)$","\\2",Text))) %>%
              dplyr::select(-Text)
            max.y <- max(xys(x)$Y_centroid)

            ## Each ROI should have
            # direction: pos or neg
            # cycle: 1-n
            # roi.type: Polygon
            # Coordinates (polygon), or parameters (X,Y, RadiusX,Radius Y) (ellipse)

            lst.rois <- lapply(seq(nrow(exported.rois)),function(i){
              tmp <- exported.rois[i,]
              if(tmp$Dir=="pos"){
                roi.dir <- "positive"
                # lcol <- 2
              }else if(tmp$Dir=="neg"){
                roi.dir <- "negative"
                # lcol <- 4
              }else{
                return(roi.dir)
              }
              roi.type <- tmp$type
              if(roi.type=="Polygon"){
                coords <- do.call(rbind,strsplit(unlist(strsplit(tmp$all_points," ")),","))
                df <- as.data.frame(apply(coords,2,as.numeric))
                names(df) <- c("x","y")
                df$y <- max.y - df$y
                # if(plot){
                #   polygon(df,lty=1,border=lcol)
                # }
              }else if(roi.type=="Ellipse"){
                df <- tmp[c("X","Y","RadiusX","RadiusY")]
                df$Y <- max.y - df$Y
                res <- 30
                theta = seq(0, 2 * pi, length = res)
                x = df$X - df$RadiusX * cos(theta)
                y = df$Y - df$RadiusY * sin(theta)
                df <- as.data.frame(cbind(x=x,y=y))
                roi.type <- "Polygon"
                # if(plot){
                #   plotrix::draw.ellipse(x=ps$X,y=ps$Y,a=ps$RadiusX,b=ps$RadiusY,border=lcol)
                # }
              }else if(roi.type=="Rectangle"){
                df <- tmp[c("X","Y","Width","Height")]
                df$Y <- max.y - df$Y
                x <- df$X + df$Width/2 * c(-1,1,1,-1,-1)
                y <- df$Y + df$Height/2 * c(-1,-1,1,1,-1)
                df <- as.data.frame(cbind(x=x,y=y))
                roi.type <- "Polygon"
              }
              roi.obj <- list(
                dir=roi.dir,
                cycle=tmp$Cycle,
                roi_type=roi.type,
                coords=df
              )
              return(roi.obj)
            })

            ## match ncycle
            ncycle <- nCycles(x)
            matched.cys <- sapply(lst.rois,function(x)x$cycle) <= ncycle
            # stop(matched.cys)
            # stop(length(lst.rois))
            lst.rois <- lst.rois[matched.cys]

            x@rois <- lst.rois
            x <- setWithinROIs(x)
            return(x)
          }
)

#' @rdname importROIs
#' @export
setMethod("importROIs", "CycifStack",
          function(x,exported.rois){
            smpls1 <- names(x)
            nm.rois <- names(exported.rois)
            if(!any(smpls1 %in% nm.rois)){
              unused <- smpls1[!smpls1 %in% nm.rois]
              stop("ROIs not available for samples: ",paste(unused,collapse=","))
            }
            for(smpl in smpls1){
              cat("Loadin ROI for ",smpl,"...\n")
              x[[smpl]] <- importROIs(x[[smpl]],exported.rois=exported.rois[[smpl]])
            }
            return(x)
          }
)

# fun: setWithinROIs Cycif ----

#' @importFrom sp point.in.polygon
#' @export
setGeneric("setWithinROIs", function(x,...) standardGeneric("setWithinROIs"))
setMethod("setWithinROIs", "Cycif",
          function(x){
            rois <- x@rois
            ncycles <- sapply(rois,function(x)x$cycle)
            coords <- xys(x)
            nc <- nCells(x)
            ncycle <- nCycles(x)

            within.rois <- sapply(seq(ncycle),function(i){
              is.used.rois <- ncycles <= i
              if(sum(is.used.rois)==0){
                within.rois <- rep(TRUE,nc)
                return(within.rois)
              }
              rois1 <- rois[is.used.rois]
              rts <- sapply(rois1,function(r)r$dir)
              if(any(rts=="positive")){
                pos.rois <- rois1[rts=="positive"]
                passed.pos.rois <- as.matrix(sapply(pos.rois,function(pr){
                  xys2 <- pr$coords
                  within.rois <- sp::point.in.polygon(coords$X,max(coords$Y)-coords$Y,xys2$x,xys2$y)==1
                }))
              }else{
                passed.pos.rois <- as.matrix(rep(TRUE,nCells(x)))
              }
              if(any(rts=="negative")){
                neg.rois <- rois1[rts=="negative"]
                passed.neg.rois <- as.matrix(sapply(neg.rois,function(nr){
                  xys2 <- nr$coords
                  within.rois <- sp::point.in.polygon(coords$X,max(coords$Y)-coords$Y,xys2$x,xys2$y)==0
                }))
              }else{
                passed.neg.rois <- as.matrix(rep(TRUE,nCells(x)))
              }
              wr <- apply(passed.pos.rois,1,any) & apply(passed.neg.rois,1,all)
              return(wr)
            })
            within.rois <- within.rois[,ncycle]

            x@within_rois <- within.rois
            return(x)
          }
)

# fun: roiFilter Cycif ----

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
setGeneric("roiFilter", function(x,...) standardGeneric("roiFilter"))
setMethod("roiFilter", "Cycif",
          function(x,rois){
            mat <- x@raw
            smpl <- x@name
            # n <- n1 <- 1000

            pos.rois <- x@rois
            ## choose positive ROI
            if(roi_type=="positive"){
              cat("Set positive ROIs.\n")
              ans <- "Y"
              while(!grepl("^[nN]",ans)){
                ns <- as.integer(readline(prompt="How many points?"))
                cat(paste0("Select ",ns," points to set a polygon\n"))
                xys1 <- locator(ns)
                lines(xys1$x[c(seq(xys1$x),1)],xys1$y[c(seq(xys1$x),1)],col=2,lty=2,lwd=2)
                check <- readline(prompt="satisfied with the ROI? (Y/N) [Y]")
                if(grepl("^[nN]",check)){
                  next
                }
                xys1$roi_type="positive"
                pos.rois <- c(pos.rois,list(xys1))
                cat("Do you want to set more positive ROIs?")
                ans <- readline(prompt="(Y/N) [Y]")
                x@rois <- pos.rois
              }
              return(x)
            }else if(roi_type=="negative"){
              cat("Set negative ROIs.\n")
              ans <- "Y"
              while(!grepl("^[nN]",ans)){
                ns <- as.integer(readline(prompt="How many points?"))
                cat(paste0("Select ",ns," points to set a polygon\n"))
                xys1 <- locator(ns)
                lines(xys1$x[c(seq(xys1$x),1)],xys1$y[c(seq(xys1$x),1)],col=4,lty=2,lwd=2)
                check <- readline(prompt="satisfied with the ROI? (Y/N) [Y]")
                if(grepl("^[nN]",check)){
                  next
                }
                xys1$roi_type="negative"
                pos.rois <- c(pos.rois,list(xys1))
                cat("Do you want to set more negative ROIs?")
                ans <- readline(prompt="(Y/N) [Y]")

                x@rois <- pos.rois
              }
            }
            return(x)
          }
)

# fun: dnaFilter Cycif ----

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
              slidePlot(x,plot_type="dna",ab=channel,mar=c(3,3,0,3),ttl="",use_rois=FALSE)
              hst <- hist_fun(x=m,ths=c(dna.ths1[i],dna.ths2[i]),brks1=brks,ttl1=ttl)

              smoothened <- hst$smoothened
              a <- hst$a

              ##
              if(!show.only){
                if(i==1){
                  cat(paste0("Filtering ",smpl,", cycle ",i,"...\n"))
                }
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
                        # stop("done")
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
                      # stop(paste(table(cexs),collapse=","))

                      ## update
                      slidePlot(x,
                                plot_type="custom",
                                custom_labs=in.rng,
                                uniq.cols=c("grey80",4,2),
                                cex=cexs,
                                mar=c(3,3,0,3),
                                ttl="",
                                use_rois=FALSE)
                      # stop("after slideplot")
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

            cat("Cell retention through all cycles are plotted.\n")
            l <- layout(matrix(c(1,2),nrow=2),heights=c(2,3))
            plotUsedCellRatio(x)
            slidePlot(x,plot_type="custom",
                      custom_labs=ret,
                      uniq.cols=uniq.cols,
                      cex=2,
                      ncells=1e4,
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

#' @rdname dnaFilter
#' @export
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


#_ -------------------------------------------------------

# fun: exprs Cycif, CycifStack ----

#' Show a raw or normalized protein expression matrix
#'
#' @param x A Cycif or CycifStack object.

#' @param type type. If "raw", a raw expression matrix is shown. If "normalized",
#'   a normalized expression matrix will be shown. The method of normalization is
#'   specified by "log" or "LogTh" in the normalize function provided previously.
#' @param silent logical. If FALSE, the method of normalization was shown in
#'  the error prompt.
#'
#' @export
setGeneric("exprs", function(x,...) standardGeneric("exprs"))

#' @export
setMethod("exprs", "Cycif", function(x,type=c("raw","log_normalized","logTh_normalized")){
  if(missing(type)){
    # cat("Error: specify type: \"raw\" or \"normalized\"\n")
    stop("'type'argument should be specified\n")
  }

  if(type=="log_normalized" && !("log_normalized" %in% slotNames(x))){
    stop("Error: slot @log_normalized not found.\n")
  }

  if(type=="logTh_normalized" && !("logTh_normalized" %in% slotNames(x))){
    stop("Error: slot @logTh_normalized not found.\n")
  }

  if(type == "raw"){
    return(x@raw)
  }else if(type == "log_normalized"){
    return(x@log_normalized)
  }else if(type == "logTh_normalized"){
    return(x@logTh_normalized)
  }
})

#' @export
setMethod("exprs", "CycifStack",function(x,type=c("raw","log_normalized","logTh_normalized")){
  if(missing(type)){
    stop("'type' should be specified\n")
  }

  all.abs <- abs_list(x)$ab
  exps <- do.call(rbind,cyApply(x,function(cy){
    tmp <- exprs(cy,type=type)
    if(any(!all.abs %in% names(tmp))){
      unused <- all.abs[!all.abs %in% names(tmp)]
      tmp1 <- data.frame(array(NA,dim=c(nrow(tmp),length(unused))))
      names(tmp1) <- unused
      tmp <- cbind(tmp,tmp1)[all.abs]
    }
    return(tmp)
  }))
  return(exps)
})

# fun: normalize Cycif ----

#' Show a raw or normalized protein expression matrix
#'
#' @param r A numeric scholar indicating a raw expression value.
#' @param x A Cycif or CycifStack object.
#' @param method Either "log" or "logTh". "log" indicates that
#'   the raw values are transformed using log1p() function. "logTh" method
#'   further trim the outliers outside outside \[trim,1-trim\] quantiles,
#'   where trm is specified by trim parameter. Outlier-trimmed, log1p-transformed
#'   raw values are further transformed such that the value indicated by the threshold,
#'   th (in the case of transform() function) or threshold(x) (in the case of
#'   normalize() function) will be set to p_thres. Values below or above the threshold
#'   are linearly transformed to values between 0 and p_thres, or values between p_thres and 1,
#'   respectively.
#' @param th A numeric scholar. A user-provided threshold for a raw expression value
#'   for each protein for each sample.
#' @export
#' @rdname normalize
setGeneric("normalize", function(x,...) standardGeneric("normalize"))

#' @rdname normalize
#' @export
setMethod("normalize", "Cycif",
          function(x,method=c("log","logTh","invlog"),trim=1e-3,p_thres=0.5){
            ## default method is logTh
            if(missing(method)){
              stop("normalize() should specify method: log, logTh, or invlog")
            }

            used.abs <- as.character(abs_list(x)$ab)

            smpl <- x@name
            raw <- x@raw
            is.used <- cumUsedCells(x)

            ## treatment is different between methods
            if(method=="log"){
              norm <- as.data.frame(
                sapply(used.abs,function(ab){
                  cycle <- abs_list(x)$cycle[abs_list(x)$ab==ab]
                  is.used.1 <- is.used[,cycle]

                  r <- raw[[ab]]
                  n <- rep(NA,length(r))

                  n[is.used.1] <- transform(r[is.used.1],method="log",trim=trim)
                  return(n)
                })
              )
              x@log_normalized <- norm
            }else if(method=="logTh"){
              gt <- paste0("gates_",names(x))
              if(!gt %in% names(abs_list(x))){
                stop('Set gates first with setGates()')
              }else{
                gates <- abs_list(x)[c("ab",gt)]
              }

              log_thres <- gates[[gt]]
              names(log_thres) <- gates$ab

              thres <- expm1(log_thres)
              used.abs <- names(which(!is.na(thres)))

              norm <- as.data.frame(
                sapply(used.abs,function(ab){
                  # for(ab in used.abs){
                  n <- rep(NA,nrow(raw))
                  if(ab %in% used.abs){
                    cycle <- abs_list(x)$cycle[abs_list(x)$ab==ab]
                    is.used.1 <- is.used[,cycle]
                    r <- raw[[ab]]
                    th <- thres[ab]
                    n[is.used.1] <- transform(r[is.used.1],method="logTh",th=th,trim=trim,p_thres=p_thres)
                  }
                  # }
                  return(n)
                })
              )
              x@logTh_normalized <- norm
            }else if(method=="invlog"){
              norm <- x@log_normalized
              raw <- expm1(norm)
              x@raw <- raw
            }

            return(x)
          }
)

# fun: transform numeric ----

# at some point I should switch to nls() to apply sigmoidal curve
#' @rdname normalize
#' @export
transform <- function(r,method=c("log","logTh","Th","invlog"),th,p_thres=0.5,trim=1e-3){
  if(missing(method)){
    stop("normalize() should specify method: log, logTh, or invlog")
  }

  # r - raw intensity value (in quantification/*.csv)

  ## log-transformation (to be replaced with log with other bases)
  if(method=="Th"){
    r1 <- lr <- r
    lth <- th
  }else if(method=="logTh"){ ## when method contains log-transformation
    r1 <- lr <- log1p(r)
    lth <- log1p(th)
  }else if(method=="log"){
    r1 <- lr <- log1p(r)
  }else if(method=="invlog"){
    r1 <- expm1(r)
  }

  ## method: either 'log'
  if(method=="log"){
    qt <- quantile(lr,c(trim,1-trim))
    r1[r1 > qt[2]] <- qt[2]
    r1[r1 < qt[1]] <- qt[1]
  }

  ## method: either 'Th' or 'logTh'
  if(method %in% c("logTh","Th")){
    if(missing(th)){
      stop("'th' should be specified when 'method=\"logTh\"' or 'method=\"Th\"'")
    } # th is an essential input from user

    ## trimming outliers
    qt <- quantile(lr,c(trim,1-trim))

    ## if lth is outside [trim, 1-trim], move lth to the bound
    if(lth < qt[1]){
      qt[1] <- lth
    }
    if(lth > qt[2]){
      qt[2] <- lth
    }

    if(any(lr > lth)){
      r1[lr > lth] <- (lr[lr > lth] - lth)/((qt[2] - lth)/(1-p_thres)) + p_thres
    }
    if(any(lr < lth)){
      r1[lr <= lth] <- (lr[lr <= lth] - qt[1])/((lth - qt[1])/p_thres)
    }
    r1[r1 < 0] <- 0
    r1[r1 > 1] <- 1
  }

  return(r1)
}

#_ -------------------------------------------------------

