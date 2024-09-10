#_ -------------------------------------------------------

# fun: importROIs Cycif ----

#' @title Import Polygon coordinates of ROIs from OMERO
#'
#' @details This function imports Polygon, Ellipse, or Rectangle ROIs from OMERO and adds them to a Cycif or CycifStack objects.
#'
#' @param x A Cycif or CycifStack object.
#' @param exported.rois A data frame containing the exported ROIs with columns:
#'   - `type`: Type of ROI ('Rectangle', 'Polygon', or 'Ellipse').
#'   - `Text`: Text describing the ROI, which includes direction and cycle information (e.g., 'pos5', 'neg7')
#'
#' @details
#' The `exported.rois` data frame should contain ROIs with valid direction ('pos' or 'neg') and cycle information.
#' Supported ROI types are 'Rectangle', 'Polygon', and 'Ellipse'. The function processes these ROIs and adds
#' them to the Cycif object. Each ROI in `lst.rois` consists of the following components:
#'   - `dir`: Direction of the ROI ('positive' or 'negative').
#'   - `cycle`: Cycle information for the ROI (numeric).
#'   - `roi_type`: Type of the ROI ('Polygon').
#'   - `coords`: If `roi_type` is 'Polygon', `coords` is a data frame containing the coordinates of the ROI (x, y).
#'     For other `roi_type` values, `coords` may contain different information based on the shape of the ROI.
#'
#' @return A Cycif object with imported ROIs.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate select
#'
#' @seealso
#' \code{\link{setWithinROIs}} \code{\link{roiFilter}}
#'
#' @rdname importROIs
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
                coords <- data.frame(do.call(rbind,strsplit(unlist(strsplit(tmp$all_points," ")),",")))
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
                x <- df$X + df$Width * c(0,1,1,0,0)
                y <- df$Y + df$Height * c(0,0,1,1,0)
                df$Y <- max.y - df$Y
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

#' @title Set within_roi slot in a Cycif object
#'
#' @description Compute if each cell is within any positive ROIs and outside any negative ROIs. The value is TRUE if that is the case, FALSE otherwise.
#' This should be used after importROIs functions are called so ROIs are set in each Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return A Cycif object with cells marked as within or outside ROIs in within_rois slot
#'
#' @details This function determines which cells in a Cycif object are within the defined ROIs. It uses information about the direction (positive or negative), cycle, and coordinates of the ROIs to classify cells as either within or outside ROIs. Cells are marked as within ROIs if they meet the following conditions:
#' - They are within the polygonal boundaries of positive ROIs.
#' - They are outside the polygonal boundaries of negative ROIs.
#
#' @importFrom sp point.in.polygon
#' @export
setGeneric("setWithinROIs", function(x,...) standardGeneric("setWithinROIs"))

#' @rdname setWithinROIs
#' @export
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
            within.rois <- rowSums(within.rois)==ncycle

            x@within_rois <- within.rois
            return(x)
          }
)

# fun: roiFilter Cycif ----

#' @title Set ROIs in a Cycif object interactively by clicking on the plot (discontinued)
#'
#' @description This function allows interactive setting of regions of interest (ROIs) in a Cycif object. ROIs can be defined as either positive or negative based on their intended use.
#'
#' @param x A Cycif object.
#' @param roi_type A character string specifying the type of ROIs to set. It can be either "positive" or "negative."
#'
#' @return A Cycif object with the defined ROIs.
#'
#' @importFrom graphics locator
#'
#' @export
#' @details This function provides an interactive way to define ROIs within a Cycif object. ROIs can be classified as either "positive" or "negative" based on their intended use. Users can set these ROIs by selecting points on the plot. For positive ROIs, cells within the polygonal boundaries are considered within the ROI. For negative ROIs, cells outside the polygonal boundaries are considered within the ROI. Users can set multiple ROIs of each type by following the prompts.
#'
#' @seealso \code{\link{importROIs}}, \code{\link{setWithinROIs}}
#'
setGeneric("roiFilter", function(x,...) standardGeneric("roiFilter"))

#' @export
#' @rdname roiFilter
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

#_ -------------------------------------------------------
# fun: dnaFilter Cycif ----

#' @title Identify available cells based on cell DNA contents on each cycle in a Cycif dataset interactively
#'
#' @description This function determines whether cells are available through multiple cycles based on the intensity of DNA stains in a Cycif object.
#' It allows interactive setting of thresholds and visualization of the filtering process.
#'
#' @param x A Cycif object with DNA stains.
#' @param ncells Maximum number of cells to analyze. Default is 1e5.
#' @param dna.thres A list containing two data frames with low and high thresholds for DNA channels, or NULL to use default thresholds.
#' @param show.only Logical, indicating whether to perform DNA filtering without interactive adjustments (TRUE) or with interactive adjustments (FALSE).
#'
#' - When set to TRUE, the function performs DNA filtering without user interaction and applies the specified thresholds. Use this option if you want to apply pre-defined DNA thresholds without further adjustments.
#'
#' - When set to FALSE, the function interactively adjusts the DNA thresholds by plotting DNA histograms and allows you to modify the thresholds based on the distribution of DNA content. You can manually adjust the thresholds to retain cells of interest while discarding outliers. Use this option for a more hands-on approach to DNA filtering.
#'
#' The default value is FALSE.
#'
#' @return A Cycif object with DNA filtering applied.
#'
#' @importFrom graphics locator par
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#' @details This function performs DNA filtering for a Cycif object,
#' allowing you to specify the maximum number of cells to analyze (`ncells`) and the thresholds for DNA channels (`dna.thres`).
#' You can also choose whether to interactively adjust the thresholds or perform DNA filtering without interaction using the `show.only` argument.
#'
#' If `show.only` is set to `FALSE`, the function interactively adjusts the thresholds by plotting DNA histograms and allowing you to modify the thresholds based on the distribution of DNA content. You can manually adjust the thresholds to retain cells of interest while discarding outliers.
#'
#' If `show.only` is set to `TRUE`, the function performs DNA filtering without user interaction and applies the specified thresholds.
#'
#' The `dna.thres` argument allows you to specify custom low and high thresholds for DNA channels. If `dna.thres` is NULL, the function uses default thresholds based on the data.
#'
setGeneric("dnaFilter", function(x,...) standardGeneric("dnaFilter"))

#' @export
#' @rdname dnaFilter
setMethod("dnaFilter", "Cycif",
          function(x,ncells=1e5,dna.thres,show.only=FALSE){
            set.seed(1)
            mat <- x@dna

            nr <- nrow(mat)
            ncells <- min(ncells,nr)
            idx <- sample(nrow(mat),ncells)
            mat1 <- mat[idx,]

            x@dna <- mat1
            x@xy_coords <- x@xy_coords[idx,]

            smpl <- x@name
            n <- n1 <- 1000

            dna.list <- names(mat1)
            mat1 <- cbind(log1p(mat1[[1]]),as.data.frame(lapply(mat1,function(x)log1p(x/mat1[[1]])))[-1])
            names(mat1) <- dna.list

            dna.thres <- x@dna_thres

            if(nrow(dna.thres)>0){
              dna.ths1 <- x@dna_thres$low
              dna.ths2 <- x@dna_thres$high
              names(dna.ths1) <- names(dna.ths2) <- dna.list
            }else{
              dna.ths1 <- sapply(mat1,min)
              dna.ths2 <- sapply(mat1,max)
              names(dna.ths1) <- names(dna.ths2) <- dna.list
            }

              if(missing(dna.thres)){
              for(i in seq(length(dna.list))){
                channel <- dna.list[i]
                m <- mat1[[channel]]

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

                  ## manually adjust the thresholds
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
                      slidePlot(x,#ncells=ncells,
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
            }else{
              dna.ths1 <- dna.thres$low
              dna.ths2 <- dna.thres$high
              names(dna.ths1) <- names(dna.ths2) <- names(mat)
            }

            ##
            mat2 <- cbind(log1p(mat[[1]]),as.data.frame(lapply(mat,function(x)log1p(x/mat[[1]])))[-1])
            names(mat2) <- names(mat)
            used.cells <- sapply(names(mat2),function(channel){
              ind <- (mat2[[channel]] >= dna.ths1[channel]) + (mat2[[channel]] > dna.ths2[channel])
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

            if(0){
              cat("Cell retention through all cycles are plotted.\n")
              l <- layout(matrix(c(1,2),nrow=2),heights=c(2,3))
              plotUsedCellRatio(x)
              slidePlot(x,plot_type="custom",
                        custom_labs=ret,
                        uniq.cols=uniq.cols,
                        cex=2,
                        ncells=1e4,
                        mar=c(3,3,0,3),ttl="")
            }
            return(x)
          }
)

# not exported
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

#' @param type type. If "raw", a raw expression matrix is shown. If "log" or "logTh",
#'   a log or logTh normalized expression matrix will be shown.
#' @param silent logical. If FALSE, the method of normalization was shown in
#'  the error prompt.
#'
#' @seealso \code{\link{normalize}} \code{\link{transform}}
#' @export
setGeneric("exprs", function(x,...) standardGeneric("exprs"))

#' @export
#' @rdname exprs
setMethod("exprs", "Cycif", function(x,type=c("raw","log","logTh")){
  if(missing(type)){
    stop("'type'argument should be specified\n")
  }

  if(type=="log" && nrow(x@log_normalized)==0){ #    !("log_normalized" %in% slotNames(x)))
    stop("Error: log_transformed data not found. normalize() first.\n")
  }

  if(type=="logTh" && nrow(x@logTh_normalized)==0){#!("logTh_normalized" %in% slotNames(x))){
    stop("Error: logTh transformed data not found.\n")
  }

  if(type == "raw"){
    return(x@raw)
  }else if(type == "log"){
    return(x@log_normalized)
  }else if(type == "logTh"){
    return(x@logTh_normalized)
  }else{
    stop("'type' argument should be one of the three: raw, log, logTh")
  }
})

#' @export
#' @rdname exprs
setMethod("exprs", "CycifStack",function(x,type=c("raw","log","logTh")){
  if(missing(type)){
    stop("'type' should be specified\n")
  }

  all.abs <- abs_list(x)$ab
  list.exps <- cyApply(x,function(cy){
    tmp <- exprs(cy,type=type)
    if(any(!all.abs %in% names(tmp))){
      unused <- all.abs[!all.abs %in% names(tmp)]
      tmp1 <- data.frame(array(NA,dim=c(nrow(tmp),length(unused))))
      names(tmp1) <- unused
      tmp <- cbind(tmp,tmp1)[all.abs]
    }
    return(tmp)
  })
  exps <- data.frame(data.table::rbindlist(list.exps))
  return(exps)
})

# fun: normalize Cycif ----

#' Normalize a Cycif or CycifStack object
#'
#' This function allows you to normalize a Cycif or CycifStack object, either using a simple logarithmic transformation or a more advanced transformation method with optional thresholding and trimming.
#'
#' @param x A Cycif or CycifStack object.
#' @param method Either "log" or "logTh". "log" indicates that the raw values are transformed using log1p() function. "logTh" method further trims the outliers outside the \[trim, 1-trim\] quantiles, where trim is specified by the trim parameter. Outlier-trimmed, log1p-transformed raw values are further transformed such that the value indicated by the threshold, th (in the case of transform() function) or threshold(x) (in the case of normalize() function) will be set to p_thres. Values below or above the threshold are linearly transformed to values between 0 and p_thres, or values between p_thres and 1, respectively.
#' @param trim A numeric scalar specifying the trim value for outlier trimming.
#' @param p_thres A numeric scalar specifying the threshold value for linear transformation.
#' @param th A numeric scholar. A user-provided threshold for a raw expression value for each protein for each sample.
#'
#' @details
#' The `normalize` function is used to preprocess protein expression data in Cycif or CycifStack objects. It provides two main normalization methods: "log" for simple logarithmic transformation and "logTh" for more advanced normalization with optional thresholding and trimming.
#'
#' For "log" method, the function applies a logarithmic transformation using log1p(), followed by optional outlier trimming using the specified trim parameter.
#'
#' For "logTh" method, in addition to logarithmic transformation and trimming, the function allows users to set a threshold (th) for each protein. Values above or below the threshold are linearly transformed to fit within the specified range, defined by p_thres.
#'
#' @seealso \code{\link{transform}}
#'
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
            is.used <- x@within_rois
            if(length(is.used)==0){
              stop("used_cells slot is empty. Set importROIs() first.")
            }
            # is.used <- cumUsedCells(x)

            ## treatment is different between methods
            if(method=="log"){
              norm <- as.data.frame(
                sapply(used.abs,function(ab){
                  cycle <- abs_list(x)$cycle[abs_list(x)$ab==ab]
                  # is.used.1 <- is.used[,cycle]
                  is.used.1 <- is.used

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
                    # is.used.1 <- is.used[,cycle]
                    is.used.1 <- is.used

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

#' Transform Numeric Data
#'
#' This function applies various transformations to numeric data based on the specified method.
#'
#' @param r A numeric vector or array.
#' @param method The transformation method to apply. Choose from "log", "logTh", "Th", or "invlog".
#'   - "log": Log-transforms the data using log1p().
#'   - "logTh": Log-transforms the data and trims outliers.
#'   - "Th": Trims outliers.
#'   - "invlog": Reverts the data back to its original scale.
#' @param th The threshold value for transformation (only for "Th" and "logTh" methods).
#' @param p_thres The threshold value for data scaling (only for "logTh" and "Th" methods).
#' @param trim The quantile-based trimming level for outlier removal (default is 1e-3).
#'
#' @return Transformed numeric data.
#'
#' @seealso \code{\link{normalize}}
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

# not exported
trim_fun <- function(x,trim_th = 1e-3){
  qts <- stats::quantile(x,c(trim_th,1-trim_th),na.rm=T)
  x[x < qts[1]] <- qts[1]
  x[x > qts[2]] <- qts[2]
  return(x)
}
#_ -------------------------------------------------------

