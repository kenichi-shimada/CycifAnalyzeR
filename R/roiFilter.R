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
  function(x,roi_type=c("positive","negative")){
    if(missing(roi_type)){
      stop("'roi_type' should be specified")
    }
    mat <- x@raw
    smpl <- x@name
    n <- n1 <- 1000

    used.cells <- x@used_cells
    if(nrow(used.cells)==0){
      stop("Run dnaFilter() first to set 'used_cells' slots")
    }

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
    o <- order(ret,decreasing=T)

    ## final after dnaFitlter
    nc <- length(unique(ret))

    uniq.cols <- colorRampPalette(brewer.pal(9,"YlGnBu"))(nc+2)[-(nc+(0:1))]
    cat("Cell retention through each cycle:\n")
    l <- layout(matrix(c(1,2),nrow=2),heights=c(2,2))
    plotUsedCellRatio(x)
    slidePlot(x,plot_type="filter",within_filter_rng=ret,
              uniq.cols=uniq.cols,
              cell.order=o,
              cex=2,ncells=1e4,
              mar=c(3,3,0,3),ttl="")

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

#' @export
setGeneric("isPassedROIs", function(x,...) standardGeneric("isPassedROIs"))
setMethod("isPassedROIs", "Cycif",
  function(x){
    rois <- x@rois
    coords <- xys(x)
    rts <- sapply(rois,function(r)r$roi_type)
    if(any(rts=="positive")){
      pos.rois <- rois[rts=="positive"]
      passed.pos.rois <- as.matrix(sapply(pos.rois,function(xys2){
        within.rois <- sp::point.in.polygon(coords$X,max(coords$Y)-coords$Y,xys2$x,xys2$y)==1
      }))
    }else{
      passed.pos.rois <- as.matrix(rep(TRUE,nCells(x)))
    }
    if(any(rts=="negative")){
      neg.rois <- rois[rts=="negative"]
      passed.neg.rois <- as.matrix(sapply(neg.rois,function(xys2){
        within.rois <- sp::point.in.polygon(coords$X,max(coords$Y)-coords$Y,xys2$x,xys2$y)==0
      }))
    }else{
      passed.neg.rois <- as.matrix(rep(TRUE,nCells(x)))
    }

    within.rois <- apply(passed.pos.rois,1,any) & apply(passed.neg.rois,1,all)
    x@within_rois <- within.rois
    return(x)
})

# pos.rois[rts=="negative"]
# pos.rois <- list(
#   list(x=1:5,y=2:6,roi_type="positive"),
#   list(x=1:5,y=2:6,roi_type="negative"),
#   list(x=1:5,y=2:6,roi_type="positive")
# )

