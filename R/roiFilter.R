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
    within.rois <- within.rois[ncycle,]

    x@within_rois <- within.rois
    return(x)
})

# pos.rois[rts=="negative"]
# pos.rois <- list(
#   list(x=1:5,y=2:6,roi_type="positive"),
#   list(x=1:5,y=2:6,roi_type="negative"),
#   list(x=1:5,y=2:6,roi_type="positive")
# )

