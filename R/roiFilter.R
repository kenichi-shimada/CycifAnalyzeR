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


# pos.rois[rts=="negative"]
# pos.rois <- list(
#   list(x=1:5,y=2:6,roi_type="positive"),
#   list(x=1:5,y=2:6,roi_type="negative"),
#   list(x=1:5,y=2:6,roi_type="positive")
# )

