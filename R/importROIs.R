#' @export
setGeneric("importROIs", function(x,...) standardGeneric("importROIs"))
setMethod("importROIs", "Cycif",
  function(x,exported.rois,plot=T){
    if(missing(exported.rois) || nrow(exported.rois) ==0){
      stop("'exported.rois' with one or more rows should be provided.")
    }
    max.y <- max(xys(x)$Y_centroid)
    lst <- list()
    for(k in seq(nrow(exported.rois))){
      tmp <- exported.rois[k,]
      roi.text <- tmp$Text
      roi.cycle <- as.numeric(sub(".+([0-9]+)$","\\1",roi.text))
      if(grepl("pos",roi.text)){
        roi.dir <- "Positive"
        lcol <- 2
      }else if(grepl("neg",roi.text)){
        roi.dir <- "Negative"
        lcol <- 4
      }else{
        return(roi.text)
      }
      roi.type <- tmp$type
      if(roi.type=="Polygon"){
        coords <- do.call(rbind,strsplit(unlist(strsplit(tmp$all_points," ")),","))
        df <- as.data.frame(apply(coords,2,as.numeric))
        names(df) <- c("x","y")
        df$y <- max.y - df$y
        if(plot){
          polygon(df,lty=1,border=lcol)
        }
      }else if(roi.type=="Ellipse"){
        ps <- tmp[c("X","Y","RadiusX","RadiusY")]
        ps$Y <- max.y - ps$Y
        if(plot){
          plotrix::draw.ellipse(x=ps$X,y=ps$Y,a=ps$RadiusX,b=ps$RadiusY,border=lcol)
        }
      }
      roi.obj <- list(
        dir=roi.dir,
        cylce=roi.cycle,
        type=roi.type,
        coords=df
      )
      lst[[k]] <- roi.obj
    }
    return(lst)
  }
)
