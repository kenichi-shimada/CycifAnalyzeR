#' @export
setGeneric("importROIs", function(x,...) standardGeneric("importROIs"))
setMethod("importROIs", "Cycif",
  function(x,exported.rois){
    if(missing(exported.rois) || nrow(exported.rois) == 0){
      stop("'exported.rois' with one or more rows should be provided.")
    }
    # exported.rois <- rois.out[[smpl]]
    exported.rois <- exported.rois %>%
      filter(type %in% c("Rectangle","Polygon","Ellipse")) %>%
      filter(grepl("^(pos|neg)([0-9]+)$",Text)) %>%
      mutate(Dir = sub("^(pos|neg)([0-9]+)$","\\1",Text)) %>%
      mutate(Cycle = as.numeric(sub("^(pos|neg)([0-9]+)$","\\2",Text))) %>%
      select(-Text)
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

    x@rois <- lst.rois
    x <- setWithinROIs(x)
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
    within.rois <- within.rois[,ncycle]

    x@within_rois <- within.rois
    return(x)
  }
)
