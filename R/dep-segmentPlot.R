setGeneric("segmentPlot", function(x,...) standardGeneric("segmentPlot"))
setMethod("segmentPlot", "Cycif",
	function(x,extent=NULL,ch="DNA0",cycle=0,pch=20,cex=0.3,add=FALSE,
		action=c("none","draw","crop"),
		params=list(shape="boxROI",center,height=1000,width=1000),clust=FALSE,
		main="",trim=0.01,int.range=NULL,col=NULL,
		...){
	### TBD ###
	action <- action[1]
	d <- x@dim[1,c("x","y")]
	xy <- x@xy_coords
	xy$Y_centroid <- d$y - xy$Y_centroid

	if(!is.null(col) & nrow(xy)!=length(col)){
		stop("nrow(xy) and length(col) should be the same length\n")
	}

	ab.cycle <- (x@abs_list %>% filter(ab == ch))$cycle
	this.cycle <- max(cycle,ab.cycle)

	## is.used -> shown cells
	if(nrow(x@used_cells)>0){
		stopifnot(nrow(x@used_cells)==nrow(x@raw))
		this.dna <- paste0("DNA",this.cycle)
		idx <- match(this.dna,names(x@used_cells))
		is.used <- apply(x@used_cells[seq(idx)],1,all)
	}else{
		is.used <- rep(TRUE,nrow(xy))
	}
	xy <- xy[is.used,]

	##
	if(!is.null(extent) && class(extent)=="Extent"){
		used <- xy$X_Centroid >= extent@xmin & xy$X_Centroid <= extent@xmax &
			xy$Y_Centroid >= extent@ymin & xy$Y_Centroid <= extent@ymax
	}

	if(!is.null(col)){
		cols <- col[is.used]
		is.used.1 <- !is.na(cols)
		cols <- cols[is.used.1]
		xy <- xy[is.used.1,]
	}else if(!clust){
		cols <- scaled_color(x,subset=is.used,ch=ch,trim=trim,int.range=int.range)
	}else{
		cols <- c(x@clusters+1)[is.used]
		cols[cols==1] <- "white"
	}

	if(action=="crop"){
		xmin <- params$center$x - params$width/2
		xmax <- params$center$x + params$width/2
		ymin <- params$center$y - params$width/2
		ymax <- params$center$y + params$width/2
		rng.x <- c(xmin,xmax)
		rng.y <- c(ymin,ymax)
		selected <- xy$X > xmin & xy$X < xmax & xy$Y > ymin & xy$Y < ymax
		xy <- xy[selected,]
		cols <- cols[selected]
	}else{
		rng.x <- c(0,d$x)
		rng.y <- c(0,d$y)
	}

	stopifnot(nrow(xy)==length(cols))

	if(add){
		return(list(x=xy$X,y=xy$Y,col=cols))
	}else{
		plot(rng.x,rng.y,asp=1,axes=F,type="n",main=main,...)
		polygon(rng.x[c(1,1,2,2)],rng.y[c(1,2,2,1)],col="black")
		points(xy,pch=pch,cex=cex,col=cols)

		# if(action %in% c("draw","crop")){
		if(action %in% "draw" & !is.null(params)){
			shape <- draw.box.coords(params)
			lines(shape,col=2,lwd=3)
		}
	}
})
