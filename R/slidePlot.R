setGeneric("slidePlot", function(x,...) standardGeneric("slidePlot"))
setMethod("slidePlot", "Cycif",
	function(x,pch=".",type=c("cell_type","smpl","exp"),col,uniq.col,...){
		xy <- xys(x)
		xy$Y_centroid <- max(xy$Y) - xy$Y
})


#' @export
setGeneric("plotAvailCellOnSlide",function(x,...) standardGeneric("plotAvailCellOnSlide"))
setMethod("plotAvailCellOnSlide", "Cycif",
	function(x,upside.down=TRUE,ncycle,mfrow=c(3,3),mar=c(0,0,4,0),legend=TRUE,
		uniq.cols=c(lost="grey80",dropped="blue",available="black",bunched="red"),legend.cex=2,...){
	stopifnot(nrow(x@used_cells)>0)

	u <- x@used_cells # not to be replaced with cumUsedCells
	nchannels <- ncol(u)

	nc.ratio <- round(statUsedCells(x)*100,1)

	xy <- xys(x)
	if(upside.down){
		xy$Y_centroid <- max(xy$Y) - xy$Y
	}

	par(mfrow=mfrow)
	v <- sapply(seq(nchannels),function(i){
		if(i==1){
			id <- rep(2,nrow(u))
		}else{
			id <- u[,i] + 1
			is.avail <- rowSums(u[,seq(i-1),drop=F]==1)==i-1
			id[!is.avail] <- 0
		}

		rs <- id + 1

		main.txt <- paste0("cycle ",i-1," (",nc.ratio[i],"% available)")
		par(mar=mar)
		plot(xy,pch='.',col=uniq.cols[rs],main=main.txt,axes=F,asp=1,
			xlab="",ylab="",...)
		box()
	})
	if(legend){
		plot.new()
		par(mar=c(1,1,1,1))
		# box()
		legend("topleft",c("available","dropped","bunched","lost previously"),
			fill=uniq.cols[c("available","dropped","bunched","lost")],
			cex=legend.cex)
	}
})
