setGeneric("slidePlot", function(x,...) standardGeneric("slidePlot"))
setMethod("slidePlot", "Cycif",
	function(x,pch=".",type=c("smpl","cell_type","exp"),ttl,ab,col,uniq.col,...){
	  stopifnot(!missing(ab))
	  n <- exprs(x,type="normalized")
	  smpl <- names(x)
	  nlev <- 50
	  if(type[1]=="smpl"){
  	  if(!ab %in% abs_list(x)$ab){
  	    stop(ab, " is not available in the sample ", names(x))
  	  }else{
    	  n.ab <- trim_fun(n[[ab]],trim_th=1e-3)
    	  rn <- range(n.ab,na.rm=T)
    	  # stop(sum(is.na(rn)))
        is.na <- is.na(n.ab)
    	  idx <- round((n.ab[!is.na] - rn[1])/diff(rn)*(nlev-1))+1
  	  }
  	 uniq.cols <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(50)

  	 cols <- rep("grey80",length(n.ab))
  	 cols[!is.na] <- uniq.cols[idx]
	  }

  xy <- xys(x)
	xy$Y_centroid <- max(xy$Y) - xy$Y
	if(missing(ttl)){
	  ttl <- paste0(smpl,", ",ab)
	}
	par(mar=c(3,3,3,3))
	plot(xy$X,xy$Y,main=ttl,asp=1,xlab="",ylab="",type="n")
	points(xy$X[is.na],xy$Y[is.na],col="grey80",pch=".")
	points(xy$X[!is.na],xy$Y[!is.na],col=cols[!is.na],pch=".")
})


#' @export
setGeneric("plotAvailCellOnSlide",function(x,...) standardGeneric("plotAvailCellOnSlide"))
setMethod("plotAvailCellOnSlide", "Cycif",
	function(x,upside.down=TRUE,ncycle,mfrow=c(3,3),mar=c(0,0,4,0),legend=TRUE,main=names(x),cex.title=1,
		uniq.cols=c(lost="grey80",dropped="blue",available="black",bunched="red"),legend.cex=2,...){
	stopifnot(nrow(x@used_cells)>0)

	u <- x@used_cells # not to be replaced with cumUsedCells
	nchannels <- ncol(u)

	nc.ratio <- round(statUsedCells(x)*100,1)

	xy <- xys(x)
	if(upside.down){
		xy$Y_centroid <- max(xy$Y) - xy$Y
	}

	par(oma=c(2,2,8,2))
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
	if(!missing(main)){
  	mtext(main, side=3, line=2, cex=cex.title, outer=TRUE)
	}
})
