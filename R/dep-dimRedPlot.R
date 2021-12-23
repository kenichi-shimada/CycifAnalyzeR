
setGeneric("dimRedPlot", function(x,...) standardGeneric("dimRedPlot"))
setMethod("dimRedPlot", "Cycif",
	function(x,ch="DNA0",cycle=0,pch=20,clust=TRUE,label=TRUE,trim=0.01,...){

	d <- x@dim[1,c("x","y")]
	xy <- x@ld_coords

	rng.x <- range(xy$x)
	rng.y <- range(xy$y)

	ggplot(xy,aes=c(x=x,y=y)) +
	geom_point(aes())
	plot(rng.x,rng.y,axes=F,type="n",xlab="",ylab="",...)
	box()
	# mtext("UMAP 1",1,1)
	# mtext("UMAP 1",2,1)

	if(!clust){
		is.used <- !is.na(x@clusters)
		cols <- scaled_color(x,subset=is.used,ch=ch,trim=trim)
		o <- order(x@raw[is.used,ch],decreasing=FALSE)
		cexs <- rep(0.5,length(o))
	}else{
		cols <- x@clusters+1
		cols <- cols[!is.na(cols)]
		cexs <- c(0.3,0.5)[(cols>1)+1]
		uniq.cols <- c("grey70",brewer.pal(11,"Set1"))
		cols <- uniq.cols[cols]

		o <- order(x@clusters+1,decreasing=FALSE)
	}

	stopifnot(nrow(xy)==length(cols))

	points(xy[o,],pch=pch,col=cols[o],cex=cexs[o])#,...)

	if(clust && label){
		clusts <- x@clusters[!is.na(x@clusters)]
		xs <- tapply(xy$x,clusts,median)[-1]
		ys <- tapply(xy$y,clusts,median)[-1]
		text(xs,ys,seq(xs),cex=2)
	}
})
