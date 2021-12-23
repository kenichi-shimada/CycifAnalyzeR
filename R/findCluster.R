# setGeneric("findClusters", function(x,...) standardGeneric("findClusters"))
# setMethod("findClusters", "Cycif",
# 	function(x,method="umap",minPts=30,
# 		plot=c("density","contour","cluster"),...){

# 	n.cells <- x@n_cells
# 	ru <- x@ld_coords

# 	# hdb_obj_real <- map(seq(40, 50, by = 10), ~ hdbscan_by_minPts(minPts = .x, dt = ru))

# 	# hdb_obj_real_res <- map_df(hdb_obj_real, get_hdbscan_result)

# 	# check_stability_score_by_minPts(hdb_obj_real_res)
# 	# check_member_prob_score_by_minPts(hdb_obj_real_res)
# 	# check_outlier_score_by_minPts(hdb_obj_real_res)


# 	cl <- largeVis::hdbscan(ru, minPts = 40)

# 	plot(ru, col=cl$cluster+1, pch=20)
# 	row.inds <- as.numeric(rownames(ru))

# 	diff.ru <- diff(range(ru))
# 	bw <- diff.ru/bw.denom

# 	dens <- KernSmooth::bkde2D(ru,bandwidth=bw,gridsize=c(gsize,gsize))
# 	names(dens) <- c("x","y","z")

# 	##
# 	imp.dens <- fields::interp.surface(dens,ru)
# 	thres <- quantile(imp.dens,dens.thres)

# 	included <- imp.dens > thres
# 	if(sum(included)>n.samples){
# 		inc <- sample(which(included),n.samples)
# 		included <- seq(included) %in% inc
# 	}
# 	db <- dbscan::hdbscan(ru[included,], minPts = minPts)

# 	# plot(ru[included,],col=db$cluster+1,pch=20)
# 	clusts <- rep(NA,n.cells)
# 	clusts.1 <- rep(0,length(included))
# 	clusts.1[included] <- db$cluster
# 	clusts[row.inds] <- clusts.1

# 	x@clusters <- clusts

# 	plot <- plot[1]

# 	if(plot=="density"){
# 		nlevels <- 100
# 		n.dens <- factor(round((imp.dens - min(imp.dens))/diff(range(imp.dens))*nlevels),
# 			levels=0:nlevels)
# 		plot(ru,col=rev(terrain.colors(101))[as.numeric(n.dens)],pch=20)
# 	}else if(plot=="contour"){
# 		contour(dens,col=terrain.colors(10))
# 	}else if(plot=="cluster"){
# 		plot(ru[!included,],pch=20,cex=.5,col="grey70",xlim=range(ru$x),ylim=range(ru$y))
# 		points(ru[included,],col=db$cluster+1,pch=20,cex=.5)
# 	}

# 	validObject(x)
# 	return(x)
# })
