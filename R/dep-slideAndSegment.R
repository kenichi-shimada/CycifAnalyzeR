
setGeneric("slidePlot", function(x,...) standardGeneric("slidePlot"))
setMethod("slidePlot", "Cycif",
	function(x,ch="DNA0",interpolate=F,action=c("none","draw","crop"),
		params=list(shape="boxROI",center=list(x=10000,y=10000),height=1000,width=1000)	,
		zlim=c(0,256),check=TRUE,...){
		action <- action[1]

		fn <- (x@tif_files %>% filter(channel==ch & dim.num==1))$cropped
		im <- raster::raster(fn)

		e <- extent(im)
		im <- crop(im,e)

		# if(is.null(params$height)) params$height <- 1000
		# if(is.null(params$width)) params$width <- 1000
		if(action=="crop"){
			stopifnot(all(c("center","height","width") %in% names(params)))
			shape <- boxROI(params)
			im <- raster::crop(im,shape)
		}

		zmin <- 0
		zmax <- 256

		zlevels <- seq(zlim[1],zlim[2],length=100)
		zcols <- colorRampPalette(c("black","white"))(99)

		lo.col <- ifelse(check,"blue","black")
		hi.col <- ifelse(check,"red","white")

		if(zlim[1]!=zmin){
			zlevels <- c(zmin,zlevels)
			zcols <- c(lo.col,zcols)
		}
		if(zlim[2]!=zmax){
			zlevels <- c(zlevels,zmax)
			zcols <- c(zcols,hi.col)
		}

		pl <- rasterVis::levelplot(im,interpolate=interpolate,
			margin=list(draw=FALSE),
			at=zlevels,
			scales=list(draw=FALSE),
			col.regions=zcols,...)
		return(pl)
	}
)

draw.box.coords <- function(params=list(center,height=1000,width=1000)){
	center <- params$center
	height <- params$height
	width <- params$width

	xmin <- center$x - width/2
	xmax <- center$x + width/2
	ymin <- center$y - height/2
	ymax <- center$y + height/2
	return(list(x=c(xmin,xmax)[c(1,1,2,2,1)],y=c(ymin,ymax)[c(1,2,2,1,1)]))
}

boxROI <- function(params=list(center,height=1000,width=1000)){
	center <- params$center
	height <- params$height
	width <- params$width

	xmin <- center$x - width/2
	xmax <- center$x + width/2
	ymin <- center$y - height/2
	ymax <- center$y + height/2
	e <- raster::extent(xmin,xmax,ymin,ymax)
	return(e)
}

setGeneric("slideHist", function(x,...) standardGeneric("slideHist"))
setMethod("slideHist", "Cycif",
	function(x,ch="DNA0",col=colorRampPalette(c("black","white"))(50),interpolate=F,
		params=list(shape="boxROI",center,height=1000,width=1000),
		segment=FALSE,...){

		fn <- (x@tif_files %>% filter(channel==ch & dim.num==dimNum))$outfile
		im <- raster::raster(fn)

		hist(im)
		if(is.null(params$height)) params$height <- 1000
		if(is.null(params$width)) params$width <- 1000
		if(action=="crop"){
			stopifnot(all(c("center","height","width") %in% names(params)))
			shape <- boxROI(params)
			im <- raster::crop(im,shape)
		}

		# par(bg="black",fg="white")
		plot(im,col=col,axes=F,interpolate=interpolate,...)
		box(col='black',lwd=2)
		lines(c(0,0,1,1,0)*ncol(im),c(0,1,1,0,0)*nrow(im),col="white")

	}
)

setGeneric("addTiff", function(x,...) standardGeneric("addTiff"))
setMethod("addTiff", "Cycif", function(x,indir,outdir,test=TRUE,split=FALSE,overwrite=FALSE){
		sn <- names(x)
		pf <- paste0(indir,"/",sn,".ome.tif")

		stopifnot(file.exists(pf))
		if(split){
			stopifnot(dir.exists(outdir))
		}

		dna.list <- names(x@dna)
		abs.list <- matrix(as.character(x@abs_list$ab),nrow=3)
		all.list <- as.vector(rbind(dna.list,abs.list))

		options(warn=-1)
    info <- capture.output(rgdal::GDALinfo(pf))
    options(warn=0)

    sds <- sub(".+\\((.+)\\).*","\\1",info[grep("^SUBDATASET_[0-9]+_NAME",info,value=F)+1])
    # sds <- gsub(" ","",sds)
    fac <- factor(sds,levels=unique(sds))

    dims <- sub("B$","",levels(fac))
    d2 <- t(sapply(strsplit(dims,split="[PL] x "),as.numeric))[,1:2]
    df2 <- data.frame(x=d2[,1],y=d2[,2],level=unique(sds),dim.num=seq(unique(sds)))

	    # stopifnot(nchannels==length(all.list))

	  if(split){
			outfiles <- paste0(outdir,"/",sn,"_",all.list,"_",as.numeric(fac),".tif")
			outfiles <- outfiles[as.numeric(fac)==1]

			if(!all(file.exists(outfiles))){
				cat("Splitting ome.tif into single channels.\n")
				cat("Read ",pf,"\n",sep="")
				img <- tiff::readTIFF(pf,native=TRUE,all=T)
				cat("Write ",length(outfiles)," files.","\n",sep="")
				for(i in seq(outfiles)){
					outfile <- outfiles[i]
					cat(i,". ",outfile,"\n",sep="")
    			if(!file.exists(outfile.cropped) || overwrite){
						tiff::writeTIFF(img[[i]],outfiles[i],compression="LZW",reduce=TRUE)
					}
				}
				cat("\n")
				rm(img)
				gc()
			}

			df <- data.frame(
				outfile=outfiles,
				channel=all.list,
				dim=fac[seq(outfiles)],
				dim.num=1,
				stringsAsFactors=F)
		}else{
			df <- data.frame(
				outfile=pf,
				channel=all.list,
				dim=dims[1],
				dim.num=1,
				stringsAsFactors=F)
			df2 <- df2[1,]
		}

		x@tif_files <- df
		x@dim <- df2

		validObject(x)
		return(x)
	}
)

setGeneric("defineExtent", function(x,...) standardGeneric("defineExtent"))
setMethod("defineExtent", "Cycif",
	function(x){
		outfiles <- (x@tif_files %>% filter(dim.num==1))$outfile
		stopifnot(all(file.exists(outfiles)))
		##
		cat("Checking the quality of the files.\n")
		for(i in seq(outfiles)){
			if(i==1){
				outfile <- outfiles[i]
				cat(i,". ",outfile,"\n",sep="")
				o <- raster::raster(outfile)
				eo1 <- eo <- extent(o)
				print(rasterVis::levelplot(o))
				n1 <- readline(prompt="Crop the image?([Y]/N)")
				while(!grepl("^[Nn]",n1)){
					print(eo1)
					print(rasterVis::levelplot(o))
					# coords <- click(o,xy=TRUE,n=1,type="p",col="red",pch=20)
					eo1@xmin <- as.numeric(readline(prompt="New xmax:"))
					if(eo1@xmin != eo@xmax){
						cat("Cropping...\n")
						o1 <- crop(o,eo1)
						print(rasterVis::levelplot(o1))
					}
					n1 <- readline(prompt="Re-define the new xmax?([Y]/N)")
				}
				if(eo@xmax!=eo1@xmin && !identical(eo,eo1)){
					eo@xmax <- eo1@xmin
				}

				# o.new <- crop(o,eo)
				lwd <- list(x=(eo@xmax-eo@xmin),y=(eo@ymax-eo@ymin))
				box.coords <- list(
					x=c(eo@xmin,eo@xmax)[c(1,1,2,2,1)],
					y=c(eo@ymin,eo@ymax)[c(1,2,2,1,1)])
				print(rasterVis::levelplot(o) + layer(panel.lines(x,y,col=4,lwd=1),data=box.coords))
			}else if(i > 1){
				n <- readline(prompt="See next channel? ([Y]/N):")
				if(!grepl("^[Nn]",n)&&file.exists(outfile <- outfiles[i])){
					cat(i,". ",outfile,"\n",sep="")
					o <- raster::raster(outfile)
					print(rasterVis::levelplot(o) + layer(panel.lines(x,y,col=4,lwd=1),data=box.coords))
				}else{
					break
				}
			}
		}
		print(eo)
		n <- readline(prompt="Are you satisfied with this new ROI?([Y]/n)")
		extent.file <- sub("_DNA0_1.tif$","_extent.rds",outfiles[1])
		if(grepl("^[Nn]",n)){
			stop("Aborted!")
		}else{
			cat("Returning extent object...\n")
			return(eo)
		}
})

setGeneric("applyNewExtent", function(x,...) standardGeneric("applyNewExtent"))
setMethod("applyNewExtent", "Cycif",
	function(x,extent,overwrite=FALSE){
		outfiles <- (x@tif_files %>% filter(dim.num==1))$outfile
		stopifnot(all(file.exists(outfiles)))

		outfiles.cropped <- sub(".tif$","_cropped.tif",outfiles)

		cat("Cropping all the tif files...\n")
		for(i in seq(outfiles)){
			outfile <- outfiles[i]
			outfile.cropped <- outfiles.cropped[i]
			cat(i,". ",outfile.cropped,"\n",sep="")
			if(!file.exists(outfile.cropped) || overwrite){
				o <- raster::raster(outfile)
				o1 <- crop(o,extent) # this create a temporary path to the cropped object
				terra::writeRaster(o1,filename=outfile.cropped,format="GTiff",options=c("PROFILE=BASELINE"),overwrite=overwrite)
			}
		}

		x@tif_files$cropped <- rep(NA,nrow(x@tif_files))
		x@tif_files$cropped[seq(outfiles)] <- outfiles.cropped
		x@dim$x[1] <- extent@xmax - extent@xmin
		x@dim$y[1] <- extent@ymax - extent@ymin
		x@dim$level[1] <- paste0(x@dim$x[1],"P x ",x@dim$y[1],"L x 1B")
		x@dim <- data.frame(
			x=extent@xmax,
			y=extent@ymax,
			level=paste0(extent@xmax,"P x ",extent@ymax,"L x 1B"),
			dim.num=1
		)

		##
		xy <- x@xy_coords
		used <- xy$X_centroid >= extent@xmin & xy$X_centroid <= extent@xmax &
			xy$Y_centroid >= extent@ymin & xy$Y_centroid <= extent@ymax

		x@raw <- x@raw[used,]
		x@dna <- x@dna[used,]
		x@xy_coords <- x@xy_coords[used,]
		x@segment_property <- x@segment_property[used,]
		x@used_cells <- x@used_cells[used,]

		validObject(x)
		return(x)
	}
)

setGeneric("trimTiff", function(x,...) standardGeneric("trimTiff"))
setMethod("trimTiff", "Cycif",
	function(x,overwrite=FALSE){

		outfiles <- (x@tif_files %>% filter(as.numeric(dim)==1))$outfile
		stopifnot(all(file.exists(outfiles)))

		cat("Checking the quality of the files.\n")
		for(i in seq(outfiles)){
			if(i==1){
				outfile <- outfiles[i]
				cat(i,". ",outfile,"\n",sep="")
				o <- raster::raster(outfile)
				eo1 <- eo <- extent(o)
				print(rasterVis::levelplot(o))
				n1 <- readline(prompt="Crop the image?([Y]/N)")
				while(!grepl("^[Nn]",n1)){
					print(eo1)
					print(rasterVis::levelplot(o))
					# coords <- click(o,xy=TRUE,n=1,type="p",col="red",pch=20)
					eo1@xmin <- as.numeric(readline(prompt="New xmax:"))
					if(eo1@xmin != eo@xmax){
						cat("Cropping...\n")
						o1 <- crop(o,eo1)
						print(rasterVis::levelplot(o1))
					}
					n1 <- readline(prompt="Re-define the new xmax?([Y]/N)")
				}
				if(eo@xmax!=eo1@xmin && !identical(eo,eo1)){
					eo@xmax <- eo1@xmin
				}

				# o.new <- crop(o,eo)
				lwd <- list(x=(eo@xmax-eo@xmin),y=(eo@ymax-eo@ymin))
				box.coords <- list(
					x=c(eo@xmin,eo@xmax)[c(1,1,2,2,1)],
					y=c(eo@ymin,eo@ymax)[c(1,2,2,1,1)])
				print(rasterVis::levelplot(o) + layer(panel.lines(x,y,col=4,lwd=1),data=box.coords))
			}else if(i > 1){
				n <- readline(prompt="See next channel? ([Y]/N):")
				if(!grepl("^[Nn]",n)&&file.exists(outfile <- outfiles[i])){
					cat(i,". ",outfile,"\n",sep="")
					o <- raster::raster(outfile)
					print(rasterVis::levelplot(o) + layer(panel.lines(x,y,col=4,lwd=1),data=box.coords))
				}else{
					break
				}
			}
		}
		print(eo)
		n <- readline(prompt="Are you satisfied with this new ROI?([Y]/n)")
		outfiles.cropped <- sub(".tif$","_cropped.tif",outfiles)
		if(grepl("^[Nn]",n)){
			stop("Aborted!")
		}else{
			cat("Cropping all the tif files...\n")
			for(i in seq(outfiles)){
				outfile <- outfiles[i]
				outfile.cropped <- outfiles.cropped[i]
				outfile1 <- sub(".tif",".rds",outfile)
				cat(i,". ",outfile1,"\n",sep="")
				if(overwrite||!file.exists(outfile1)){
					o <- raster::raster(outfile)
					o1 <- crop(o,eo) # this create a temporary path to the cropped object
					terra::writeRaster(o1,filename=outfile.cropped,format="GTiff",options=c("PROFILE=BASELINE"),overwrite=FALSE)
					# o1@file@name <- outfile.cropped
					# saveRDS(o1,file=outfile1)
				}
			}
		}

		x@tif_files$cropped <- rep(NA,nrow(x@tif_files))
		x@tif_files$cropped[seq(outfiles)] <- outfiles.cropped
		x@dim$x[1] <- eo@xmax - eo@xmin
		x@dim$y[1] <- eo@ymax - eo@ymin
		x@dim$level[1] <- paste0(x@dim$x[1],"P x ",x@dim$y[1],"L x 1B")

		validObject(x)
		return(x)
	}
)
