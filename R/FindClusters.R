#' @export
setGeneric("FindClusters", function(x,...) standardGeneric("FindClusters"))

#' @rdname FindClusters
#' @export
setMethod("FindClusters", "matrix",
          function(x,k.param,
                   initial.membership=NULL,node.sizes=NULL,resolution=0.8,algorithm=1,
                   rsc,with.labels=FALSE,...){   
  ## fin neighbors
  g <- Seurat::FindNeighbors(
    object = x,
    k.param = k.param)
  
  cls <- Seurat::FindClusters(
    object = g$snn,
    initial.membership = initial.membership,
    node.sizes = node.sizes,
    resolution = resolution,
    algorithm = algorithm)[[1]]
  
  return(cls)
  }
)

#' @rdname FindClusters
#' @export
setMethod("FindClusters", "Cycif",
  function(x,ld_name,k.param = 20,
           initial.membership=NULL,node.sizes=NULL,resolution=0.8,algorithm=1,...){
    if(missing(ld_name)){
      stop("'ld_name' should be specified.")
    }else if(!ld_name %in% ld_names(cs1)){
      stop("Specified 'ld_name' does not exist.")
    }
    
    ## subsetting the expression matrix
    ld <- ld_coords(cs1,ld_name)
    used.cts <- ld@used.cts
    this.abs <- ld@used.abs
    
    e <- exprs(cs1,type="logTh_normalized")
    e1 <- data.matrix(e[is.used,this.abs])
    cls <- FindClusters(e1,
                        k.param = k.param,
                        initial.membership = initial.membership,
                        node.sizes = node.sizes,
                        resolution = resolution,
                        algorithm = algorithm)
    x@ld_coords[[ld_name]]@clusters <- cls
    return(x)
  }
)    

#' @rdname FindClusters
#' @export
setMethod("FindClusters", "CycifStack",
  function(x,ld_name,k.param = 20,
           initial.membership=NULL,node.sizes=NULL,resolution=0.8,algorithm=1,...){
    if(missing(ld_name)){
      stop("'ld_name' should be specified.")
    }else if(!ld_name %in% ld_names(cs1)){
      stop("Specified 'ld_name' does not exist.")
    }
    
    ## subsetting the expression matrix
    ld <- ld_coords(cs1,ld_name)
    used.cts <- ld@used.cts
    this.abs <- ld@used.abs
    
    e <- exprs(cs1,type="logTh_normalized")
    e1 <- data.matrix(e[is.used,this.abs])
    cls <- FindClusters(e1,
                 k.param = k.param,
                 initial.membership = initial.membership,
                 node.sizes = node.sizes,
                 resolution = resolution,
                 algorithm = algorithm)
    x@ld_coords[[ld_name]]@clusters <- cls
    return(x)
  }
)    
# rsc,with.labels=FALSE,
#       uniq.cols <- RColorBrewer::brewer.pal(8,"Pastel1")
#       cols <- uniq.cols[as.numeric(cls)]
#       plot(ld@ld_coords,col=cols,pch=20)
#       mids <- sapply(ld@ld_coords,function(x){
#         tapply(x,cls,median)
#       })
#       text(mids,as.character(seq(nrow(mids))),cex=2)
#     }else{
#       uniq.cols <- RColorBrewer::brewer.pal(8,"Dark2")
#       cols <- uniq.cols[as.numeric(cls)]
#       plot(ld@ld_coords,col=cols)
#     }
#      
#     o <- order(cls,decreasing=T)
#     tmp <- oe1 <- e1[o,]
#     ocls <- RColorBrewer::brewer.pal(8,"Set1")[cls[o]]
#     smpls.all <- sub("\\..+","",rownames(oe1))
#     pts.f <- factor(smpls.all,levels=sort(unique(smpls.all)))
#     
#     rownames(tmp) <- rep("",nrow(tmp))
#     heatmap3(tmp,Colv=NA,Rowv=NA,scale="none",RowSideColors=rsc)
#             
# })
# 
# pd <- pData(cs1) %>% 
# select(id,Patient.ID,gBRCA.status,BOR,TimePoint)
# npts <- length(unique(pd$Patient.ID))
# nbrca <- length(unique(pd$gBRCA.status))
# nbor <- length(unique(pd$BOR))
# ntp <- length(unique(pd$TimePoint))
# pd <- pd %>%
# mutate(pts.col = brewer.pal(npts,"Set3")[factor(Patient.ID)]) %>%
# mutate(brca.col = c("black","lightyellow")[factor(gBRCA.status)]) %>%
# mutate(bor.col = brewer.pal(nbor,"Set1")[c(1,3,2)][factor(BOR)]) %>%
# mutate(tp.col = brewer.pal(ntp,"Accent")[factor(TimePoint)])
# 
# pd1 <- data.frame(
# clust = ocls,
# id = smpls.all
# ) %>%
# left_join(pd,by="id") %>%
# rename(smpl = id)
# pd2 <- pd1 %>%  select(clust,pts.col,brca.col,bor.col,tp.col)
# names(pd2) <- c("Clust","Patient","BRCA","BOR","TimePoint")
# pd2 <- as.matrix(pd2)