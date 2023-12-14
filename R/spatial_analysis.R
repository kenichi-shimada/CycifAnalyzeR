# plan
# - ratioTumorStroma(): change to 'generic function'
# - make defineTumorBed() available for CycifStack obj
# - make slots for output of defineTumorBed() in Cycif object
# - make slots for output of ratioTumorStroma() so this function doesn't have to run defineTumorBed()
# spatialStat obj (temtative)
# for cycif
# - lst (should have a better name)
# for cycifstack
# - ratio.ts

#_ -------------------------------------------------------

# fun: defineTumorBed Cycif ----
#' Identify Tumor bed based on dbscan clusters
#' @export
setGeneric("defineTumorBed", function(x,...) standardGeneric("defineTumorBed"))

#' @rdname defineTumorBed
#'
#' @param x A Cycif obj
#' @param strict logical. if cell_types should be strict.
#' @param dth radius for dbscan (eps)
#' @param min.pts minimum points for dbscan
#' @param parallel logical. The default is TRUE.
#' @param n.cores number of cores.
#'
#' @importFrom dbscan dbscan
#' @importFrom dbscan frNN
#' @importFrom parallel mclapply
#' @importFrom concaveman concaveman
#' @importFrom sp point.in.polygon
#' @export
#' Violin plots to show protein expressions
#' @export
setGeneric("defineTumorBed", function(x,...) standardGeneric("defineTumorBed"))

#' @export
setMethod("defineTumorBed", "Cycif",
  function(x,strict=FALSE,dth,parallel=TRUE,min.pts=3,n.cores=1,ct_name="default",verbose=TRUE){
    require(dbscan)
    require(RColorBrewer)
    require(concaveman)
    require(sp)
    if(verbose){
      cat("processing defineTumorBed() for ",names(x),"...\n")
    }

    ## coordinates on slides
    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid

    ## cell types
    this.cts <- cell_types(x,strict=strict,ct_name=ct_name)$cell_types
    if(any(levels(this.cts)!="NA")){
      this.cts <- factor(this.cts,levels=c(levels(this.cts),"NA"))
    }
    this.cts[is.na(this.cts)] <- "NA"

    ## identify Tumor and non-Tumor cells within ROIs
    is.tumor <- !is.na(this.cts) & this.cts == "Tumor"
    xy.tumor <- xy[is.tumor,]

    ### find optimal eps (dth) for adjacent cells - from all the cells (not only tumors)
    if(missing(dth)){
      majl <- x@segment_property$MajorAxisLength
      m1sd <- mean(majl) + sd(majl)*2
      dth <- m1sd * 2 # 35
    }

    ## dbscan clustering - only tumor cells
    dbs1 <- dbscan(xy.tumor, eps=dth, minPts = min.pts, weights = NULL, borderPoints = TRUE)
    cls <- dbs1$cluster

    ## Find neighbors using frNN
    frt <- frNN(xy.tumor,eps=dth)
    fr.ids <- sapply(seq(frt$id),function(i){
      tmp <- frt$id[[i]]
      tmp <- tmp[tmp > i] ## this gaurantees unique pairs
    })

    ## Find within tumors
    is.in <- edge.pts <- list()

    nls <- unique(cls)
    nls <- sort(nls[nls!=0])

    if(parallel){
      ## edge.points
      edge.pts <- parallel::mclapply(seq(nls),function(i){
        xyt <- xy.tumor[cls==nls[i],]
        edge.pt <- as.data.frame(concaveman(as.matrix(xyt), concavity = .8, length_threshold = 0))
        return(edge.pt)
      },mc.cores=n.cores)
      is.in <- parallel::mclapply(seq(edge.pts),function(i){
        is.in <- point.in.polygon(xy$X,xy$Y,edge.pts[[i]]$V1,edge.pts[[i]]$V2)==1
        return(is.in)
      },mc.cores=n.cores)
    }else{
      for(i in seq(nls)){
        xyt <- xy.tumor[cls==nls[i],]
        edge.pts[[i]] <- as.data.frame(concaveman(as.matrix(xyt), concavity = .8, length_threshold = 0))
        is.in[[i]] <- point.in.polygon(xy2$X,xy2$Y,edge.pts[[i]]$V1,edge.pts[[i]]$V2)==1
      }
    }

    is.in.1 <- rowSums(do.call(cbind,is.in))>0

    lst <- list(
      dist.th = dth,
      min.pts = min.pts,
      xy = xy, ## all cells
      within.tumor.bed = is.in.1, ## all cells
      edge.pts = edge.pts, ## for non-empty clusters
      cell.types = this.cts, ## all cells - tumors can be retrieved here
      clusters = cls, # only tumors
      neighbors = fr.ids # only tumors
    )
    return(lst)
  }
)

setMethod("defineTumorBed", "CycifStack",
  function(x,strict=FALSE,dth,parallel=TRUE,min.pts=3,n.cores=1,ct_name="default",verbose=TRUE){
    if(missing(dth)){
      dths <- cyApply(x,function(cy){
        majl <- cy@segment_property$MajorAxisLength
        m1sd <- mean(majl) + sd(majl)
        dth <- m1sd * 2 # 35
        return(dth)
      },simplify=T)
      dth <- median(dths) # roughly around 35 ptx = 23 um
    }
    lst <- cyApply(x,function(cy)defineTumorBed(cy,strict=strict,dth=dth,parallel=parallel,min.pts=min.pts,n.cores=n.cores,ct_name=ct_name,verbose=verbose))
    return(lst)
  }
)

#_ -------------------------------------------------------
# fun: ratioTumorStroma Cycif,CycifStack ----
#  ,strict=FALSE,dth,parallel=TRUE,min.pts=3,n.cores=1,ct_name="default",verbose=TRUE
ratioTumorStroma <- function(x,stroma.cts,...){
  if(is(x,"Cycif")){
    n <- names(x)
    cat("processing ratioTumorStroma() for ",n,"...\n")
    lst0 <- defineTumorBed(x,...)

    tab <- table(cell_types(x)$cell_types,lst0$within.tumor.bed)

    if(missing(stroma.cts)){
      uniq.cts <- rownames(tab)
      stroma.cts <- uniq.cts[!uniq.cts %in% c("Tumor","outOfROI")]
    }

    tab1 <- tab[stroma.cts,]
    sum.str <- colSums(tab1)
    r.all <- sum.str["TRUE"]/sum(sum.str)

    r.each <- tab1[,"TRUE"]/rowSums(tab1)
    r.ts <- r.each/r.all

    return(r.ts)
  }else if(is(x,"CycifStack")){
    ratio <- cyApply(x,ratioTumorStroma)
    r.ts <- data.frame(do.call(rbind,ratio))
    return(r.ts)
  }else{
    stop("this function is defined only for a Cycif object")
  }
}

#_ -------------------------------------------------------
# fun: dist2tumor (immature) ----

dist2tgt <- function(this.cts, xy, seed.tum=1,n.tumors=100,summary=TRUE){
  if(sum(this.cts=="Tumor",na.rm=T)==0){
    stop("some of this.cts should be 'Tumor'")
  }

  has.cts <- !is.na(this.cts)
  uniq.cts <- levels(this.cts)

  lst.xys <- lapply(uniq.cts,function(ct)xy[!is.na(this.cts) & this.cts==ct,])
  names(lst.xys) <- uniq.cts

  ## sampling tumor cellsa
  set.seed(seed.tum)
  n.tumors.all <- nrow(lst.xys$Tumor)
  n.tumors <- min(n.tumors,n.tumors.all)
  idx1 <- sample(n.tumors.all,n.tumors,replace=FALSE)

  xys.tumors <- lst.xys$Tumor[idx1,]

  tgts <- names(lst.xys)
  tgts <- tgts[tgts != "Tumor"]
  dist2tumor <- sapply(tgts,function(tgt){
    if(nrow(lst.xys[[tgt]]) == 0){
      return(c(nn.dist=Inf))
    }
    lst.gds <- RANN::nn2(data=lst.xys[[tgt]],query=xys.tumors,k=1,searchtype="standard") #distance from every tumor to cd8t
    if(summary){
      med.d2tum <- median(lst.gds$nn.dists)
      return(med.d2tum)
    }else{
      return(lst.gds$nn.dists)
    }
  })
  return(dist2tumor)
}

# fun: stat.dist2tumor (immature) ----

stat.dist2tgt <- function(this.cts,xy,
                         seed.tum=1,n.tumors=100,summary=TRUE,
                         n.iter=1e4,seed.iter=123){

  i1 <- which(this.cts!="outOfROI")
  this.cts <- this.cts[i1]
  xy <- xy[i1,]

  ## distance - original
  cts.orig <- this.cts

  orig.d <- dist2tgt(cts.orig,xy,n.tumors=n.tumors,seed.tum=seed.tum,summary=summary)

  ## sampling
  ds <- do.call(rbind,parallel::mclapply(seq(n.iter),function(i){
    set.seed(i)
    cts <- sample(this.cts)
    ds <- dist2tgt(cts,xy,n.tumors=n.tumors,seed.tum=seed.tum,summary=summary)
    return(ds)
  },mc.cores=15))

  med.ds <- apply(ds,2,median)

  folds <- orig.d/med.ds
  return(folds)
}

#_ -------------------------------

# fun: area (for density) ----
#' @importFrom dbscan dbscan
#' @importFrom parallel mclapply
#' @importFrom concaveman concaveman
#' @importFrom sp point.in.polygon
#' @export
compute.area <- function(x,dth,unit=c("mm2"),plot=TRUE,strict=FALSE,
                         ct_name="default",fn){
  if(!is(x,"Cycif")){
    stop("x should be a Cycif object")
  }

  ## coordinates
  xy <- xys(x)
  ymax <- max(xy$Y_centroid)
  xy$Y_centroid <- ymax - xy$Y_centroid
  xy.sp <- sp::SpatialPoints(xy)

  ##
  ## cell types
  this.cts <- cell_types(x,strict=strict,ct_name=ct_name)$cell_types
  levs <- levels(this.cts)
  levs <- levs[!levs %in% c("NA","outOfROI")]
  this.cts <- factor(this.cts,levels=levs)

  ## find cells within rois
  wr <- x@within_rois
  xy1 <- xy[wr,]
  xy.sp1 <- xy.sp[wr,]
  this.cts1 <- this.cts[wr]

  ## define clusters of tumor chunk using dbscan
  dbs1 <- dbscan::dbscan(xy1, eps=dth, minPts = 3, weights = NULL, borderPoints = TRUE) # min 2 cells together
  cls3 <- dbs1$cluster

  ## Merge overlapped clustersFind tumor borders
  nls3 <- unique(cls3) # unique cluster label
  nls3 <- sort(nls3[nls3!=0])

  mcc <- min(15,length(nls3))
  concs3 <- parallel::mclapply(seq(nls3),function(i){
    xyt <- xy1[cls3==nls3[i],]
    conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = .8, length_threshold = dth))
    names(conc) <- c("X","Y")
    return(conc)
  },mc.cores=mcc)

  # Function to check if a point is inside a polygon using 'point.in.polygon()'
  is_point_inside_polygon <- function(point, polygon) {
    sp::point.in.polygon(point$X, point$Y, polygon$X, polygon$Y)
  }

  if(length(concs3)>1){
    # Identify overlapping clusters and merge their points
    overlaps <- sapply(concs3,function(cluster_i){
      sapply(concs3,function(cluster_j){
        any(is_point_inside_polygon(cluster_j, cluster_i))
      })
    })
    diag(overlaps) <- 0

    # update clusters
    ol <- as.data.frame(which(overlaps==1,arr.ind=T))
    names(ol) <- c("child","parent")
    ol1 <- ol

    nls3.updated <- nls3
    names(nls3.updated) <- nls3

    while(nrow(ol1)){
      pa <- ol1$pa[1]
      ch <- ol1$ch[1]
      nls3.updated[ch] <- pa
      ol1$parent[ol1$parent==ch] <- pa
      ol1$child[ol1$child==ch] <- pa
      ol1 <- ol1[-1,]
    }

    idx <- cls3 >0
    cls3.updated <- cls3
    cls3.updated[idx] <- nls3.updated[cls3[idx]]
    nls3.updated <- unique(nls3.updated)

    concs3.updated <- parallel::mclapply(seq(nls3.updated),function(i){
      xyt <- xy1[cls3.updated==nls3.updated[i],]
      conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = .8, length_threshold = dth))
      names(conc) <- c("X","Y")
      return(conc)
    },mc.cores=7)
    names(concs3.updated) <- nls3.updated
  }else{
    cls3.updated <- cls3
    nls3.updated <- nls3
    concs3.updated <- parallel::mclapply(seq(nls3.updated),function(i){
      xyt <- xy1[cls3.updated==nls3.updated[i],]
      conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = .8, length_threshold = dth))
      names(conc) <- c("X","Y")
      return(conc)
    },mc.cores=7)
    names(concs3.updated) <- nls3.updated
  }

  if(plot){
    png(fn,width=700,height=700)

    par(mar=c(5,5,3,3))
    sp::plot(xy.sp1,col=cls3.updated+1,pch=c(2,20)[(cls3.updated!=0)+1],
         # asp=1,xlim=range(l1$x),ylim=range(l1$y),
         cex=.5,
         main="tumors, clustered (dbscan)")
    sapply(concs3.updated,function(conc)lines(conc))
    dev.off()
  }

  ## compute area
  areas <- sapply(concs3.updated,function(coords){
    polygon <- sp::Polygon(coords)
    return(polygon@area)
  })

  sum.area <- sum(areas) * (0.65 * (10^-3))^2 # unit: mm^2
  return(sum.area)
}

#_ -------------------------------

# fun: computeRCN ----
#' @importFrom dbscan frNN
#' @importFrom data.table := rbindlist
#' @importFrom dplyr %>%
#' @export
computeRCN <- function(x,r,unit=c("pixel","um"),
                     focused.cts,
                     cts.in.rcn,
                     n.sampling,
                     seed=123){ # only for Cycif, and CycifStack
  if(is(x,"Cycif")){
    if(missing(r)){
      stop("radius (r) should be specified")
    }else if(unit[1]=="um"){
      r <- r/0.65
      unit <- "pixel"
    }

    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid
    cts <- cell_types(x)

    wr <- x@within_rois
    xy1 <- xy[wr,]
    cts1 <- cts[wr,]
    lev.cts <- levels(cts1$cell_types)
    lev.cts <- lev.cts[lev.cts != "outOfROI"]

    if(missing(focused.cts)){
      focused.cts <- lev.cts
    }

    if(missing(cts.in.rcn)){
      cts.in.rcn <- lev.cts
    }

    ## frnn object - all within_rois samples
    frnn <- dbscan::frNN(xy1,eps=r,bucketSize=10) # 265793 cells
    frnn.ids <- frnn$id
    n.frnn <- sapply(frnn.ids,length)
    n0 <- which(n.frnn==0)
    idx.df <- data.frame(center_id=rep(seq(frnn.ids),n.frnn),
                         neighbor_id=.Internal(unlist(frnn.ids, FALSE, FALSE))) %>%
      data.table::as.data.table()

    ##
    df1 <- cell_types(x) %>%
      cbind(exprs(x,type="log")) %>%
      filter(cell_types != "outOfROI") %>%
      # tibble::rowid_to_column("center_id") %>%
      data.table::data.table()

    rn <- levels(df1$cell_types)
    csts <- x@cell_types$default@cell_state_def[rn,]

    lst.csts <- lapply(csts,function(x1){
      rn[!is.na(x1) & x1 == "CAN"]
    })
    cst.abs <- names(csts)

    for(ab in cst.abs){
      this.cts <- lst.csts[[ab]]
      df1 <- df1[!cell_types %in% this.cts ,(ab) := NA]
    }

    ##
    df1 <- cbind(idx.df,df1[idx.df$neighbor_id,])
    dt <-  df1[ ,lapply(.SD, mean, na.rm=T), by=center_id,.SDcols=cst.abs]
    dt0 <- data.table::as.data.table(cbind(n0,array(NA,c(length(n0),length(cst.abs)))))
    names(dt0) <- names(dt)

    dt <- data.table::rbindlist(list(dt,dt0))[order(center_id),]
    has.neighbors <- !seq(nrow(dt)) %in% n0

    ## within positive ROIs + has neighbors
    selected.ids1 <- which(cts1$cell_types %in% focused.cts &
                          sapply(frnn$id,length)>0 &
                          sapply(frnn$id,function(ids){
                            any(cts1$cell_types[ids] %in% cts.in.rcn)
                          }) &
                          has.neighbors
                        )

    ## sampling
    set.seed(seed)
    n.this <- length(selected.ids1)
    n.cts <- min(n.sampling,n.this)

    selected.ids2 <- sample(selected.ids1,n.cts) # idx after ROI filter

    ## selected ids among avialble focused celltypes (eg tumor cells)
    sid <- seq(frnn$id) %in% selected.ids2


    ## cell type frequency (not absolute count)
    rcn.freq <- t(sapply(frnn$id,function(ids){ # [c(83538,101787)]
      ct <- cts1$cell_types[ids]
      tab <- table(ct)
      ntab <- tab/sum(tab)
      return(ntab)
    }))
    rcn.count <- t(sapply(frnn$id,function(ids){ # [c(83538,101787)]
      ct <- cts1$cell_types[ids]
      tab <- table(ct)
      return(tab)
    }))

    return(list(within.rois=wr,
                focused.cts=focused.cts,
                cts.in.rcn=cts.in.rcn,
                n.cells.selected=n.cts,
                frnn=frnn,
                exp=dt,
                is.selected=sid,
                rcn.count=rcn.count,
                rcn.freq=rcn.freq))

  }else if(is(x,"CycifStack")){
    cat("Get neighbors ...\n")
    frnn1 <- cyApply(x,function(cy){
      cat(names(cy),"\n")
      computeRCN(x=cy,r=r,unit=c("pixel","um"),
        focused.cts=focused.cts,
        cts.in.rcn=cts.in.rcn,
        n.sampling=n.sampling,
        seed=seed)
    })

    cat("Restructure data ...\n")
    ## assemble frnn
    lst.frnn <- lapply(frnn1,function(fr)fr$frnn)

    ### dists
    dists <- do.call(c,lapply(lst.frnn,function(frnn)frnn$dist))

    ### ids
    n.frnns <- sapply(lst.frnn,function(frnn)length(frnn$id)) # 1325874, all cells, excluding outOfROIs
    n.frnns.pre <- c(0,cumsum(n.frnns)[-length(n.frnns)])
    names(n.frnns.pre) <- names(x)

    ## frnn.ids, frnn.ids1, frnn.tum.ids - list of neighboring cells ids for tumors used for the frnn analysis
    frnn.ids <- lapply(names(x),function(nm){
      x <- lst.frnn[[nm]]
      n.prior <- n.frnns.pre[nm]
      this.ids <- x$id
      new.ids <- lapply(this.ids,function(id){
        new.id <- id + n.prior
        return(new.id)
      })
      return(new.ids)
    })
    frnn.ids1 <- do.call(c,frnn.ids) ## 1325874, now all data are combined - and the indices are after excluding outOfROIs

    ### eps, sort
    eps <- unique(sapply(lst.frnn,function(frnn)frnn$eps))
    sort <- unique(sapply(lst.frnn,function(frnn)frnn$sort))

    frnn <- list(dist=dists,id=frnn.ids1,eps=eps,sort=sort)

    ## within.rois
    within.rois <- do.call(c,lapply(frnn1,function(fr)fr$within.rois))

    ## focused.cts
    focused.cts <- unique(as.vector(sapply(frnn1,function(fr)fr$focused.cts)))

    ## cts.in.rcn
    cts.in.rcn <- frnn1[[1]]$cts.in.rcn

    ## n.cells.selected
    n.cells.selected <- sapply(frnn1,function(fr)fr$n.cells.selected)
    v.smpls <- rep(names(n.cells.selected),n.cells.selected)

    ## is.selected
    mclustda <- list()
    mclustda$sele <- mclustda$all <- list()
    mclustda$sele$is.used <- do.call(c,lapply(frnn1,function(fr)fr$is.selected)) #67119 -> 66904, how did I do that?

    ##
    exps <- data.table::rbindlist(lapply(frnn1,function(fr)as.data.frame(fr$exp)))
    rcn.count <- data.table::rbindlist(lapply(frnn1,function(fr)as.data.frame(fr$rcn.count)))
    rcn.freq <- data.table::rbindlist(lapply(frnn1,function(fr)as.data.frame(fr$rcn.freq)))

    return(list(within.rois=within.rois, ## all cells (including outOfROIs)
                focused.cts=focused.cts, ## 1 or a few
                cts.in.rcn=cts.in.rcn, ## 1 or a few
                n.cells.selected=n.cells.selected,
                v.smpls=v.smpls,
                frnn=frnn,
                mclustda=mclustda, #is.selected=is.selected, <= included in mclustda
                exp=exps,
                rcn.count=rcn.count,
                rcn.freq=rcn.freq))
  }else{
    stop("x should be Cycif or CycifStack objects")
  }
}

#_ -------------------------------

# fun: rcnClust ----
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @export
rcnClust <- function(frnn,
                     g=50,
                     seed=123,
                     sort.by="CD8T",
                     sort.type=c("freq","count"),
                     sort.smpls=c("all","selected"),
                     data.type=c("ct_exp","ct"),
                     extrapolate=FALSE,
                     mc.cores=1){
  mclustda <- frnn$mclustda

  exps <- as.matrix(frnn$exp)[,-1]
  exps.imp <- imputeData(exps)

  this.cts <- cts.in.rcn

  mat.count.all <- as.matrix(frnn$rcn.count)[,this.cts]
  mat.freq.all <- t(apply(as.matrix(frnn$rcn.freq)[,this.cts],1,function(x)x/sum(x)))

  is.selected <- mclustda$sele$is.used

  if(missing(data.type)){
    data.type="ct_exp"
  }

  if(data.type=="ct"){
    mat.count.sele <- mat.count.all[is.selected,]
    mat.freq.sele <- mat.freq.all[is.selected,]
  }else if(data.type=="ct_exp"){
    mat.count.sele <- cbind(mat.count.all,exps.imp)[is.selected,]
    mat.freq.sele <- cbind(mat.freq.all,exps.imp)[is.selected,]
  }

  "ct_exp"

  ## clustering & classification
  set.seed(seed)

  cat("Clustering with Mclust ...\n")
  mem.ori <- mclust::Mclust(data=mat.freq.sele,G=g,modelNames="EII")$classification

  cat("Training MclustDA ...\n")
  mc1 <- mclust::MclustDA(data=mat.freq.sele,class=mem.ori,
                          G=g,modelNames="EII",modelType = "EDDA") # 100, EII
  mem.sele <- factor(predict(mc1)$classification)
  g1 <- nlevels(mem.sele)

  ## mean counts & frequencies

  mean.count.sele <- sapply(seq(g1),function(i)colMeans(mat.count.sele[mem.sele==i,]))
  mean.freq.sele <- sapply(seq(g1),function(i)colMeans(mat.freq.sele[mem.sele==i,]))
  colnames(mean.count.sele) <- colnames(mean.freq.sele) <- seq(g1)

  if(extrapolate){
    ## extrapolate clusters
    is.available <- !apply(mat.freq.all,1,function(x)any(is.na(x))) # 8587
    mat.count.all1 <- mat.count.all[is.available,]
    mat.freq.all1 <- mat.freq.all[is.available,]

    cat("Applying MclustDA to the entire data ...\n")
    mc.idx <- sort(rep(seq(mc.cores),length=nrow(mat.freq.all1)))

    mem.all <- parallel::mclapply(seq(mc.cores),function(i){
      predict(mc1,newdata=mat.freq.all1[mc.idx==i,frnn$cts.in.rcn])$classification
    },mc.cores=mc.cores)
    mem.all <- do.call(c,mem.all)

    ## update labels
    mem.all <- factor(mem.all)
    g1 <- nlevels(mem.all)
    levels(mem.all) <- seq(g1)

    ## mean counts & frequencies
    mean.count.all <- sapply(seq(g1),function(i)colMeans(mat.count.all1[mem.all==i,]))
    mean.freq.all <- sapply(seq(g1),function(i)colMeans(mat.freq.all1[mem.all==i,]))

    colnames(mean.count.all) <- colnames(mean.freq.all) <- seq(g1)
  }

  ## Sort clusters based on frequency of a cell type (CD8T by default)
  if(missing(sort.type)){
    sort.type <- "freq"
  }

  if(missing(sort.smpls)){
    if(extrapolate){
      sort.smpls <- "all"
    }else{
      sort.smpls <- "selected"
    }
  }

  if(sort.type == "freq" & sort.smpls == "all"){
    mean.freq <- mean.freq.all
  }else if(sort.type == "freq" & sort.smpls == "selected"){
    mean.freq <- mean.freq.sele
  }else if(sort.type == "count" & sort.smpls == "all"){
    mean.freq <- mean.count.all
  }else if(sort.type == "count" & sort.smpls == "selected"){
    mean.freq <- mean.count.sele
  }

  if(sort.by=="CD8T"){
    o <- order(mean.freq["CD8T",],decreasing=T)
    mean.freq.sele <-  mean.freq.sele[,o]
    mean.count.sele <-  mean.count.sele[,o]
    mem.sele <- as.numeric(factor(as.character(mem.sele),levels=as.character(seq(ncol(mean.freq))[o])))
    colnames(mean.count.sele) <- colnames(mean.freq.sele) <- paste0("Rcn",seq(g1))

    mclustda$sele <- list(is.used = mclustda$sele$is.used,
                          mem = mem.sele,
                          mean.freq = mean.freq.sele,
                          mean.count = mean.count.sele)

    if(extrapolate){
      mean.freq.all <-  mean.freq.all[,o]
      mean.count.all <-  mean.count.all[,o]
      mem.all <- as.numeric(factor(as.character(mem.all),levels=as.character(seq(ncol(mean.freq))[o]))) ## mem redefined
      colnames(mean.count.all) <- colnames(mean.freq.all) <- paste0("RCN",seq(g1))

      mclustda$all <- list(is.used = is.available,
                           mem = mem.all,
                           mean.freq = mean.freq.all,
                           mean.count = mean.count.all)
    }
  }else{
    stop("Define 'sort.by' first")
  }

  ##
  mclustda$g <- g1

  frnn$mclustda <- mclustda

  return(frnn)
}


#_ -------------------------------
# fun: meanExpRCN ----
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist
#' @importFrom tibble rowid_to_column
#' @importFrom dplyr %>%
#' @export
meanExpRCN <- function(x,
                       frnn,
                       focused.cts="Tumor",
                       cts.in.rcn=levels(cell_types(x)$cell_types)[1:10],
                       per.ct=TRUE,
                       extrapolate=TRUE){

  g <-   frnn$mclustda$g

  lin.abs <- names(x@cell_types$default@cell_lineage_def[-(1:2)])
  cst.abs <- names(x@cell_types$default@cell_state_def)
  all.abs <- unique(c(lin.abs,cst.abs))

  if(extrapolate){
    mclustda <- frnn$mclustda$all
  }else{
    mclustda <- frnn$mclustda$sele
  }

  # clusts <- factor(mclustda$mem)
  clusts <- factor(paste0("RCN",mclustda$mem),labels=paste0("RCN",seq(g)))
  df1 <- cell_types(x) %>%
    filter(frnn$within.rois) %>%
    tibble::rowid_to_column("idx") %>% ## idx is numbered within 'within.rois'
    dplyr::left_join(
      data.frame(tcn = clusts) %>%
        mutate(idx=which(mclustda$is.used)),by="idx") %>%
    arrange(idx) %>%
    dplyr::left_join(exprs(x,type="log") %>%
                filter(frnn$within.rois) %>%
                tibble::rowid_to_column("idx"),by="idx") %>%
    rename(idx.all="idx")## 1325874, within ROIs

  idx.nonna.tum <-!is.na(df1$tcn) & df1$cell_types %in% focused.cts # 471910 / 1325874
  df.tum <- df1 %>%
    filter(idx.nonna.tum) %>%
    tibble::rowid_to_column("idx.tum") # 471910

  ## lst.frnn: convert ids in each sample to ids in all samples
  frnn.ids <- frnn$frnn$id[!is.na(df1$tcn)] # 1325311: 563 don't have proper RCNs.
  frnn.tum.ids <- frnn$frnn$id[idx.nonna.tum] # 471910

  ## all unique ids per tcn cluster
  if(per.ct){
    lst.mean.exp <- tapply(frnn.tum.ids,df.tum$tcn,function(this.ids){
      this.ids1 <- unique(sort(unlist(this.ids)))
      df.tmp <- df1[this.ids1,] %>%
        group_by(cell_types) %>%
        summarize_at(vars(!!!syms(all.abs)),mean,na.rm=TRUE)
      return(df.tmp)
    })
    df.exp <- as.data.frame(data.table::rbindlist(lapply(seq(lst.mean.exp),function(i){
      lst.mean.exp[[i]] %>%
        mutate(tcn=paste0("RCN",i))
    }))) %>%
      mutate(tcn = factor(tcn,levels=paste0("RCN",seq(g))))
    mes <- lapply(all.abs,function(ab){
      df.exp %>%
        select(cell_types,tcn,!!sym(ab)) %>%
        tidyr::spread(tcn,!!sym(ab))
    })
    names(mes) <- all.abs
    return(mes)
  }else{
    lst.mean.exp <- tapply(frnn.tum.ids,df.tum$tcn,function(this.ids){
      this.ids1 <- unique(sort(unlist(this.ids)))
      df.tmp <- df1[this.ids1,] %>%
        group_by() %>%
        summarize_at(vars(!!!syms(all.abs)),mean,na.rm=TRUE)
      return(df.tmp)
    })

    mes <- as.data.frame(do.call(rbind,lst.mean.exp))
    rownames(mes) <- names(lst.mean.exp)
    colnames(mes) <- all.abs
    return(mes)
  }
}

#_ -------------------------------
# metaRcnClust(x,frnn){
#   if(!is(x,"Cycif")){
#     stop("x should be a cycif object")
#   }
# }



# compute nubmer of tumor cells, T cells, etc
# spatial information - connectivity

#_ -------------------------------

# connectivity <- function(cs7=cs7,
#                        frnn=frnn,
#                        focused.cts="Tumor",
#                        cts.in.tcn=levels(df1$cell_types)[1:10],
#                        per.ct=TRUE,
#                        =TRUE){
#
#
#   if(extrapolate){
