#_ -------------------------------

# fun: computeCN ----

#' Compute Recurrent Cell Neighbors (RCN) for Cycif or CycifStack Objects
#'
#' This function computes the Recurrent Cell Neighbors (RCN) for Cycif or CycifStack objects.
#' RCN measures the relative frequencies of neighboring cell types around each cell within specified
#' radius 'r'. The RCN analysis can be performed on a single Cycif object or across a CycifStack object.
#'
#' @param x A Cycif or CycifStack object.
#' @param r The radius within which neighboring cells are considered (in 'unit').
#' @param unit The unit of measurement for the radius ('pixel' or 'um'). If 'um' is specified,
#' the radius 'r' will be converted to pixels based on the assumed resolution (0.65 um per pixel).
#' @param cts.in.center A character vector specifying the cell types around which RCN is computed. If not is specified, all available cts are used.
#' @param cts.in.rcn A character vector specifying the cell types to consider when computing RCN values. If not is specified, all avaialble cts are used.
#' @param n.sampling The number of cells to randomly sample for RCN analysis.
#' @param seed The random seed for reproducibility.
#'
#' @return A CellNeighborhood object containing the following components:
#' - 'within.rois': A logical vector indicating whether each cell is within a region of interest (ROI). The length is the same as the number of cells in the dataset.
#' - 'cts.in.rcn': A character vector specifying the cell types considered when computing cell neighbors, cell type frequency, and expressions.
#' - 'n.cells.selected': The number of cells selected for RCN analysis, which is the smaller of `n.sampling` and the number of cells within ROIs.
#' - 'frnn': A list with recurrent Neighborhood information, including 'dist', 'id', 'eps', and 'sort'. The length of 'dist' and 'id' is the same as 'n.cells.selected'.
#' - 'cn_exp': A data frame containing expression data for selected cells.
#' - 'is.selected': A logical vector indicating whether each cell is selected for RCN analysis. The sum of the vector is the same as 'n.cells.selected'.
#' - 'rcn.count': A data frame containing the counts of neighboring cell types.
#' - 'rcn.freq': A data frame containing the relative frequencies of neighboring cell types.
#'
#' @details The RCN analysis is performed as follows:
#' 1. The function first identifies the cells that are within the ROIs.
#' 2. It then computes the Recurrent Neighborhood (frNN) for the selected cells using the 'dbscan::frNN' function.
#' 3. It then computes the RCN values for each cell type based on the relative frequencies of neighboring cell types.
#' 4. The function returns a list containing the RCN values for each cell type.
#' The RCN analysis can be performed on a single Cycif object or across a CycifStack object.
#' If the input is a CycifStack object, the RCN analysis is performed on each Cycif object in the stack.
#'
#' @seealso \code{\link{cyApply}}, \code{\link{dbscan::frNN}}
#'
#' @importFrom dbscan frNN
#' @importFrom data.table := rbindlist setDT melt
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tibble rowid_to_column
#' @importFrom dplyr na_if
#'
#' @rdname computeCN
#' @export
setGeneric("computeCN", function(x,...) standardGeneric("computeCN"))

#' @rdname computeCN
#' @export
setMethod("computeCN", "Cycif",
  function(x,r_um=20,k=20,
           type=c("knn","frnn"),
           used.cts,
           n.sampling=1000,
           ct_name="default",
           seed=123,
           off_target = NULL, #list(Tumor_panCK = c("GrB","PD1"), Fibro = "gH2AX")
           off_target_mode = c("na","zero")
           ){ # only for Cycif, and CycifStack
    ## find cells within rois - roi is a circle with a fixed radius (r_um)
    if(missing(type)){
      type <- "frnn"
    }else if(!type %in% c("knn","frnn")){
      stop("type must be either 'knn' or 'frnn'")
    }

    if(type=="frnn"){
      cat("compute frnn\n")
      r <- r_um/0.65 # unit converted to pixel
    }else{
      cat("compute knn\n")
    }

    smpl <- names(x)
    cts <- cell_types(x,ct_name=ct_name)

    ## coordinates: x => xy
    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid
    # xy.sf <- st_as_sf(xy, coords = c("X_centroid", "Y_centroid"), crs = NA)

    wr <- within_rois(x)
    ## ---- cell types to consider (used.cts) - "OutOfROI" is removed for this subset of cells ----
    if(missing(used.cts)){
      lev.cts <- levels(cts$cell_types)
      lev.cts <- lev.cts[lev.cts != "outOfROI"]
      used.cts <- lev.cts
    }

    is.used <- cts$cell_types %in% used.cts

    # ---- expression matrix per cell type per CellNeighborhood ----
    xy1 <- xy[is.used,]
    cts1 <- cts[is.used,]
    df1 <- exprs(x,type="log")[is.used,]

    ## convert xy1, cts1, df1 to data.table
    data.table::setDT(xy1)
    data.table::setDT(cts1)
    data.table::setDT(df1)

    ## ---- compute frNN for all cells within ROIs and convert them to frNN object.  ----
    if(type=="frnn"){
      nn <- dbscan::frNN(xy1,eps=r,bucketSize=10)
      nn.ids <- nn$id <- lapply(seq_along(nn$id), function(i) c(i, nn$id[[i]]))
    }else if(type=="knn"){
      nn <- dbscan::kNN(xy1,k=k,bucketSize=10)
      nnids <- cbind(seq(nrow(nn$id)),nn$id)
      colnames(nnids) <- 0:k
      nn.ids <- nn$id <- apply(nnids,1,function(x)x,simplify=FALSE)
      nn$dist <- apply(nn$dist,1,function(x)c(0,x),simplify=FALSE)
    }

    # Add a row number column to xy1 and cts1 for joining
    xy1[, rn := .I]
    cts1[, rn := .I]
    df1[, rn := .I]

    ## find out which cell types express which cell state markers
    csts <- x@cell_types[[ct_name]]@cell_state_def
    if(colnames(csts)[1] != "cell_types"){
      # stop("The first column of the cell state marker definition table must be 'cell_types'")
      csts <- csts %>% rownames_to_column("cell_types")
    }
    csts <- csts %>% filter(cell_types %in% used.cts)
    protein_columns <- names(csts)
    protein_columns <- protein_columns[protein_columns != "cell_types"]

    # Convert csts to data.table and melt it to long format
    csts_dt <- data.table::melt(data.table::setDT(csts), id.vars = "cell_types", variable.name = "ab", value.name = "expression")

    # Filter for cell types that do not express each antibody (i.e., expression is NA)
    csts_dt <- csts_dt[is.na(expression)]

    # Join xy1, cts1, and df1 by row number
    combined_df <- xy1[cts1, on = "rn"][df1, on = "rn"]

    # [discontinued] Loop through each antibody and set expression to NA for cell types that do not express it
    if(0){ ## decided not to turn them as NAs as I want to compute expression per CN, irrespective of cell type composition.
      for (ab1 in as.character(unique(csts_dt$ab))) {
        non_expressing_cts <- csts_dt[ab == ab1, cell_types]
        combined_df[.(non_expressing_cts), (ab1) := NA, on = .(cell_types)]
      }
    }

    ## ### NEW: Drop off-target signals BEFORE neighborhood aggregation
    if (!is.null(off_target) && length(off_target) > 0) {
      # Normalize names to character vectors
      # Keep only valid (existing) cell types and protein columns
      valid_cts <- unique(combined_df$cell_types)
      valid_abs <- protein_columns

      # Build a two-column data.table of (cell_types, ab) pairs to zero/NA-out
      ot_pairs <- rbindlist(lapply(names(off_target), function(ct) {
        if (!(ct %in% valid_cts)) return(NULL)
        abs_here <- intersect(as.character(off_target[[ct]]), valid_abs)
        if (length(abs_here) == 0) return(NULL)
        data.table::data.table(cell_types = ct, ab = abs_here)
      }), use.names = TRUE, fill = TRUE)

      if (nrow(ot_pairs) > 0) {
        # Efficiently set each specified antibody to NA (or 0) for rows of that cell type
        # Use data.table::set for speed in a loop over antibodies (column-wise)
        if (off_target_mode == "na") {
          for (ab1 in unique(ot_pairs$ab)) {
            ct_to_null <- ot_pairs[cell_types %in% valid_cts & ab == ab1, unique(cell_types)]
            if (length(ct_to_null) == 0) next
            idx <- combined_df$cell_types %in% ct_to_null
            # set NA
            set(combined_df, which(idx), ab1, NA_real_)
          }
        } else { # "zero"
          for (ab1 in unique(ot_pairs$ab)) {
            ct_to_zero <- ot_pairs[cell_types %in% valid_cts & ab == ab1, unique(cell_types)]
            if (length(ct_to_zero) == 0) next
            idx <- combined_df$cell_types %in% ct_to_zero
            set(combined_df, which(idx), ab1, 0.0)
          }
        }
      }
    }

    # tapply(combined_df$PD1,combined_df$cell_types,mean,na.rm=T) # sanity check

    # Create a neighborhood data table from nn
    nmat <- data.table::data.table(
      cell_id = rep(seq_along(nn$id), lengths(nn$id)),
      neighbor_id = unlist(nn$id, use.names = FALSE)
    )

    # Join with combined_df to get cell type and expression data for each neighbor
    neighborhood <- nmat[
      combined_df, on = .(neighbor_id = rn), nomatch = 0
    ]

    # Calculate the average expression per cell type in each neighborhood
    # Assuming protein columns in combined_df are named "protein1", "protein2", etc.

    exp_per_ct_cn <- neighborhood[, lapply(.SD, mean, na.rm = TRUE), by = .(cell_id, cell_types), .SDcols = protein_columns]
    exp_per_cn <- neighborhood[, lapply(.SD, mean, na.rm = TRUE), by = .(cell_id), .SDcols = protein_columns]

    # For exp_per_cn: ensure sorted by cell_id
    data.table::setorder(exp_per_cn, cell_id)


    if(type=="knn"){
      nn$eps <- -1
    }else if (type=="frnn"){
      nn$k <- -1
    }

    ##  ---- convert nn to nn object (not necessary?)  ----
    nn1 <- new("NN",
               type = type,
               dist = nn$dist,
               id = nn$id,
               k = nn$k,
               eps = nn$eps,
               sort = nn$sort)

    ##  ---- within positive ROIs + has neighbors ----
    n.nn <- lengths(nn1@id) # number of neighbors - the same as sum(wr)

    selected.ids1 <- which(n.nn>0 &
                           sapply(nn.ids,function(ids){
                             any(cts1$cell_types[ids] %in% used.cts)
                           })) # cell indices are after ROI filter

    n.this <- length(selected.ids1)
    n.cts <- min(n.sampling,n.this)

    ## sampling
    set.seed(seed)

    # Randomly sample n.cts cells from selected.ids1 (idx after ROI)
    selected.ids2 <- sample(selected.ids1,n.cts)

    ## selected ids among available focused celltypes (eg tumor cells)
    is.selected <- seq(nn.ids) %in% selected.ids2

    ## cell type count and frequency in each RCN
    rcn.freq <- t(sapply(nn1@id,function(ids){
      ct <- cts1$cell_types[ids]
      tab <- table(ct)
      ntab <- tab/sum(tab)
      return(ntab)
    }))

    rcn.count <- t(sapply(nn1@id,function(ids){
      ct <- cts1$cell_types[ids]
      tab <- table(ct)
      return(tab)
    }))

    # For rcn.count: make sure rownames are 1..N in the same order
    rownames(rcn.freq) <- seq_len(nrow(rcn.freq))
    rownames(rcn.count) <- seq_len(nrow(rcn.count))

    if(type=="knn"){
      r <- sapply(nn$dist,function(x)x[k+1])
    }
    rcn.dens <- rcn.count/(pi*r^2)

    cn <- new("CellNeighborhood",
              within.rois=wr, # logical, within ROIs
              used.cts=used.cts, # character, cell types to consider
              n.cells.selected=n.cts, # integer, number of cells selected
              is.selected = is.selected,
              smpls = smpl,
              nn=nn1,
              n.neighbors = n.nn,
              exp.per.ct.cn=exp_per_ct_cn,
              exp.per.cn=exp_per_cn,
              rcn.count=rcn.count,
              rcn.dens=rcn.dens,
              rcn.freq=rcn.freq)
    return(cn)
})

#' @rdname computeCN
#' @export
setMethod("computeCN", "CycifStack",
  function(x,r_um = 20,
           used.cts,
           n.sampling,
           seed=123){ # only for Cycif, and CycifStack

    cat("Get neighbors ...\n")
    nn1 <- cyApply(x,function(cy){
      cat(names(cy),"\n")
      computeCN(x=cy,r=r,unit=c("pixel","um"),
        used.cts=used.cts,
        n.sampling=n.sampling,
        seed=seed)
    })

    cat("Restructure data ...\n")

    ## within.rois
    within.rois <- unlist(lapply(nn1,function(fr)fr@within.rois)) # same as nCells() for each sample

    ## n.cells.selected
    n.cells.selected <- sapply(nn1,function(fr)fr@n.cells.selected)

    # is.selected
    is.selected <- unlist(sapply(nn1, function(fr)fr@is.selected))

    if(sum(is.selected) != sum(n.cells.selected)){
      stop("is.selected and n.cells.selected are not consistent")
    }

    ## used.cts
    used.cts <- nn1[[1]]@used.cts

    ## smpls
    n.smpls <- sapply(nn1,function(fr)sum(fr@within.rois))
    smpls <- rep(names(n.smpls),n.smpls)

    ## nn
    lst.nn <- lapply(nn1,function(fr)fr@nn)

    ### nn@dists
    nn.dists <- do.call(c,lapply(lst.nn,function(fr)fr@dist))

    ### nn@id - combine indices so they can specify selected cells in the entire dataset
    n.nns <- sapply(lst.nn,function(nn)length(nn@id)) # 1325874, all cells, excluding outOfROIs
    n.nns.pre <- c(0,cumsum(n.nns)[-length(n.nns)])
    names(n.nns.pre) <- names(x)

    ## nn.ids, nn.ids1, nn.tum.ids - list of neighboring cells ids for tumors used for the nn analysis
    nn.ids <- lapply(names(x),function(nm){
      x <- lst.nn[[nm]]
      n.prior <- n.nns.pre[nm]
      this.ids <- x@id
      new.ids <- lapply(this.ids,function(id){
        new.id <- id + n.prior
        return(new.id)
      })
      return(new.ids)
    })
    nn.ids1 <- do.call(c,nn.ids) ## 1325874, now all data are combined - and the indices are after excluding outOfROIs

    ### nn@eps
    ### nn@sort
    eps <- unique(sapply(lst.nn,function(nn)nn@eps))
    sort <- unique(sapply(lst.nn,function(nn)nn@sort))

    ### assemble nn
    nn1 <- new("NN",
               type = type,
               dist = nn$dist,
               id = nn$id,
               k = nn$k,
               eps = nn$eps,
               sort = nn$sort)

    # n.neighbors
    n.neighbors <- lengths(nn@id)

    # exp.per.ct.cn
    exp.per.ct.cn <- data.table::rbindlist(lapply(nn1,function(nn)nn@exp.per.ct.cn))

    # exp.per.cn
    exp.per.cn <- data.table::rbindlist(lapply(nn1,function(nn)nn@exp.per.cn))

    ## rcn.count
    rcn.count <- as.matrix(data.table::rbindlist(lapply(nn1,function(fr)as.data.frame(fr@rcn.count))))

    ## rcn.freq
    rcn.freq <- as.matrix(data.table::rbindlist(lapply(nn1,function(fr)as.data.frame(fr@rcn.freq))))

    ## is.selected
    mclustda <- list()
    mclustda$sele <- mclustda$all <- list()

    cn <- new("CellNeighborhood",
              within.rois=within.rois,
              n.cells.selected=n.cells.selected,
              is.selected=is.selected,
              used.cts=used.cts,
              smpls=smpls,
              nn=nn,
              exp.per.ct.cn=exp.per.ct.cn,
              exp.per.cn=exp.per.cn,
              rcn.count=rcn.count,
              rcn.freq=rcn.freq,
              mclustda = mclustda)

    return(cn)
  }
)

#_ -------------------------------

# fun: setDist ----

#' set distance to tumorBorder for NN objects
#'
#' @param x A NN object.
#' @param value A numeric vector specifying the distance to tumor border for each cell.
#'
#' @export
setGeneric("setDist", function(x,...) standardGeneric("setDist"))

#' @rdname setDist
#' @export
setMethod("setDist", "CellNeighborhood",
  function(x,value){
    x@dist2tumorBorder <- value # only for Cycif, and CycifStack
    return(x)
  }
)

#_ -------------------------------
# fun: tcnClust ----

#' Cluster and Sort Recurrent Cell Neighbors (RCN) for CyCIF or CyCIFStack Objects
#'
#' This function clusters and sorts the Recurrent Cell Neighbors (RCN) for CyCIF or CyCIFStack objects.
#' It clusters cells based on their RCN profiles, sorts clusters based on the specified cell type,
#' and optionally extrapolates the clustering to the entire dataset.
#'
#' @param nn An object containing RCN information, typically obtained from 'computeCN'.
#' @param g The number of clusters to create.
#' @param seed The random seed for reproducibility.
#' @param sort.by The cell type to sort clusters by (e.g., "CD8T").
#' @param sort.type The type of data to use for sorting: "freq" (relative frequencies) or "count" (counts).
#' @param sort.smpls The subset of samples to use for sorting: "all" (entire dataset) or "selected" (selected cells).
#' @param data.type The type of data to use for clustering: "ct_exp" (cell type and expression data) or "ct" (cell type data only).
#' @param extrapolate Whether to extrapolate clusters to the entire dataset (TRUE) or not (FALSE).
#' @param mc.cores The number of CPU cores to use for parallel processing.
#'
#' @return An updated 'nn' object with clustering and sorting information.
#'
#' @details
#' The `tcnClust` function uses the provided `nn` object to perform clustering and classification of cells based on their neighborhood relationships. It allows you to specify the number of clusters (`g`), the cell type to sort clusters by (`sort.by`), and other clustering parameters.
#' The clustering process results in the classification of cells into distinct clusters, and the function provides information about these clusters, including the mean frequencies, counts, and more.
#' By specifying different options for `sort.by`, `sort.type`, and `sort.smpls`, you can customize the sorting behavior of clusters based on cell types and data types.
#' Additionally, you can choose to extrapolate clusters to the entire dataset using the `extrapolate` argument, which can be helpful for analyzing the overall dataset.
#'
#' @seealso \code{\link{computeCN}}
#'
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @importFrom dbscan  frNN kNN
#' @importFrom data.table as.data.table rbindlist
#' @importFrom parallel mclapply
#' @export
setGeneric("tcnClust", function(nn,...) standardGeneric("tcnClust"))

#' @rdname tcnClust
#' @export
setMethod("tcnClust","data.frame",
  function(nn,
           g=50,
           seed=123,
           sort.by="CD8T",
           sort.type=c("freq","count"),
           sort.smpls=c("all","selected"),
           data.type=c("ct_exp","ct"),
           extrapolate=FALSE,
           mc.cores=1){
  mclustda <- nn@mclustda

  exps <- as.matrix(nn@exp)[,-1]
  exps.imp <- imputeData(exps)

  this.cts <- cts.in.rcn

  mat.count.all <- nn$rcn.count[,this.cts]
  mat.freq.all <- t(apply(nn$rcn.freq[,this.cts],1,function(x)x/sum(x)))

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
      predict(mc1,newdata=mat.freq.all1[mc.idx==i,nn$cts.in.rcn])$classification
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

  nn$mclustda <- mclustda

  return(nn)
})


#_ -------------------------------

# fun: rcnClust ----

#' Cluster and Sort Recurrent Cell Neighbors (RCN) for CyCIF or CyCIFStack Objects
#'
#' This function clusters and sorts the Recurrent Cell Neighbors (RCN) for CyCIF or CyCIFStack objects.
#' It clusters cells based on their RCN profiles, sorts clusters based on the specified cell type,
#' and optionally extrapolates the clustering to the entire dataset.
#'
#' @param cn An object containing CN information, typically obtained from 'computeCN'.
#' @param g The number of clusters to create.
#' @param seed The random seed for reproducibility.
#' @param sort.by The cell type to sort clusters by (e.g., "CD8T").
#' @param sort.type The type of data to use for sorting: "freq" (relative frequencies) or "count" (counts).
#' @param sort.smpls The subset of samples to use for sorting: "all" (entire dataset) or "selected" (selected cells).
#' @param data.type The type of data to use for clustering: "ct_exp" (cell type and expression data) or "ct" (cell type data only).
#' @param extrapolate Whether to extrapolate clusters to the entire dataset (TRUE) or not (FALSE).
#' @param mc.cores The number of CPU cores to use for parallel processing.
#'
#' @return An updated 'nn' object with clustering and sorting information.
#'
#' @seealso \code{\link{computeCN}}
#'
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @importFrom data.table as.data.table rbindlist
#' @importFrom parallel mclapply
#' @export
setGeneric("rcnClust", function(cn,...) standardGeneric("rcnClust"))

#' @rdname rcnClust
#' @export
setMethod("rcnClust","CellNeighborhood",
          function(cn,
                   g=50,
                   seed=123,
                   sort.by="dist",
                   sort.type=c("freq","count"),
                   sort.smpls=c("all","selected"),
                   data.type=c("ct_exp","ct"),
                   extrapolate=FALSE,
                   mc.cores=1){

            mclustda <- cn@mclustda
            is.selected <- cn@is.selected # n1 = 1325874
            smpls <- cn@smpls # n1

            ## exps
            exps <- as.matrix(cn@exp.per.cn)[,-1] # nrow = n1
            exps.imp <- imputeData(exps)

            ## rcn.freq
            rcn.freq <- cn@rcn.freq[,cts.in.rcn] # nrow = n1
            rf <- t(apply(rcn.freq,1,function(x)x/sum(x)))

            ## dist
            dist <- cn@dist2tumorBorder # length = n1

            ## data type
            if(missing(data.type)){
              data.type="ct_exp"
            }

            if(data.type=="ct"){
              df <- rcn.freq
            }else if(data.type=="ct_exp"){
              ## combine the data
              # df <- cbind(exps.imp,rcn.freq,dist)
              df <- cbind(exps.imp,rcn.freq)
            }

            df.sele <- df[is.selected,]
            norm.df.sele <- scale(df.sele,center=TRUE,scale=TRUE)

            ## clustering & classification
            set.seed(seed)

            cat("Clustering with Mclust ...\n")
            mem.ori <- mclust::Mclust(data=norm.df.sele,G=g,modelNames="EII")$classification

            if(extrapolate){
              cat("Training MclustDA ...\n")
              mc1 <- mclust::MclustDA(data=norm.df.sele,class=mem.ori,
                                      G=g,modelNames="EII",modelType = "EDDA") # 100, EII
              mem.sele <- factor(predict(mc1)$classification)
              g1 <- nlevels(mem.sele)
            }else{
              mem.sele <- mem.ori
              g1 <- nlevels(mem.sele)
            }

            ## mean counts & frequencies
            df.sele1 <- cbind(dist[is.selected],df.sele)
            colnames(df.sele1)[1] <- "dist"
            mean.df.sele <- sapply(seq(g1),function(i)colMeans(df.sele1[mem.sele==i,]))
            o <- order(mean.df.sele["dist",],decreasing=T)

            ## update labels
            mem.sele1 <- factor(as.numeric(factor(mem.sele,levels=levels(mem.sele)[o])))
            mean.df.sele1 <- sapply(seq(g1),function(i)colMeans(df.sele1[mem.sele1==i,]))

            boxplot(df.sele1[,"dist"] ~ mem.sele1,pch=NA)
            boxplot(df.sele1[,"pTBK1"] ~ mem.sele1,pch=NA)

            par(mfrow=c(2,1))
            boxplot(df.sele1[,"cCaspase3"] ~ mem.sele1,pch=NA)
            boxplot(df.sele1[,"CD8T"] ~ mem.sele1,pch=NA)
            # boxplot(df.sele1[,"pTBK1"] ~ mem.sele1,pch=NA)
            par(mfrow=c(2,1))
            boxplot(df.sele1[,"cCaspase3"] ~ mem.sele1,pch=NA)
            boxplot(df.sele1[,"BCLXL"] ~ mem.sele1,pch=NA)

            plot(t(mean.df.sele1)[,c("cCasepase3","BCLXL")])

            colnames(mean.df.sele1) <- seq(g1)

            heatmap3(mean.df.sele1,Rowv=NA,Colv=NA,scale="row",balanceColor = TRUE)
            heatmap3(mean.df.sele1[1:9,],Rowv=NA,Colv=NA,scale="row",balanceColor = TRUE)
            heatmap3(mean.df.sele1[10:18,],Rowv=NA,Colv=NA,scale="none",balanceColor = FALSE)

            if(extrapolate){
              ## extrapolate clusters
              is.available <- !apply(df,1,function(x)any(is.na(x))) # 8587
              df1 <- df[is.available,]

              cat("Applying MclustDA to the entire data ...\n")
              mc.idx <- sort(rep(seq(mc.cores),length=nrow(df1)))

              mem.all <- parallel::mclapply(seq(mc.cores),function(i){
                predict(mc1,newdata=df1[mc.idx==i,cn@cts.in.rcn])$classification
              },mc.cores=mc.cores)
              mem.all <- do.call(c,mem.all)

              ## update labels
              mem.all <- factor(mem.all)
              g1 <- nlevels(mem.all)
              levels(mem.all) <- seq(g1)

              ## mean counts & frequencies
              mean.df.all <- sapply(seq(g1),function(i)colMeans(df1[mem.all==i,]))
              colnames(mean.freq.all) <- seq(g1)
            }

            ## Sort clusters based on frequency of a cell type (CD8T by default)
            if(extrapolate){
              sort.smpls <- "all"
            }else{
              sort.smpls <- "selected"
            }

            if(sort.type == "freq" & sort.smpls == "all"){
              mean.freq <- mean.freq.all
            }else if(sort.type == "freq" & sort.smpls == "selected"){
              mean.freq <- mean.freq.sele
            }

            if(sort.by=="dist"){
              o <- order(mean.freq[sort.by,],decreasing=T)
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
            }

            ##
            mclustda$g <- g1

            cn$mclustda <- mclustda

            return(cn)
          })
#_ -------------------------------
# fun: meanExpRCN ----

#' @title Compute Mean Expression Profiles per RCN Cluster
#'
#' @description
#' This function computes the mean expression profiles for specified cell types
#' or features within Recurrent Cell Neighborhood (RCN) clusters. It takes a data frame,
#' typically containing expression data, and computes the mean expression values
#' for each RCN cluster based on the provided RCN information from 'computeCN'.
#' The function allows you to focus on specific cell types and features and can
#' extrapolate the clustering results to the entire dataset if needed.
#'
#' @param x A data frame containing expression data, typically from a CyCIF or similar dataset.
#' @param nn An object containing RCN information, typically obtained from 'computeCN'.
#' @param cts.in.center A character vector specifying the cell typesaround which RCN was computed (e.g., "Tumor").
#' @param cts.in.rcn A character vector specifying the cell types to include in the RCN analysis.
#' @param per.ct A logical value indicating whether to compute mean expression profiles per RCN cluster (TRUE) or for the entire dataset (FALSE).
#' @param extrapolate A logical value indicating whether to extrapolate clustering results to the entire dataset (TRUE) or not (FALSE).
#'
#' @return A list of data frames containing mean expression profiles for specified cell types or features within RCN clusters.
#'
#' @details
#' The 'meanExpRCN' function calculates the mean expression profiles for the specified cell types or features within RCN clusters. The function works as follows:
#' - It takes a data frame 'x' containing expression data and an object 'nn' containing RCN information obtained from the 'computeCN' function.
#' - You can specify the 'cts.in.center' argument to select specific cell types to focus on during the analysis.
#' - The 'cts.in.rcn' argument allows you to specify the cell types to include in the RCN analysis.
#' - If 'per.ct' is set to TRUE, the function computes mean expression profiles per RCN cluster; otherwise, it computes mean expression profiles for the entire dataset.
#' - The 'extrapolate' argument determines whether to extrapolate clustering results to the entire dataset based on RCN information.
#' The function returns a list of data frames containing mean expression profiles for the specified cell types or features within RCN clusters.
#'
#' @seealso \code{\link{computeCN}}, \code{\link{tcnClust}}
#'
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist
#' @importFrom tibble rowid_to_column
#' @importFrom dplyr %>% arrange left_join summarize_at mutate_at mutate group_by summarize
#' @importFrom tidyr spread
#' @importFrom ggplot2 ggplot aes geom_point geom_line sym syms
#'
#' @export
setGeneric("meanExpRCN", function(x,...) standardGeneric("meanExpRCN"))

#' @rdname meanExpRCN
#' @export
setMethod("meanExpRCN","data.frame",
          function(x,
                   nn,
                   cts.in.center="Tumor",
                   cts.in.rcn=levels(cell_types(x)$cell_types)[1:10],
                   per.ct=TRUE,
                   extrapolate=TRUE){

    g <-   nn$mclustda$g

    lin.abs <- names(x@cell_types$default@cell_lineage_def[-(1:2)])
    cst.abs <- names(x@cell_types$default@cell_state_def)
    all.abs <- unique(c(lin.abs,cst.abs))

    if(extrapolate){
      mclustda <- nn$mclustda$all
    }else{
      mclustda <- nn$mclustda$sele
    }

    # clusts <- factor(mclustda$mem)
    clusts <- factor(paste0("RCN",mclustda$mem),labels=paste0("RCN",seq(g)))
    df1 <- cell_types(x) %>%
      filter(nn$within.rois) %>%
      tibble::rowid_to_column("idx") %>% ## idx is numbered within 'within.rois'
      dplyr::left_join(
        data.frame(tcn = clusts) %>%
          mutate(idx=which(mclustda$is.used)),by="idx") %>%
      dplyr::arrange(idx) %>%
      dplyr::left_join(exprs(x,type="log") %>%
                  filter(nn$within.rois) %>%
                  tibble::rowid_to_column("idx"),by="idx") %>%
      rename(idx.all="idx")## 1325874, within ROIs

    idx.nonna.tum <-!is.na(df1$tcn) & df1$cell_types %in% cts.in.center # 471910 / 1325874
    df.tum <- df1 %>%
      filter(idx.nonna.tum) %>%
      tibble::rowid_to_column("idx.tum") # 471910

    ## lst.nn: convert ids in each sample to ids in all samples
    nn.ids <- nn$nn$id[!is.na(df1$tcn)] # 1325311: 563 don't have proper RCNs.
    nn.tum.ids <- nn$nn$id[idx.nonna.tum] # 471910

    ## all unique ids per tcn cluster
    if(per.ct){
      lst.mean.exp <- tapply(nn.tum.ids,df.tum$tcn,function(this.ids){
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
      lst.mean.exp <- tapply(nn.tum.ids,df.tum$tcn,function(this.ids){
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
)
