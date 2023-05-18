#_ -------------------------------------------------------

# utils Cycif ----

#' @export
setMethod("names", "Cycif", function(x) x@name)

#' @export
setGeneric("abs_list", function(x) standardGeneric("abs_list"))
setMethod("abs_list", "Cycif", function(x) x@abs_list)

#' @export
setGeneric("dna", function(x) standardGeneric("dna"))
setMethod("dna", "Cycif", function(x) x@dna)

#' @export
setGeneric("nCells", function(x) standardGeneric("nCells"))
setMethod("nCells", "Cycif", function(x) x@n_cells)

#' @export
setGeneric("xys", function(x) standardGeneric("xys"))
setMethod("xys", "Cycif", function(x) x@xy_coords)

#' @export
setGeneric("segProp", function(x) standardGeneric("segProp"))
setMethod("segProp", "Cycif", function(x) x@segment_property)

#' @export
setGeneric("dna_thres", function(x) standardGeneric("dna_thres"))
setMethod("dna_thres", "Cycif", function(x) x@dna_thres)

#' @export
setGeneric("used_cells", function(x) standardGeneric("used_cells"))
setMethod("used_cells", "Cycif", function(x) x@used_cells)

#' @export
setGeneric("within_rois", function(x) standardGeneric("within_rois"))
setMethod("within_rois", "Cycif", function(x) x@within_rois)

#' #' #' Replace the values of within_rois for a Cycif object
#' #' #' @param object A cycif object.
#' #' #' @param value logical. If TRUE, the area within ROIs is concerned for the downstream
#' #' #' @export
#' #' within_rois <- function(object, value) {
#' #'   if(!is(object,"Cycif")){
#' #'     stop("this function applies only to a Cycif object")
#' #'   }
#' #'   object@within_rois <- value
#' #'   return(object)
#' #' }
#'
#' #' Replace the values of within_rois for a Cycif object
#' #'
#' #' @param obj A \code{\link{Cycif}} object.
#' #' @param value A named list of data.frames, with names corresponding to
#' #'   the names of \code{object}'s channels.
#' #'
#' #' @return An updated \code{\link{Cycif}} object.
#' #'
#' #' @usage
#' #' within_rois(x) <- value
#' #' @export
#' setReplaceMethod("within_rois",
#'   signature(obj = "Cycif", value = "logical"),
#'     function(obj, value) {
#'       obj@name <- value
#'       return(obj)
#'   })

#' @export
setGeneric("nCycles", function(x) standardGeneric("nCycles"))

#' @rdname nCycles
#' @export
setMethod("nCycles", "Cycif", function(x) x@n_cycles)

#' Set nCycles - note this should run before cell type calling
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export
setGeneric("nCycles<-", function(x,...,value) standardGeneric("nCycles<-"))
setMethod("nCycles<-", "Cycif", function(x,value){

  if(!is.numeric(value)){
    stop("value should be numeric")
  }
  value <- as.integer(value)

  current.nc <- nCycles(x)
  if(any(value > current.nc)){
    stop("New max n_cycles should be less than or equal to the current max n_cycles.")
  }

  ## n_cycles, abs_list, raw, log_normalized, logTh_normalized, dna,
  ## dna_thres
  x@n_cycles <- value

  x@abs_list <- abs_list(x) %>% filter(cycle <= value)
  this.abs <- as.character(abs_list(x)$ab)

  x@raw <- x@raw[this.abs]
  if(all(this.abs %in% names(x@log_normalized))){
    x@log_normalized <- x@log_normalized[this.abs]
  }
  if(all(this.abs %in% names(x@logTh_normalized))){
    x@logTh_normalized <- x@logTh_normalized[this.abs]
  }

  x@dna <- x@dna[paste0("DNA",seq(value))]

  if(nrow(x@used_cells)>0){
    x@used_cells <- x@used_cells[,seq(value)]
  }

  if(nrow(x@dna_thres)>0){
    x@dna_thres <- x@dna_thres[paste0("DNA",seq(value)),]
  }

  if(length(x@rois)>0){
    ncys <- sapply(x@rois,function(y)y$cycle)
    x@rois <- x@rois[ncys <= value]
  }
  return(x)
})

#' Show gates in a Cycif/CycifStack/CellTypeCycif/CellTypeCicyfStack object
#'
#' @rdname gates
#' @export
setGeneric("gates", function(x) standardGeneric("gates"))

#' @rdname gates
#' @export
setMethod("gates", "Cycif", function(x)x@cell_type@gates)

#_ -------------------------------------------------------

# fun: parse_ft path ----
parse_ft <- function(path){
  if(!file.exists(path)){
    stop("path for feature table doesn't exist. Make sure to provide the absolute path to the csv file")
  }

  ft_filename <- sub(".+/","",path)
  root_dir <- sub(paste0("/",ft_filename),"",path)

  ftype0 <- grepl("\\.csv$",ft_filename)
  ftype1 <- grepl("^unmicst-",ft_filename)
  ftype2 <- grepl("--unmicst",ft_filename)

  if(!ftype0){
    stop("feature table file should be in csv format")
  }else if(!ftype1 && !ftype2){
    stop("feature table file format is incorrect")
  }else{
    f1 <- strsplit(ft_filename,"[-_]")[[1]]
    if(ftype1){
      smpl <- f1[2]
      mask_type <- f1[3]
    }else if(ftype2){
      smpl <- f1[1]
      mask_type <- f1[4]
    }

    fp <- new("file_paths",
              sample_name = smpl,
              root_dir = root_dir,
              mask_type=mask_type,
              registration="",
              segmentation="",
              quantification=path
    )

    return(fp)
  }
}


# fun: find_mcmicro_output_path path ----
find_mcmicro_output_path <- function(smpl.path,mask_type=c("cellRing","cell")){
  smpl <- sub("^.+/","",smpl.path)

  ## registration
  reg.file <- dir(file.path(smpl.path,"registration"),pattern="\\.ome.tif")
  reg.path <- file.path(smpl.path,"registration",reg.file)

  ## quantification
  qua.files <- dir(file.path(smpl.path,"quantification"),pattern="\\.csv$")
  msks <- sub(".+.+_([^_]+)\\.csv","\\1",qua.files)
  idx <- which.min(match(msks,mask_type))
  mask_type <- msks[idx]
  qua.file <- qua.files[idx]
  qua.path <- file.path(smpl.path,"quantification",qua.file)

  ## segmentation
  seg.files <- dir(file.path(smpl.path,"segmentation",paste0("unmicst-",smpl)))
  seg.file <- seg.files[sub(".ome.tif","",seg.files) %in% mask_type]
  seg.path <- file.path(smpl.path,"segmentation",paste0("unmicst-",smpl),seg.file)

  ## return a file_paths obj
  new("file_paths",
      sample_name = smpl,
      root_dir = smpl.path,
      mask_type=mask_type,
      registration=reg.path,
      quantification=qua.path,
      segmentation=seg.path
  )
}

#_ -------------------------------------------------------

# fun: constructor Cycif ----
#' Instantiate and show a Cycif object
#'
#' @param ft_filename quantification file (.csv)
#' @param path path to the dir containing ft_filename
#' @param mask_type a common mask_type of each file name or channel. If applicable, the mask_type is
#'   removed from the channel names.
#' @param object a Cycif object
#' @rdname Cycif
#'
#' @export
Cycif <- function(
    filename,
    path="",
    mask_type=c("cellRing","cell"),
    mcmicro=FALSE,
    use_scimap=FALSE){
  filename <- file.path(path,filename)
  if(missing(filename)){
    stop("'filename' is missing.")
  }
  if(!file.exists(file.path(path,ft_filename))){
    stop("the provided path doesn't exist")
  }

  ## file_paths, mask_type, name
  if(mcmicro){
    smpl.path <- filename
    fp <- find_mcmicro_output_path(smpl.path=smpl.path,mask_type=mask_type)
  }else{
    mt <- sub(".+_(.+)\\.csv","\\1",filename)
    if(!any(mt %in% mask_type)){
      stop("The feature table path doesn't match the 'mask_type'.")
    }
    ft_path <- filename
    fp <- parse_ft(ft_path)
  }

  ft_path <- fp@quantification
  mask_type <- fp@mask_type

  name <- fp@sample_name

  ## n_cycles, n_cells, abs_list
  txt <- utils::read.csv(ft_path)
  n.cells <- nrow(txt)

  cn <- names(txt)
  is.ch <- !cn %in% c("CellID","X_centroid","Y_centroid","Area","MajorAxisLength","MinorAxisLength",
                      "Eccentricity","Solidity","Extent","Orientation")

  ch.list <- sub(mask_type,"",cn[is.ch])
  if(any(grepl("centroid",ch.list))){
    ch.list <- ch.list[!grepl("centroid",ch.list)]
  }

  i.dna <- grep("^DNA",ch.list)
  ab.list <- ch.list[-i.dna]

  n.abs <- length(ab.list)

  n.cycles <- (n.abs/3)

  abs_list <- data.frame(ab=ab.list,
                         cycle=rep(seq(n.cycles),each=3),
                         channel=rep(1:3,times=n.cycles),
                         stringsAsFactors=FALSE)

  dna.list <- ch.list[i.dna]

  ## raw and dna
  sub.txt <- txt[is.ch]
  names(sub.txt) <- ch.list

  raw <- sub.txt[ab.list]
  dna <- sub.txt[dna.list]
  if(dna.list[1]=="DNA0"){
    dna.list <- paste0("DNA",seq(dna.list))
    names(dna) <- dna.list
  }

  ## xy_coords
  xy_coords <- txt[c("X_centroid","Y_centroid")]

  ## segment_property
  segment_property <- txt[c("Area","MajorAxisLength",
                            "MinorAxisLength", "Eccentricity","Solidity",
                            "Extent","Orientation")]

  obj <- new("Cycif",
      name = name,
      file_paths = fp,
      mask_type = mask_type,
      n_cycles = n.cycles,
      n_cells = n.cells,
      abs_list = abs_list,
      raw = raw,
      dna = dna,
      xy_coords = xy_coords,
      segment_property = segment_property)

  ## anndata
  if(use_scimap){
    if(!reticulate::py_module_available("scimap")){
      stop("'scimap' in python is not available.\n",
           "Install it by pip install PyQt5 scimap, and\n",
           "specify the python path in RETICULATE_PYTHON in .Renviron.\n",call.=FALSE)
    }else{
      adata = sm$pp$mcmicro_to_scimap(ft_path)
      obj@adata = adata
    }
  }
  return(obj)
}


# fun: setValidity Cycif ----
setValidity("Cycif", function(object) {
  if (!all(sapply(object@rois, inherits, "roi"))) {
    "All elements of roi_list must be of class 'roi'"
  }
})

# fun: show Cycif ----
#' @rdname Cycif
#' @export
setMethod("show", "Cycif", function(object) {
  m <- matrix(as.character(abs_list(object)$ab),nrow=3)
  rownames(m) <- c("R","G","B")
  colnames(m) <- paste0("Cycle",seq(nCycles(object)))
  m <- data.frame(m)

  nd <- nrow(dna_thres(object)) > 0
  nln <- nrow(object@log_normalized) > 0
  nltn <- nrow(object@logTh_normalized) > 0
  nct <- length(object@cell_type@cell_types) > 0
  # cat(m)
  cat("[",is(object)[[1]], " object]\n\n",
      "Name: ", object@name, "\n",
      "nCells:  ", object@n_cells, "\n",
      "nCycles: ", object@n_cycles, "\n\n",
      "DNA filtered: ", nd, "\n",
      "Log-normalized: ", nln, "\n",
      "LogTh-normalized: ", nltn, "\n",
      "CellType-idenitified: ", nct, "\n\n",
      "Abs: \n",sep="")
  print(m)
})



#_ -------------------------------------------------------

# fun: list2Cycif list ----
#' Converting a list to a Cycif object
#'
#' @param x A list object
#' @return A Cycif object
#'
#' @export
setGeneric("list2Cycif", function(x) standardGeneric("list2Cycif"))

#' @export
setMethod("list2Cycif", "list",function(x){
  slots.in.Cycif <- slotNames(getClassDef("Cycif"))
  stopifnot(all(names(x) %in% slots.in.Cycif))
  n.samples <- sum(sapply(x,length))
  xs <- do.call(c,lapply(x,function(cs){
    if(is(cs,"Cycif")){
      out <- cs@samples
    }else if(is(cs,"Cycif")){
      out <- cs
    }
  }))

  mask_type <- xs[[1]]@mask_type
  nms <- 	names(xs) <- sapply(xs,names)
  n_cells <- 	sapply(xs,nCells)
  n_cycles <- sapply(xs,nCycles)
  max_cycles <- max(n_cycles)

  idx <- which(n_cycles == max_cycles)[1]
  abs_list <- abs_list(xs[[idx]])
  new("Cycif",
     n_samples = n.samples,
     names = nms,
     abs_list = abs_list,
     n_cycles = n_cycles,
     max_cycles = max_cycles,
     n_cells = n_cells,
     samples = xs,
     mask_type = mask_type
  )
})

#_ -------------------------------------------------------

# fun: cumUsedCells Cycif ----
#' Show available cells at each cycle
#'
#' @param x A Cycif object
#' @param cumulative logical. If TRUE, available cells through cycles 1 to n are computed.
#'   If FALSE, available cells at each cycle are computed.
#' @param ratio Logical. If TRUE, show ratio of available and unavailable cells
#'  at each cycle. If FALSE, actual actual number is provided.
#'
#' @usage
#' cumUsedCells(x,within.rois=TRUE)
#'
#' @return cumUsedCells() returns a matrix where rows correspond to cells
#'   and columns correspond to cycles. statUsedCells() returns a table sumarizing
#'   the numbers of available and unavailable cells at each cycle.
#'
#' @export
setGeneric("cumUsedCells", function(x,...) standardGeneric("cumUsedCells"))
setMethod("cumUsedCells", "Cycif",
          function(x,within.rois=TRUE){
            u <- x@used_cells
            w <- x@within_rois

            u <- sapply(seq(ncol(u)),function(i){
              id <- rowSums(u[,seq(i),drop=F]==1)==i
              return(id)
            })

            if(within.rois & length(w)>0){
              u <- u & w
            }

            return(u)
          }
)

# fun: statUsedCells Cycif ----

#' @export
setGeneric("statUsedCells", function(x,...) standardGeneric("statUsedCells"))

#' Summarize available cells at each cycle
#' @rdname statUsedCells
#'
#' @param cumulative logical. If TRUE, return the number of cells available from cycle 1 to N. FALSE, return the number of cells in each cycle. Default is TRUE.
#' @param ratio logical. If TRUE, return the ratio of available cells over all the cells. Default is TRUE.
#' @param within.rois logical. If TRUE, compute available cells within pre-specified ROIs.
#'
#' @usage
#' statUsedCells(x, cumulative=TRUE, ratio=TRUE, within.rois=TRUE)
#'
#' @order 2
#' @export
setMethod("statUsedCells", "Cycif",
          function(x,cumulative=TRUE,ratio=TRUE,within.rois=TRUE){
            stopifnot(nrow(x@used_cells)>0)
            if(cumulative){
              mat <- cumUsedCells(x,within.rois=within.rois)
            } else{
              mat <- x@used_cells == 1 # 0 - dropped, 1 - alive, 2 - bunched
            }
            tab <- apply(mat,2,function(l){
              f <- factor(l,levels=c("TRUE","FALSE"))
              table(f)
            })
            if(ratio){
              nc.true <- tab["TRUE",]
              tab <- nc.true/nc.true[1]
            }
            return(tab)
          }
)

#_ -------------------------------------------------------
