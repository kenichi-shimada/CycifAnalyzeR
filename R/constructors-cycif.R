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
    f1 <- sub("\\.csv","",strsplit(ft_filename,"[-_]")[[1]])
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

#' @title Create a Cycif object
#'
#' @description
#' The `Cycif` function creates a Cycif object from a given feature table file.
#'
#' @param ft_filename Character string, the name of the feature table file (a comma-separated file).
#' @param path Character string, the path to the feature table file. Default is ".".
#' @param mask_type Character vector specifying the type of mask used, options are "cellRing" or "cell".
#' @param mcmicro Logical, indicating whether the feature table comes from MCMICRO.
#' @param use_scimap Logical, indicating whether to use scimap python package for additional data processing. The default is FALSE (Note: scimap is used through reticulate but it's use is not fully implemented yet).
#'
#' @return
#' A Cycif object containing information extracted from the feature table.
#'
#' @details
#' The `Cycif` function reads a feature table file and extracts information such as the number of cycles,
#' the number of cells, the list of antibodies, raw protein and DNA intensities, xy coordinates, and segment properties.
#' It then creates a Cycif object with the extracted information.
#'
#' @seealso
#' \code{\link{Cycif}}, \code{\link{CycifStack}}, \code{\link{CellTypes}}
#'
#' @export
Cycif <- function(
    ft_filename,
    path=".",
    mask_type=c("cellRing","cell"),
    mcmicro=FALSE,
    use_scimap=FALSE){
  filename <- file.path(path,ft_filename)
  if(missing(filename)){
    stop("'filename' is missing.")
  }
  if(!file.exists(filename)){
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

  ## cycles
  cycles <- rep(1,length(ch.list))
  cycles[i.dna] <- seq(i.dna)
  j <- 1
  for(i in seq(cycles)){
    if(cycles[i] == j + 1){
      j <- j+1
    }
    cycles[i] <- j
  }

  n.cycles <- max(cycles)

  ## abs_list
  ab.list <- ch.list[-i.dna]
  cycles <- cycles[-i.dna]

  abs_list <- data.frame(ab=ab.list,
                         cycle=cycles,
                         stringsAsFactors=FALSE)

  ## dna.list
  dna.list <- ch.list[i.dna]

  ## raw and dna
  sub.txt <- txt[is.ch]
  names(sub.txt) <- ch.list

  raw <- sub.txt[ab.list]
  dna <- sub.txt[dna.list]
  names(dna) <- paste0("DNA",seq(dna.list))

  ## xy_coords
  xy_coords <- txt[c("X_centroid","Y_centroid")]

  ## segment_property
  segment_property <- txt[c("Area","MajorAxisLength",
                            "MinorAxisLength", "Eccentricity","Solidity",
                            "Extent","Orientation")]

  ## adata
  obj <- new("Cycif",
      name = name,
      file_paths = fp,
      mask_type = mask_type,
      n_cycles = n.cycles,
      n_cells = n.cells,
      abs_list = abs_list,
      raw = raw,
      dna = dna,
      cell_types=list(),
      xy_coords = xy_coords,
      segment_property = segment_property)

  use_scimap=FALSE
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
# setValidity("Cycif", function(object) {
#   if (!is(object, "Cycif")) {
#     stop("Invalid object: 'object' is not of class 'Cycif'.")
#   }
#   if (!is.null(object@adata)) {
#     if (!inherits(object@adata, "anndata._core.anndata.AnnData")) {
#       stop("Invalid object: 'adata' slot is not of the expected class.")
#     }
#   }
#   if (!all(sapply(object@rois, inherits, "roi"))) {
#     "All elements of roi_list must be of class 'roi'"
#   }
# })
