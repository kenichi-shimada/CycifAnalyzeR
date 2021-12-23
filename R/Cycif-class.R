#' @include CellType-class.R
#' @export
setClass("Cycif",
         slots = c(
           ## using qts
           name = "character",
           suffix = "character",

           abs_list = "data.frame",
           used_abs = "vector",

           raw = "data.frame",
           normalized = "data.frame",
           normalize.method = "character",
           dna = "data.frame",

           n_cycles = "numeric",
           n_cells = "numeric",

           xy_coords = "data.frame",
           segment_property = "data.frame",

           cell_type = "CellTypeCycif",
           threshold = "numeric", ## this is going to be discontinued

           used_cells = "matrix", # matrix
           dna_thres = "data.frame",
           area_filter= "vector",

           ld_coords = "data.frame",
           clusters = "numeric", # numeric

           tif_files = "data.frame", # full path
           dim = "data.frame", # size of the tif_file
           rect = "list",

           segmentation_masks = "character",

           commands = "list"
         )
)

#' @export
Cycif <- function(filename,suffix="_cellMask") {
  stopifnot(file.exists(filename))

  name <- sub("unmicst-(.+)\\.csv","\\1",filename)

  txt <- read.csv(filename)
  n.cells <- nrow(txt)

  cn <- names(txt)
  is.ch <- grepl(suffix,cn)
  ch.list <- sub(suffix,"",cn[is.ch])
  i.dna <- grep("^DNA",ch.list)
  ab.list <- ch.list[-i.dna]

  n.abs <- length(ab.list)
  n.cycles <- (n.abs/3)
  abs_list <- data.frame(ab=ab.list,
                         cycle=rep(seq(n.cycles)-1,each=3),
                         channel=rep(1:3,times=n.cycles),
                         stringsAsFactors=TRUE)

  abs_params <- data.frame(ab=ab.list,
                           low=rep(-Inf,n.abs),
                           high=rep(Inf,n.abs),
                           use=rep(TRUE,n.abs))

  dna.list <- ch.list[i.dna]

  ## raw
  sub.txt <- txt[is.ch]
  names(sub.txt) <- ch.list

  raw <- sub.txt[ab.list]
  dna <- sub.txt[dna.list]

  xy_coords <- txt[c("X_centroid","Y_centroid")]
  # xy_coords$Y_centroid <- max(xy_coords$Y_centroid) - xy_coords$Y_centroid

  segment_property <- txt[c("Area","MajorAxisLength",
                            "MinorAxisLength", "Eccentricity","Solidity",
                            "Extent","Orientation")]

  new("Cycif",
      name = name,
      suffix = suffix,
      abs_list = abs_list,
      raw = raw,
      dna = dna,
      n_cycles = n.cycles,
      n_cells = n.cells,
      xy_coords = xy_coords,
      segment_property = segment_property)
}

#' @export
setMethod("show", "Cycif", function(object) {
  m <- matrix(abs_list(object)$ab,nrow=3)
  rownames(m) <- c("R","G","B")
  colnames(m) <- paste0("Cycle",seq(nCycles(object))-1)
  m <- data.frame(m)
  # cat(m)
  cat(is(object)[[1]], "\n",
      "Name: ", object@name, "\n",
      "nCells:  ", object@n_cells, "\n",
      "nCycles: ", object@n_cycles, "\n",
      "Abs: \n",sep="")
  print(m)
})
