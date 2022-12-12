#' Instantiate and show a Cycif object
#'
#' @include Cycif-class.R
#' @param filename quantification file (.csv)
#' @param suffix a common suffix of each file name or channel. If applicable, the suffix is
#'   removed from the channel names.
#' @param object a Cycif object
#' @rdname Cycif
#' @export
Cycif <- function(filename,path,suffix="_cellRing") {
  stopifnot(file.exists(file.path(path,filename)))

  if(grepl(suffix,filename)){
    name <- sub(paste0("unmicst-(.+)",suffix,"\\.csv"),"\\1",filename)
  }else{
    name <- sub(paste0("unmicst-(.+)\\.csv"),"\\1",filename)
  }

  names(name) <- NULL

  if(!missing(path)){
    filename <- file.path(path,filename)
  }

  txt <- read.csv(filename)
  n.cells <- nrow(txt)

  cn <- names(txt)
  #is.ch <- grepl(suffix,cn)
  is.ch <- !cn %in% c("CellID","X_centroid","Y_centroid","Area","MajorAxisLength","MinorAxisLength",
    "Eccentricity","Solidity","Extent","Orientation")
  ch.list <- sub(suffix,"",cn[is.ch])
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

#' @rdname Cycif
#' @export
setMethod("show", "Cycif", function(object) {
  m <- matrix(abs_list(object)$ab,nrow=3)
  rownames(m) <- c("R","G","B")
  colnames(m) <- paste0("Cycle",seq(nCycles(object)))
  m <- data.frame(m)
  # cat(m)
  cat(is(object)[[1]], "\n",
      "Name: ", object@name, "\n",
      "nCells:  ", object@n_cells, "\n",
      "nCycles: ", object@n_cycles, "\n",
      "Abs: \n",sep="")
  print(m)
})
