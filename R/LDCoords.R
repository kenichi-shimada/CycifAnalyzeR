LDCoords <- function(type,smpls,used.abs,used.cts,
                     n_cells_per_smpl,n_cells_total,
                     ld_coords,is_used){
  if(!type %in% c("PCA","tSNE","UMAP")){
    stop("type should be PCA, tSNE, or UMAP.")
  }
  new("LDCoords",
      type = type,

      smpls = smpls,
      used.abs = used.abs,
      used.cts = used.cts,

      n_cells_per_smpl = n_cells_per_smpl,
      n_cells_total = n_cells_total,

      ld_coords = ld_coords,
      is_used = is_used)
}

#' @rdname LDCoords
#' @export
setMethod("show", "LDCoords", function(object) {
  type <- object@type

  n.smpls <- length(object@smpls)
  used.abs <- object@used.abs
  n.abs <- length(used.abs)
  abs <- paste(used.abs,collapse=",")
  used.cts <- object@used.cts
  n.cts <- length(object@used.cts)
  cts <- paste(used.cts,collapse=",")

  n_cells_smpl <- object@n_cells_per_smpl
  n_cells_total <- object@n_cells_total

  cat("[",is(object)[[1]], " object]\n\n",
      "Type: ", object@type, "\n\n",
      "abs (",n.abs,") :",abs,"\n",
      "cts (",n.cts,") :",cts,"\n\n",
      "# cells:  ", object@n_cells_per_smpl, "\n",
      "# cycles: ", object@n_cells_total, "\n",
      sep="")

})
