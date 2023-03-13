#' @rdname LDCoords
#' @export
LDCoords <- function(ld_type,norm_type,smpls,used.abs,used.cts,
                     n_cells_per_smpl,n_cells_total,
                     ld_coords,is_used,
                     cts_params,ld_params,
                     ld_call,clust_call){
  if(!ld_type %in% c("PCA","tSNE","UMAP")){
    stop("ld_type should be PCA, tSNE, or UMAP.")
  }
  new("LDCoords",
      ld_type = ld_type,
      norm_type = norm_type,

      smpls = smpls,
      used.abs = used.abs,
      used.cts = used.cts,

      n_cells_per_smpl = n_cells_per_smpl,
      n_cells_total = n_cells_total,

      ld_coords = ld_coords,
      is_used = is_used,
      cts_params = cts_params,
      ld_params = ld_params,
      ld_call=ld_call,
      clust_call=clust_call
      )
}

#' @rdname LDCoords
#' @export
setMethod("show", "LDCoords", function(object) {
  n.smpls <- length(object@smpls)
  used.abs <- object@used.abs
  n.abs <- length(used.abs)
  abs <- paste(used.abs,collapse=",")
  used.cts <- object@used.cts
  n.cts <- length(object@used.cts)
  cts <- paste(used.cts,collapse=",")

  n_cells_per_smpl <- object@n_cells_per_smpl
  n_cells_total <- object@n_cells_total

  cat("[",is(object)[[1]], " object]\n\n",
      "Type: ", object@ld_type, "\n",
      "Normalization: ",object@norm_type, "\n\n",
      "cts (",n.cts,") : ",cts,"\n",
      "abs (",n.abs,") : ",abs,"\n\n",
      "# samples : ", n.smpls,"\n",
      "# cells per smpl :\t", n_cells_per_smpl, "\n",
      "# cells in total :\t", sum(n_cells_total), "\n\n",sep="")
  if(length(n_cells_total)>1){
    print(n_cells_total)
  }
})
