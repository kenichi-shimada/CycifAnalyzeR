#' Instanciate and show a CellTypeDef object
#'
#' @param filename Excel sheet that contains two spreadsheets, named 'cell lineage' and 'cell state', respectively.
#' @param object A CellTypeDef object
#'
#' @export
CellTypeDef <- function(filename) {
  stopifnot(file.exists(filename))

  sheet.names <- openxlsx::getSheetNames(filename)

  if(!all(c("cell lineage","cell state") %in% sheet.names)){
    stop("Cell type and cell state should be defined in two spreadsheets named 'cell type' and 'cell state', respectively.")
  }

  clineage <-  openxlsx::readWorkbook(filename,sheet="cell lineage",colNames=TRUE,rowNames=TRUE)
  cstate <-  openxlsx::readWorkbook(filename,sheet="cell state",colNames=TRUE,rowNames=TRUE)

  r1 <- rownames(clineage)
  r2 <- rownames(cstate)
  c1 <- colnames(clineage)
  c2 <- colnames(cstate)

  if(!identical(r1,r2)){
    stop("row.names in 'cell lineage' and 'cell state' should be identical.")
  }else{
    lineages <- r1
    lineages <- gsub(", ","_",lineages)
    rownames(clineage) <- rownames(cstate) <- lineages
  }

  c12 <- intersect(c1,c2)
  if(length(c12)>0){
    stop("Abs defined as both lineage markers and state markers: ",paste(c12,collapse=", "))
  }else{
    markers <- list(lineage=c1,state=c2)
  }

  new("CellTypeDef",
      cell_lineage_df = clineage,
      cell_state_df = cstate,
      lineages = lineages,
      markers = markers
  )
}

#' @rdname CellTypeDef
#' @export
setMethod("show", "CellTypeDef", function(object) {
  lins <- lineages(object)
  n.lins <- length(lins)

  mks <- markers(object)
  n.mks <- length(unlist(mks))
  n.lm <- length(mks$lineage)
  n.sm <- length(mks$state)

  used.abs <- used_abs(object)
  n.uab <- length(used.abs)

  abs.list <- uniq_abs(object)

  ## first three examples
  lin.txt <- ifelse(n.lins >0,paste0(" (",paste(lins[seq(min(n.lins,3))],collapse=", "),",...)"),"")
  lm.txt <- ifelse(n.lm >0,paste0(" (",paste(mks$lineage[seq(min(n.lm,3))],collapse=", "),",...)"),"")
  sm.txt <- ifelse(n.sm >0,paste0(" (",paste(mks$state[seq(min(n.sm,3))],collapse=", "),",...)"),"")
  uab.txt <- ifelse(n.uab >0,paste0(" (",paste(used.abs[seq(min(n.uab,3))],collapse=", "),",...)"),"")

  cat(is(object)[[1]], "\n",
      " lineages:  ", n.lins, lin.txt,"\n",
      " markers:   ", n.mks, "\n",
      "   lineage: ", n.lm, lm.txt,"\n",
      "   state:   ", n.sm, sm.txt,"\n",
      " used_abs:  ", n.uab, uab.txt,"\n",
      " uniq_abs:  ",
      sep="")
  if(nrow(abs.list)>0){
    cat("\n")
    print(abs.list)
  }else{
    cat("not set\n")
  }
})
