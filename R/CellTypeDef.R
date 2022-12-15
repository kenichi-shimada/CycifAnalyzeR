CellTypeDef <- function(x,lineage_df,state_df){
  if(class(x)=="Cyclf"){
    used.abs <- as.character(abs_list(x)$ab)
  }else if(class(x)=="CycifStack"){
    used.abs <- as.character(uniq_abs(x)$ab)
  }else{
    stop("1st argument should be either a Cycif or CycifStack object")
  }
  if(missing(lineage_df) || missing(state_df)){
    stop("both lineage and state definitions should be provided")
  }
  if(class(lineage_df) != "data.frame"){
    stop("cell lineage definition should be a data.frame")
  }
  if(class(state_df) != "data.frame"){
    stop("cell state definition should be a data.frame")
  }
  mks <- unique(c(colnames(lineage_df)[-(1:2)],colnames(state_df)))
  if(!all(mks %in% used.abs)){
    unknown <- mks[!mks %in% used.abs]
    stop(paste0("some markers in the defs not found: ",paste(unknown,collapse=",")))
  }
  new("CellTypeDef",
      cell_lineage_df = lineage_df,
      cell_state_mks = state_df,
      markers = mks
  )
}
