CellTypeGraph <- function(ctype,plot=F){
  ctgraph <- ctype[c("Parent","Child")]
  g <- igraph::graph_from_data_frame(ctgraph)
  l <- igraph::layout_as_tree(g)
  levs <- l[,2]
  names(levs) <- names(igraph::V(g))
  ctlevs <- rev(tapply(names(levs),levs,identity))
  if(plot){
    plot(g,layout=l)
  }
  return(ctlevs)
}


