CellTypeGraph <- function(cy,ctype,plot=F){
  require(igraph)
  ctgraph <- ctype[c("Parent","Child")]  
  g <- graph_from_data_frame(ctgraph)
  l <- layout_as_tree(g)
  levs <- l[,2]
  names(levs) <- names(V(g))
  ctlevs <- rev(tapply(names(levs),levs,identity))
  if(plot){
    plot(g,layout=l)
  }
  return(ctlevs)
}


