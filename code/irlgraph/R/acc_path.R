#'@name acc_path
#'@title Generate a shortest path raster from a specific graph node or set of graph coordinates
#'@export
#'@importFrom igraph shortest.paths
#'@param graph igraph object
#'@param snode numeric node number
#'@param scoord numeric coordinates of length 2
#'@param costsurf RasterLayer cost surface

acc_path <- function(graph, costsurf, snode = NULL, scoord = NULL){
  
  if(is.null(snode) & all(is.null(scoord))){
    stop("Must supply either a starting node or coordinates")
  }
  
  if(is.null(snode)){
    snode <- raster::cellFromXY(costsurf, scoord)
  }
  
  #ADD A CHECK TO MAKE SURE SNODE EXISTS IN GRID!!!
  #snode %in% graph$cells 
  
  spaths <- igraph::shortest.paths(graph$graph, v = snode)
  
  result <- as(costsurf, "RasterLayer")
  result[] <- NA
  result[] <- spaths
  result
}
