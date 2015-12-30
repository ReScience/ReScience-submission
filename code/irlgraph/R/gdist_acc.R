#'@name gdist_acc
#'@title Generate an accumulated cost surface using a regular landscape graph and the gdistance transition class
#'@description The igraph package is used to contruct a regular landscape graph and calculate the accumulated cost of movement from each starting coordinate to points in a cost raster. 
#'@param costsurf Raster cost surface
#'@param scoord matrix 2 column starting coordinate
#'@param snode numeric starting node number
#'@param geocorrect logical use gdistance::geoCorrection() on transition matrix?
#'@importFrom gdistance transition transitionMatrix
#'@importFrom igraph graph.adjacency shortest.paths E
#'@importFrom Matrix cBind rBind
#'@importFrom raster cellFromXY
#'@export
#'@examples
#'dm <- as.matrix(read.delim(system.file("extdata/etherington20120code/cost-surface20x20.txt",
#' package = "irlgraph"), skip = 6, na.strings = "-9999", header = FALSE, sep = " "))
#'
#'costsurf <- raster::raster(nrows=dim(dm)[1],ncols=dim(dm)[2],
#'resolution=1,xmn=0, xmx = dim(dm)[1], ymn = 0, ymx = dim(dm)[2])
#'#neccessary to set resolution
#'
#'costsurf[] <- dm
#'scoord <- matrix(c(10.5, 10.5, 15.5, 15.5), ncol = 2, byrow = TRUE)
#'result <-  gdist_acc(costsurf, scoord)
#'
#'sp::plot(result[[1]]$result)
#'points(scoord[1,1], scoord[1,2])
#'
#'sp::plot(result[[2]]$result)
#'points(scoord[2,1], scoord[2,2])

gdist_acc <- function(costsurf, scoord = NULL, snode = NULL, geocorrect = TRUE){
  
  gdist_acc_surf <- function(costsurf, scoord, snode, geocorrect){
    
    scoord <- matrix(scoord, ncol = 2)
    
    ctrans <- gdistance::transition(costsurf, function(x) 1/mean(x), directions = 16)

    if(geocorrect == TRUE){
      ctrans <- gdistance::geoCorrection(ctrans, "c")
    }
    
    if(is.null(snode) & all(is.null(scoord))){
      stop("Must supply either a starting node or coordinates")
    }
  
    if(is.null(snode)){
      snode <- raster::cellFromXY(costsurf, scoord)
    }
    #print(paste("selected node number:", snode))
  
    mtrans <- gdistance::transitionMatrix(ctrans)
    mtrans <- Matrix::rBind(mtrans,rep(0, nrow(mtrans)))
    mtrans <- Matrix::cBind(mtrans,rep(0, nrow(mtrans)))
    
    startNode <- nrow(mtrans) #extra node to serve as origin
    adjP <- cbind(rep(startNode, times = length(snode)), snode)
    mtrans[adjP] <- Inf
    adjacencyGraph <- igraph::graph.adjacency(mtrans, mode="directed", weighted=TRUE)
    
    igraph::E(adjacencyGraph)$weight <- 1/igraph::E(adjacencyGraph)$weight		
    shortestPaths <- igraph::shortest.paths(adjacencyGraph, v = startNode)[-startNode]
    
    result <- as(costsurf, "RasterLayer")
    result <- raster::setValues(result, shortestPaths)	
    list(result = result, graph = adjacencyGraph)
  }
  
  lapply(1:nrow(scoord), function(x) gdist_acc_surf(costsurf = costsurf, scoord = scoord[x,], snode = snode, geocorrect = geocorrect))
}