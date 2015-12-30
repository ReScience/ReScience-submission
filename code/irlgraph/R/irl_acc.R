#'@name irl_acc
#'@title Generate an accumulated cost surface from an irregular landscape graph
#'@description Wraps the irl_graph, acc_path, and irl_acc functions
#'@export
#'@param dm matrix 2d cost values
#'@param poicoords matrix 2 column point-of-interest coordinates
#'@param grainprop numeric proportion of cells to be included as "grain cells"
#'@param cutoff numeric between 0 and 1. Determines the threshold at which a cell is determined "very important"
#'@param costsurf raster cost surface
#'@param scoord matrix matrix 2 column starting coordinates
#'@param snode numeric starting node
#'@param irregular logical skip irregular node subsetting routine?
#'@param warn logical print warnings?
#'@description This high level function allows for generation of an accumulated cost surface from a single function call
#'@examples \dontrun{
#'set.seed(123) #make reproducible
#'
#'dm <- as.matrix(read.delim(system.file(
#'"extdata/etherington20120code/cost-surface20x20.txt",
#' package = "irlgraph"), skip = 6, na.strings = "-9999", header = FALSE, sep = " "))
#'
#'costsurf <- raster::raster(nrows=dim(dm)[1],ncols=dim(dm)[2],resolution=1,xmn=0,
#' xmx = dim(dm)[1], ymn = 0, ymx = dim(dm)[2]) #neccessary to set resolution
#'costsurf[] <- dm
#'
#'poicoords <- matrix(c(10.5, 10.5), ncol = 2)
#'
#'result <- irl_acc(dm, poicoords = poicoords, grainprop = 0.25,
#' costsurf = costsurf, scoord = matrix(c(19.5, 0.5, 11.5, 11.5),
#' ncol = 2, byrow = TRUE))
#'
#'result <- irl_acc(dm, poicoords = poicoords, grainprop = 0.25,
#' costsurf = costsurf, scoord = matrix(c(11.5, 11.5),
#' ncol = 2, byrow = TRUE))
#'}

irl_acc <- function(dm, poicoords = NA, cutoff = 0, grainprop = 0.25, costsurf, scoord, snode = NULL, irregular = TRUE, warn = TRUE){
  
  graph <- irl_graph(dm, poicoords = poicoords, grainprop = grainprop, irregular = irregular, cutoff = cutoff)
  
  acc_surf <- function(graph, scoord, snode, costsurf, warn){
    path <- acc_path(graph = graph, scoord = scoord, snode = snode, costsurf = costsurf)
    impute_na(path, costsurf, graph, warn = warn)
  } 
  
  lapply(1:nrow(scoord), function(x) acc_surf(graph, scoord =  scoord[x,], snode = NULL, costsurf = costsurf, warn = warn))
  
}