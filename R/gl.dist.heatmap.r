#' Represent a distance matrix as a heatmap 
#'
#' The script plots a heat map to represent the distances in the distance or dissimilarity matrix
#'
#' @param d -- name of the distance matrix [required]
#' @param ncolors -- number of colors to display [default 5]
#' @param rank -- if TRUE, then the distance matrix will be ordered, otherwise order will be displayed as given [default TRUE]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @import gclus
#' @return NULL
#' @export
#' @author Francois Gillet, August 2009, modified by Arthur Georges (Post to https://groups.google.com/d/forum/dartr)
#' @examples
#'    dist <- dist(testset.gl)
#'    gl <- gl.dist.heatmap(dist)

gl.dist.heatmap <- function(dist, ncolors=5, rank=TRUE, v=2){
  
  x <- dist
  
  # ERROR CHECKING
  
  if(class(x)!="dist") {
    cat("Fatal Error: distance matrix required!\n"); stop("Execution terminated\n")
  }

  if (v < 0 | v > 5){
    cat("    Warning: verbosity must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    v <- 2
  }
  
  if (ncolors < 0){
    cat("    Warning: ncolors must be a positive integer, set to 5\n")
    ncolors <- 5
  }
  
  if (max(dist) > 1) {
    cat("    Warning: matrix contains distances greater than 1, rescaling\n")
    dist <- dist/max(dist)
  } 
  
  # FLAG SCRIPT START
  
  if (v >= 1) {
    cat("Starting gl.plot.heatmap: Displaying distance matrix\n")
  }

  if (rank) {
    spe.color = dmat.color(1-D, rainbow(nc))
    spe.o = order.single(1-D)
    speo.color = spe.color[spe.o,spe.o]
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
               main="Ordered Distance Matrix")
  }
  else {
    spe.color = dmat.color(1-D, byrank=FALSE, rainbow(nc))
    spe.o = order.single(1-D)
    speo.color = spe.color[spe.o,spe.o]
    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
               main="Unordered Distance Matrix")
  }
  
  return()

}