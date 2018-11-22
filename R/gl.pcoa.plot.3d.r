#' 3D interactive plot of the results of a PCoA ordination
#'
#' This script takes output from the ordination undertaken using gl.pcoa() and plots the individuals in 3D space. The visualisation
#' can be rotated with the mouse to examine the structure.
#'
#' The factor scores are taken from the output of gl.pcoa(), an object of class glPca, and the population assignments from the original data file
#' and plots the specimens in a 3D plot.
#'
#' Axes can be specified from the ordination, provided they are within the range of the nfactors value provided to gl.pcoa().
#'
#' This script is essentially a wrapper for function pca3d \{pca3d\} maintained by January Weiner.
#' 
#' @param x -- name of the glPca object containing the factor scores and eigenvalues [required]
#' @param gl -- name of the genlight object from which the PCoA was generated
#' @param title -- a title for the plot [default "PCoA"]
#' @param xaxis -- identify the x axis from those available in the ordination (xaxis <= nfactors) [default 1]
#' @param yaxis -- identify the y axis from those available in the ordination (yaxis <= nfactors) [default 2]
#' @param zaxis -- identify the z axis from those available in the ordination (zaxis <= nfactors) [default 3]
#' @param shape -- shape of the points, one of sphere, tetrahaedron or cube [default "sphere"]
#' @param radius -- size of the points [default 2]
#' @param legend -- one of bottomright, bottom, bottomleft, left, topleft, top, topright, right, center [default "bottom"]
#' @return An interactive 3D plot of the ordination in a separate window
#' @export
#' @importFrom pca3d pca3d
#' @import adegenet
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' pcoa <- gl.pcoa(testset.gl, nfactor=5)
#' gl.pcoa.plot.3d(pcoa, testset.gl, xaxis=1, yaxis=2, zaxis=3)

gl.pcoa.plot.3d <- function(x, gl, title= "PCoA", xaxis=1, yaxis=2, zaxis=3,
                            shape="sphere", radius=2, legend="topright") {

  # Extract the coordinates in a form suitable for pca3d
    coords <- cbind(x$scores[,xaxis],x$scores[,yaxis],x$scores[,zaxis])
  
  # Convert the eigenvalues to percentages
   s <- sum(x$eig)
   e <- round(x$eig*100/s,1)
  # Create labels for the axes
   xlab <- paste0("P", xaxis, " (",e[xaxis],"%)")
   ylab <- paste0("P", yaxis, " (",e[yaxis],"%)")
   zlab <- paste0("P", zaxis, " (",e[zaxis],"%)")
  # Create a title
   #t <- paste("PCoA plot of Axes", xaxis, yaxis,"and",zaxis)
  # Set the row labels to the population names
   row.names(x$scores) <- as.character(pop(gl))
  # Plot 
   pca3d(coords, shape=shape, radius=radius, group=row.names(x$scores), legend=legend, 
         axe.titles=c(xlab,ylab,zlab))
}

