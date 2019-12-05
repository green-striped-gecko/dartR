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
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL, plots an interactive 3D plot of the ordination in a separate window
#' @export
#' @importFrom pca3d pca3d
#' @importFrom methods is
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' library(rgl)  #needed for the example
#' pcoa <- gl.pcoa(testset.gl, nfactor=5)
#' gl.pcoa.plot.3d(pcoa, testset.gl, xaxis=1, yaxis=2, zaxis=3)

# Last amended 3-Feb-19

gl.pcoa.plot.3d <- function(x, gl, title= "PCoA", xaxis=1, yaxis=2, zaxis=3,  shape="sphere", radius=2, legend="topright", verbose=2) {

# TIDY UP FILE SPECS

  #outfilespec <- file.path(outpath, outfile)
  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(!is(gl,"genlight")) {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(gl)) | is.na(length(pop(gl))) | length(pop(gl)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(gl) <- array("pop1",dim = nLoc(gl))
      pop(gl) <- as.factor(pop(gl))
    }

# FUNCTION SPECIFIC ERROR CHECKING

# DO THE JOB
                            
  # Extract the coordinates in a form suitable for pca3d
    if (verbose >= 2) {cat("  Extracting coordinates of PCoA solution\n")}
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
    if (verbose >= 2) {cat("  Plotting three specified axes\n")}
   pca3d(coords, shape=shape, radius=radius, group=row.names(x$scores), legend=legend, 
         axe.titles=c(xlab,ylab,zlab))

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(NULL)
}

