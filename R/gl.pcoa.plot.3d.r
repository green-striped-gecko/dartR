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
#' @param glPca -- name of the glPca object containing the factor scores and eigenvalues [required]
#' @param x -- name of the genlight object or fd object from which the PCoA was generated
#' @param title -- a title for the plot [default "PCoA"]
#' @param xaxis -- identify the x axis from those available in the ordination (xaxis <= nfactors) [default 1]
#' @param yaxis -- identify the y axis from those available in the ordination (yaxis <= nfactors) [default 2]
#' @param zaxis -- identify the z axis from those available in the ordination (zaxis <= nfactors) [default 3]
#' @param shape -- shape of the points, one of sphere, tetrahaedron or cube [default "sphere"]
#' @param radius -- size of the points [default 2]
#' @param legend -- one of bottomright, bottom, bottomleft, left, topleft, top, topright, right, center [default "bottom"]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return NULL, plots an interactive 3D plot of the ordination in a separate window
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' library(rgl)  #needed for the example
#' pcoa <- gl.pcoa(testset.gl, nfactor=5)
#' gl.pcoa.plot.3d(pcoa, testset.gl, xaxis=1, yaxis=2, zaxis=3)


gl.pcoa.plot.3d <- function(glPca, 
                            x, 
                            title= "PCA", 
                            xaxis=1, 
                            yaxis=2, 
                            zaxis=3,  
                            shape="sphere", 
                            radius=2, 
                            legend="topright", 
                            verbose=NULL) {
# CHECK IF PACKAGES ARE INSTALLED
  pkg <- "pca3d"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package",pkg," needed for this function to work. Please   install it.")
  } else {

  
  # TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if(class(x)=="genlight"){
    if (is.null(verbose)){ 
      if(!is.null(x@other$verbose)){ 
        verbose <- x@other$verbose
      } else { 
        verbose <- 2
      }
    }
  }
  if(class(x)=="fd"){
    x <- x$gl
    if (is.null(verbose)){
      verbose <- 2
    }  
  }
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nLoc(x))
      pop(x) <- as.factor(pop(x))
    }

# FUNCTION SPECIFIC ERROR CHECKING

# DO THE JOB
                            
  # Extract the coordinates in a form suitable for pca3d
    if (verbose >= 2) {cat("  Extracting coordinates of PCoA solution\n")}
    coords <- cbind(glPca$scores[,xaxis],glPca$scores[,yaxis],glPca$scores[,zaxis])
  
  # Convert the eigenvalues to percentages
   s <- sum(glPca$eig)
   e <- round(glPca$eig*100/s,1)
  # Create labels for the axes
   xlab <- paste0("P", xaxis, " (",e[xaxis],"%)")
   ylab <- paste0("P", yaxis, " (",e[yaxis],"%)")
   zlab <- paste0("P", zaxis, " (",e[zaxis],"%)")
  # Create a title
   #t <- paste("PCoA plot of Axes", xaxis, yaxis,"and",zaxis)
  # Set the row labels to the population names
   row.names(glPca$scores) <- as.character(pop(x))
  # Plot
    if (verbose >= 2) {cat("  Plotting three specified axes\n")}
   pca3d::pca3d(coords, shape=shape, radius=radius, group=row.names(glPca$scores), legend=legend, 
         axe.titles=c(xlab,ylab,zlab))

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(NULL)
  }
}
