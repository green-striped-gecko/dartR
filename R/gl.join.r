#' Combines two genlight objects 
#'
#' This function combines two genlight objects and their associated metadata. The history associated with the two genlight objects is cleared
#' from the new genlight object. The individuals/samples must be the same in each genlight object.
#' 
#' The function is typically used to combine datasets from the same service where the files have been split because of size limitations. The
#' data is read in from multiple csv files, then the resultant genlight objects are combined.
#' 
#' @param x1 -- name of the first genlight object [required]
#' @param x2 -- name of the first genlight object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A new genlight object
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' x1 <- testset.gl[,1:100]
#' x1@other$loc.metrics <-  testset.gl@other$loc.metrics[1:100,]
#' nLoc(x1)
#' x2 <- testset.gl[,101:150]
#' x2@other$loc.metrics <-  testset.gl@other$loc.metrics[101:150,]
#' nLoc(x2)
#' gl <- gl.join(gl1, gl2, verbose=2)
#' nLoc(gl)

# Last amended 17-Sep-19

gl.join <- function(x1, x2, verbose=2) {

# TIDY UP FILE SPECS

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
  
  if (class(x1)=="genlight"){
    if(verbose >= 2){cat("  First genlight object detected\n")}
    if (class(x2)=="genlight"){
      if(verbose >= 2){cat("  Second genlight object detected\n")}
    } else{
      cat("  Fatal Error: two genlight objects required! Only one is supplied\n"); stop("Execution terminated\n")
    }
  } else {
    cat("  Fatal Error: two genlight objects required!\n"); stop("Execution terminated\n")
  }

# SCRIPT SPECIFIC ERROR CHECKING
  
  # Check that names and ind.metadata are the same and in the same order
  if (!identical(indNames(x1),indNames(x2))){
    cat("  Fatal Error: the two genlight objects do not have data for the same individuals in the same order\n")
    stop("Execution terminated\n")
  }
  if (!is.null(x1@other$ind.metrics)){
  if (!identical(x1@other$ind.metrics,x1@other$ind.metrics)){
    cat("  Fatal Error: the two genlight objects do not have identical metadata for the same individuals\n")
    stop("Execution terminated\n")
  }
  }
  if (!is.null(x1@other$latlon)){
  if (!identical(x1@other$latlon,x1@other$latlon)){
    cat("  Fatal Error: the two genlight objects do not have latlong data for the same individuals\n")
    stop("Execution terminated\n")
  }
  }
  
# DO THE JOB

  if (verbose >= 2){
    cat("  Concatenating two genlight objects\n")
  }
  
  if (verbose >=3) {
    cat("    Number of individuals:",nInd(x1),"\n")
    cat("    First genlight object has",nLoc(x1),"loci\n")
    cat("    Second genlight object has",nLoc(x2),"loci\n")
  }
  
  x <- cbind(x1,x2)
  if (!is.null(x1@other$loc.metrics)){
    x@other$loc.metrics <- rbind(x1@other$loc.metrics,x2@other$loc.metrics)
  } else {
    cat("  Warning: Input genlight objects and output genlight object lacks locus metrics\n")
  } 
  if (!is.null(x1@other$ind.metrics)){
    x@other$ind.metrics <- x1@other$ind.metrics
  } else {
    cat("  Warning: Input genlight objects and output genlight object lacks individual metrics\n")
  }
  if (!is.null(x1@other$latlon)){
    x@other$latlon <- x1@other$latlon
  } else {
    cat("  Warning: Input genlight objects and output genlight object lacks latlong data\n")
  }  
  
  if (verbose >=3) {
    cat("    Number of individuals:",nInd(x1),"\n")
    cat("    Combined genlight object has",nLoc(x),"loci\n")
  }

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("\nCompleted:",funname,"\n")
  }
    
    return(x)
}

