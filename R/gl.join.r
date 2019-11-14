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
#' gl <- gl.join(x1, x2, verbose=2)
#' nLoc(gl)

gl.join <- function(x1, x2, verbose=NULL) {

# TIDY UP FILE SPECS
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# FLAG SCRIPT START
  # set verbosity
  if (is.null(verbose)) verbose=2
 
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose >= 1){
    cat("Starting",funname,"\n")
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x1)!="genlight" | class(x2)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
    if (all(x1@ploidy == 1 )){
      if(verbose==2){cat("  Processing Presence/Absence (SilicoDArT) data in genlight object 1\n")}
      data.type1 <- "SilicoDArT"
    } else if (all(x1@ploidy == 2)){
      if(verbose==2){cat("  Processing SNP data in genlight object 1\n")}
      data.type1 <- "SNP"
    } else {
      stop ("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data) in genlight object 1")
    }

    if (all(x2@ploidy == 1)){
      if(verbose==2){cat("  Processing Presence/Absence (SilicoDArT) data in genlight object 2\n")}
      data.type2 <- "SilicoDArT"
    } else if (all(x2@ploidy == 2)){
      if(verbose==2){cat("  Processing SNP data in genlight object 2\n")}
      data.type2 <- "SNP"
    } else {
      stop ("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data) in genlight object 2")
    }

  if (data.type1 != data.type2){
    stop("Fatal Error: Cannot join two genlight objects, one with SNP data, the other with Tag Presence/Absence data (SilicoDArT)")
  }
  
  if (is.null(x1) | is.null(x2)){
    stop("Fatal Error: Two genlight objects of the same type must be provided")
  }
  
# SCRIPT SPECIFIC ERROR CHECKING
  
  # Check that names and ind.metadata are the same and in the same order
  if (!identical(indNames(x1),indNames(x2))){
    stop("Fatal Error: the two genlight objects do not have data for the same individuals in the same order\n")
  }
  if (!is.null(x1@other$ind.metrics)){
  if (!identical(x1@other$ind.metrics,x1@other$ind.metrics)){
    stop("  Fatal Error: the two genlight objects do not have identical metadata for the same individuals\n")
  }
  }
  if (!is.null(x1@other$latlon)){
  if (!identical(x1@other$latlon,x1@other$latlon)){
    stop("  Fatal Error: the two genlight objects do not have latlong data for the same individuals\n")
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
  
  # Recalculate metrics
  
  x@other$loc.metrics.flags$monomorphs <- FALSE
  x <- gl.recalc.metrics(x,verbose=min(c(verbose,1)))
  if (verbose>=2) cat("  Locus metrics recalculated\n")
  
  # Check monomorphs
  
  tmp <- gl.filter.monomorphs(x, verbose = 0)
  if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
    cat("  Warning: new genlight object contains monomorphic loci\n")
  }

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("\nCompleted:",funname,"\n")
  }
    
    return(x)
}

