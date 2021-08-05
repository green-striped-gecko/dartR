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
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
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

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
# STANDARD ERROR CHECKING
  
# CHECK DATATYPE 
  data.type1 <- utils.check.datatype(x1, verbose=verbose)
  data.type2 <- utils.check.datatype(x2, verbose=verbose)
  
    if((data.type1 != "SNP" && data.type1 != "SilicoDArT")|| (data.type2 != "SNP" && data.type2 != "SilicoDArT")){
      stop (error("Fatal Error: Genlight objects must both be either SNP data or SilicoDArT data (fragment P/A data)"))
    }

    if(data.type1 != data.type2){
      stop (error("Fatal Error: Genlight objects must both be either SNP data or SilicoDArT data (fragment P/A data), not a mixture\n"))
    } 
    if (data.type1 == 1 && data.type2 == 1){
      if(verbose==2){cat(report("  Processing Presence/Absence (SilicoDArT) data\n"))}
    } else if (data.type1 == 2 && data.type2 == 2){
      if(verbose==2){cat(report("  Processing SNP data \n"))}
    } 

  if (is.null(x1) | is.null(x2)){
    stop(error("Fatal Error: Two genlight objects must be provided"))
  }
  
# SCRIPT SPECIFIC ERROR CHECKING
  
  # Check that names and ind.metadata are the same and in the same order
  if (!identical(indNames(x1),indNames(x2))){
    stop(error("Fatal Error: the two genlight objects do not have data for the same individuals in the same order\n"))
  }
  if (!is.null(x1@other$ind.metrics)){
    if (!identical(x1@other$ind.metrics,x1@other$ind.metrics)){
      stop(error("Fatal Error: the two genlight objects do not have identical metadata for the same individuals\n"))
    }
  }
  if (!is.null(x1@other$latlon)){
    if (!identical(x1@other$latlon,x1@other$latlon)){
      stop(error("Fatal Error: the two genlight objects do not have latlon data for the same individuals\n"))
    }
  }
  
# DO THE JOB

  if (verbose >= 2){
    cat(report("  Concatenating two genlight objects,",substitute(x1),"and",substitute(x2),"\n"))
  }
  
  if (verbose >=3) {
    cat("    Number of individuals:",nInd(x1),"\n")
    cat("    First genlight object",substitute(x1),"has",nLoc(x1),"loci\n")
    cat("    Second genlight object",substitute(x2),"has",nLoc(x2),"loci\n")
  }
  
# Join the two genlight objects
  x <- cbind(x1,x2)
  
# Join the locus metrics, if they exist
  if (verbose >= 2){
    cat(report("  Concatenating the locus metrics\n"))
  }
  if (!is.null(x1@other$loc.metrics) & !is.null(x2@other$loc.metrics)){
    x@other$loc.metrics <- rbind(x1@other$loc.metrics,x2@other$loc.metrics)
  } else {
    cat(warn("  Warning: Input genlight objects and/or output genlight object lacks locus metrics\n"))
  } 
  
# Add the ind metrics, assuming they are the same in both genlight objects 
  if (verbose >= 2){
    cat(report("  Adding the individual metrics\n"))
  }
  if (!is.null(x1@other$ind.metrics)){
    x@other$ind.metrics <- x1@other$ind.metrics
  } else if (!is.null(x2@other$ind.metrics)){
      x@other$ind.metrics <- x2@other$ind.metrics
  } else {    
    cat(warn("  Warning: Input genlight objects and/or output genlight object lacks individual metrics\n"))
  }
  
# Add the lat lon metrics, assuming they are the same in both genlight objects  
  if (verbose >= 2){
    cat(report("  Adding the latlons if they exist\n"))
  }
  if (!is.null(x1@other$latlon)){
    x@other$latlon <- x1@other$latlon
  } else if (!is.null(x2@other$latlon)){
    x@other$latlon <- x2@other$latlon
  } else {
    cat(warn("  Warning: Input genlight objects and/or output genlight object lacks latlon data\n"))
  }  
  
# Add the loc metrics flags, set to 1 only if 1 in both genlight objects  
  if (verbose >= 2){
    cat(report("  Setting the locus metrics flags\n"))
  }
  if (!is.null(x1@other$loc.metrics.flags) & !is.null(x2@other$loc.metrics.flags)){
    x@other$loc.metrics.flags$AvgPIC <- x1@other$loc.metrics.flags$AvgPIC*x2@other$loc.metrics.flags$AvgPIC
    x@other$loc.metrics.flags$OneRatioRef <- x1@other$loc.metrics.flags$OneRatioRef*x2@other$loc.metrics.flags$OneRatioRef
    x@other$loc.metrics.flags$OneRatioSnp <- x1@other$loc.metrics.flags$OneRatioSnp*x2@other$loc.metrics.flags$OneRatioSnp
    x@other$loc.metrics.flags$PICRef <- x1@other$loc.metrics.flags$PICRef*x2@other$loc.metrics.flags$PICRef
    x@other$loc.metrics.flags$PICSnp <- x1@other$loc.metrics.flags$PICSnp*x2@other$loc.metrics.flags$PICSnp
    x@other$loc.metrics.flags$CallRate <- x1@other$loc.metrics.flags$CallRate*x2@other$loc.metrics.flags$CallRate
    x@other$loc.metrics.flags$maf <- x1@other$loc.metrics.flags$maf*x2@other$loc.metrics.flags$maf
    x@other$loc.metrics.flags$FreqHets <- x1@other$loc.metrics.flags$FreqHets*x2@other$loc.metrics.flags$FreqHets
    x@other$loc.metrics.flags$FreqHomRef <- x1@other$loc.metrics.flags$FreqHomRef*x2@other$loc.metrics.flags$FreqHomRef
    x@other$loc.metrics.flags$FreqHomSnp <- x1@other$loc.metrics.flags$FreqHomSnp*x2@other$loc.metrics.flags$FreqHomSnp
    x@other$loc.metrics.flags$monomorphs <- x1@other$loc.metrics.flags$monomorphs*x2@other$loc.metrics.flags$monomorphs
    x@other$loc.metrics.flags$OneRatio <- x1@other$loc.metrics.flags$OneRatio*x2@other$loc.metrics.flags$OneRatio
    x@other$loc.metrics.flags$PIC <- x1@other$loc.metrics.flags$PIC*x2@other$loc.metrics.flags$PIC
  } else {
    cat(warn("  Warning: Input genlight objects and/or output genlight object lacks metrics flags. Flags set to zero\n"))
    x@other$loc.metrics.flags$AvgPIC <- 0
    x@other$loc.metrics.flags$OneRatioRef <- 0
    x@other$loc.metrics.flags$OneRatioSnp <- 0
    x@other$loc.metrics.flags$PICRef <- 0
    x@other$loc.metrics.flags$PICSnp <- 0
    x@other$loc.metrics.flags$CallRate <- 0
    x@other$loc.metrics.flags$maf <- 0
    x@other$loc.metrics.flags$FreqHets <- 0
    x@other$loc.metrics.flags$FreqHomRef <- 0
    x@other$loc.metrics.flags$FreqHomSnp <- 0
    x@other$loc.metrics.flags$monomorphs <- 0
    x@other$loc.metrics.flags$OneRatio <- 0
    x@other$loc.metrics.flags$PIC <- 0
  } 
  
# Create the history repository, taking the base from X1 if it exists
  if (verbose >= 2){
    cat(report("  Adding the history\n"))
  }
  if (is.null(x@other$history)) {
    x@other$history <- list(match.call())
  } else {
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()   
  }

  if (verbose >=3) {
    cat("    Number of individuals:",nInd(x1),"\n")
    cat("    Combined genlight object has",nLoc(x),"loci\n")
  }
  
# FLAG SCRIPT END

  if (verbose >= 1) {
    cat(report("Completed:",funname,"\n"))
  }
    
  return(x)
}

