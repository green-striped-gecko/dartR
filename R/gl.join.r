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

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Juliette"

# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x1@other$verbose)){ 
      verbose <- x1@other$verbose
    } else { 
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
  
  if(class(x1)!="genlight" | class(x2)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
    if (all(x1@ploidy == 1 )){
      if(verbose==2){cat("  Processing Presence/Absence (SilicoDArT) data in genlight object 1,",substitute(x1),"\n")}
      data.type1 <- "SilicoDArT"
    } else if (all(x1@ploidy == 2)){
      if(verbose==2){cat("  Processing SNP data in genlight object 1,",substitute(x1),"\n")}
      data.type1 <- "SNP"
    } else {
      stop ("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data) in genlight object 1,",substitute(x1))
    }

    if (all(x2@ploidy == 1)){
      if(verbose==2){cat("  Processing Presence/Absence (SilicoDArT) data in genlight object 2,",substitute(x2),"\n")}
      data.type2 <- "SilicoDArT"
    } else if (all(x2@ploidy == 2)){
      if(verbose==2){cat("  Processing SNP data in genlight object 2,",substitute(x2),"\n")}
      data.type2 <- "SNP"
    } else {
      stop ("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data) in genlight object 2,",substitute(x2))
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
    cat("  Concatenating two genlight objects,",substitute(x1),"and",substitute(x2),"\n")
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
    cat("  Concatenating the locus metrics\n")
  }
  if (!is.null(x1@other$loc.metrics) & !is.null(x2@other$loc.metrics)){
    x@other$loc.metrics <- rbind(x1@other$loc.metrics,x2@other$loc.metrics)
  } else {
    cat("    Warning: Input genlight objects and/or output genlight object lacks locus metrics\n")
  } 
  
# Add the ind metrics, assuming they are the same in both genlight objects 
  if (verbose >= 2){
    cat("  Adding the individual metrics\n")
  }
  if (!is.null(x1@other$ind.metrics)){
    x@other$ind.metrics <- x1@other$ind.metrics
  } else if (!is.null(x2@other$ind.metrics)){
      x@other$ind.metrics <- x2@other$ind.metrics
  } else {    
    cat("    Warning: Input genlight objects and/or output genlight object lacks individual metrics\n")
  }
  
# Add the lat lon metrics, assuming they are the same in both genlight objects  
  if (verbose >= 2){
    cat("  Adding the latlongs if they exist\n")
  }
  if (!is.null(x1@other$latlon)){
    x@other$latlon <- x1@other$latlon
  } else if (!is.null(x2@other$latlon)){
    x@other$latlon <- x2@other$latlon
  } else {
    cat("    Warning: Input genlight objects and/or output genlight object lacks latlong data\n")
  }  
  
# Add the loc metrics flags, set to 1 only if 1 in both genlight objects  
  if (verbose >= 2){
    cat("  Setting the locus metrics flags\n")
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
    cat("    Warning: Input genlight objects and/or output genlight object lacks metrics flags. Flags set to zero\n")
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
    cat("  Adding the history\n")
  }
  if (is.null(x@other$history)) {
    x@other$history <- list(match.call())
  } else {
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()   
  }
  
# Create a verbosity flag, set to the max for X1 and X2 
  if (verbose >= 2){
    cat("  Carrying forward the verbosity\n")
  }
  if (!is.null(x1@other$verbose) & !is.null(x1@other$verbose)){
    x@other$verbose <- max(x1@other$verbose,x2@other$verbose)
  } else if (!is.null(x1@other$verbose)){
    x@other$verbose <- x1@other$verbose
  } else if (!is.null(x2@other$verbose)){
    x@other$verbose <- x2@other$verbose
  } else {  
    cat("    Warning: Input genlight objects and output genlight object lacks verbosity setting, verbosity set to 2\n")
    x@other$verbose <- 2
  }
  
  if (verbose >=3) {
    cat("    Number of individuals:",nInd(x1),"\n")
    cat("    Combined genlight object has",nLoc(x),"loci\n")
  }
  
# FLAG SCRIPT END

  if (verbose >= 1) {
    cat("Completed:",funname,"\n")
  }
    
  return(x)
}

