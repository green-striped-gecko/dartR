#' Checks a gl object to see if it complies with dartR expectations, and amends to comply if necessary
#'
#' A genelight object used by dartR has a number of requirements that allow functions within the package to operate
#' correctly. The genlight object comprises
#' 
#' (a) The SNP genotypes or Tag Presence/Absence data (SilicoDArT);
#' (b) An associated dataframe (gl@other$loc.metrics) containing the locus metrics (e.g. Call Rate, Repeatability, etc);
#' (c) An associated dataframe (gl@other$ind.metrics) containing the individual/sample metrics (e.g. sex, latitude, longitude, etc);
#' (d) A specimen identity field (indNames(gl)) with the unique labels applied to each individual/sample;
#' (e) A population assignment (popNames) for each individual/specimen;
#' (f) Flags that indicate whether or not calculable locus metrics have been updated.
#' 
#' This function will check to see that the genlight object conforms to expectation in regard to the above requirements,
#' and if it does not, will rectify it.
#' 
#' @param x -- name of the input genlight object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A genlight object that conforms to the expectations of dartR
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' x <- gl.compliance.check(testset.gl)

gl.compliance.check <- function (x, verbose=NULL) {

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
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
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }
 
    if (all(x@ploidy == 1)){
      if (verbose >= 2) cat("  Processing Presence/Absence (SilicoDArT) data\n")
      data.type <- "SilicoDArT"
    } else if (all(x@ploidy == 2)){
      if (verbose >= 2) cat("  Processing a SNP dataset\n")
      data.type <- "SNP"
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }

# DO THE JOB
  
  # Check that the data exist, and that they are restricted to the appropriate values
  
  if (data.type == "SNP"){
    mat <- as.matrix(x)
    scores <- c(0,1,2,NA)
    if (verbose >= 2){cat("  Checking coding of SNPs\n")}
    if (max(mat) %in% scores){
      if (verbose >= 1){cat("    SNP data scored NA, 0, 1 or 2 confirmed\n")}
    } else {
      if (verbose >= 1){cat("    Error: SNP data must be scored NA, 0 or 1 or 2, revisit data input\n")}
    }
  } else {
    mat <- as.matrix(x)
    scores <- c(0,1,NA)
    if (verbose >= 2){cat("  Checking coding of Tag P/A data\n")}
    if (max(mat) %in% scores){
      if (verbose >= 1){cat("    Tag P/A data (SilicoDArT) scored 1, 0 (present or absent) confirmed\n")}
    } else {
      if (verbose >= 1){cat("    Error: Tag P/A data (SilicoDArT) must be scored 0 for absent or 1 for present, revisit data input\n")}
    }
  }
  
  # Check for the locus metrics flags, and create if they do not exist.
  
  if (verbose >= 2){cat("  Checking locus metrics flags\n")}
  if(is.null(x@other$loc.metrics.flags)){
    if (verbose >= 1){cat("    Creating locus metrics flags\n")}
    x <- utils.reset.flags(x,set=FALSE,verbose=0)
  } else {
    if (verbose >= 1){cat("    Locus metrics flags confirmed\n")}
  }
  
  # Check that the verbosity flag exists, and if not, create it
  
  if (verbose >= 2){cat("  Checking verbosity flag\n")}
  if(is.null(x@other$verbose)){
    if (verbose >= 1){cat("    Creating verbosity flag\n")}
    x <- utils.reset.verbosity(x,value=2,verbose=0)
  } else {
    if (verbose >= 1){cat("    Verbosity flag confirmed\n")}
  }
  
  # Calculate locus metrics
  
  if (verbose >= 2){cat("  Recalculating locus metrics\n")}
  x <- gl.recalc.metrics(x,verbose=verbose)
  
  # Check that the number of values in the loc.metrics dataframe is the same as the number of loci

  if (nLoc(gl)!=nrow(gl@other$loc.metrics)) { 
    cat("  The number of rows in the loc.metrics table does not match the number of loci! This is potentially a major problem if there is a mismatch of the loci with the metadata. Trace back to identify the cause.\n")  
  }
  
  # Check that the individual metrics exist, and if not, create the df
  
  if (verbose >= 2){cat("  Checking for individual metrics\n")}
  if(is.null(x@other$loc.metrics)){
    if (verbose >= 1){cat("    Creating a slot for individual metrics\n")}
    x@other$ind.metrics$id <- indNames(x)
  } else {
    if (verbose >= 1){cat("    Individual metrics confirmed\n")}
  }
  
  # Check that the population variable exists, and if it does not, create it with a single population 'A'
  
  if (verbose >= 2){cat("  Checking for population assignments\n")}
  if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
    if (verbose >= 1){cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
    pop(x) <- array("pop1",dim = nInd(x))
    pop(x) <- as.factor(pop(x))
  } else {
    if (verbose >= 1){cat("  Population assignments confirmed\n")}
  }
  
  # ADD TO HISTORY
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()  
  
# FLAG SCRIPT END
  
  if (verbose > 0) {
    cat("Completed:", funname, "\n")
  }

   return(x)
}
