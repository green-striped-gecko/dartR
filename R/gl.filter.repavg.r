#' Filter loci in a genlight \{adegenet\} object based on average repeatability of alleles at a locus
#'
#' SNP datasets generated by DArT have in index, RepAvg, generated by reproducing the data independently for 30\% of loci.
#' RepAvg is the proportion of alleles that give a repeatable result, averaged over both alleles for each locus.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param threshold -- threshold value below which loci will be removed [default 0.99]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return     Returns a genlight object retaining loci with a Repavg greater than the specified threshold deleted.
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.report.repavg(testset.gl)
#' result <- gl.filter.repavg(testset.gl, threshold=0.95, verbose=3)

# Last amended 3-Feb-19

gl.filter.repavg <- function(x, threshold=0.99, verbose=2) {

# TIDY UP FILE SPECS

  funname <- match.call()[[1]]
  cat("Warning: Depreciated. Use gl.filter.repeatability\n")

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Work around a bug in adegenet if genlight object is created by subsetting
      if (nLoc(x)!=nrow(x@other$loc.metrics)) { stop("The number of rows in the loc.metrics table does not match the number of loci in your genlight object!")  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

# FUNCTION SPECIFIC ERROR CHECKING

# DO THE JOB

  if (verbose > 2) {cat("  Note: RepAvg is a DArT statistic reporting repeatability averaged across alleles for each locus. \n\n")}
  
  n0 <- nLoc(x)
  if (verbose > 2) {cat("  Initial no. of loci =", n0, "\n")}

  if(class(x)=="genlight") {
    # Remove SNP loci with RepAvg < threshold
    if (verbose > 1){cat("    Removing loci with RepAvg <",threshold,"\n")}
    x2 <- x[, x@other$loc.metrics["RepAvg"]>=threshold]
    # Remove the corresponding records from the loci metadata
    x2@other$loc.metrics <- x@other$loc.metrics[x@other$loc.metrics["RepAvg"]>=threshold,]
    if (verbose > 2) {cat ("  No. of loci deleted =", (n0-nLoc(x2)),"\n")}
    
  } else if (class(x)=="genind") {
    x2 <- x[,(colSums(is.na(tab((x))))/nInd(x))<(1-threshold)]
    idx <- which((colSums(is.na(tab((x))))/nInd(x))<(1-threshold))
    x2@other$loc.metrics <- x@other$loc.metrics[c(idx),]
    
    # Remove SNP loci with RepAvg < threshold
    if (verbose > 1){cat("    Removing loci with Reproducibility <",threshold,"\n")}
    x2 <- x[, x@other$loc.metrics["Reproducibility"]>=threshold]
    # Remove the corresponding records from the loci metadata
    x2@other$loc.metrics <- x@other$loc.metrics[x@other$loc.metrics["Reproducibility"]>=threshold,]
    if (verbose > 2) {cat ("No. of loci deleted =", (n0-nLoc(x2)),"\n")}
  } else {
    cat("  Fatal Error: genlight or genind objects required for call rate filtering!\n"); stop()
  }
  
  # REPORT A SUMMARY
  if (verbose > 2) {
    cat("  Summary of filtered dataset\n")
    cat(paste("    Repeatability >=",threshold,"\n"))
    cat(paste("    No. of loci:",nLoc(x2),"\n"))
    cat(paste("    No. of individuals:", nInd(x2),"\n"))
    cat(paste("    No. of populations: ", length(levels(factor(pop(x2)))),"\n"))
  }  
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  #add to history
  nh <- length(x2@other$history)
  x2@other$history[[nh + 1]] <- match.call()  
  return(x2)
  
}
