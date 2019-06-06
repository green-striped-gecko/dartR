#' Convert a genlight object to format suitable for input to Bayescan
#'
#' The output text file contains the snp data and relevant BAyescan command lines to guide input.
#' 

#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default bayescan.txt]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2bayescan(testset.gl)

#' @references Foll M and OE Gaggiotti (2008) A genome scan method to identify selected loci appropriate for both dominant and codominant markers: A Bayesian perspective. Genetics 180: 977-993.

gl2bayescan <- function(x, outfile="bayescan.txt", outpath=tempdir(), verbose=2) {

# TIDY UP FILE SPECS

  outfilespec <- file.path(outpath, outfile)
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

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# DO THE JOB

  if (verbose >= 2) {cat(paste("Extacting SNP data and creating records for each individual\n"))}
  
# Prepare the data  
  mat <- gl.percent.freq(x, verbose=verbose)
  mat <- mat[order(mat$popn),]

# Create the bayescan input file  
  if (verbose >= 2) {cat(paste("Writing text input file for Bayescan",outfilespec,"\n"))}
  sink(outfilespec)
  
  cat(paste0("[loci]=",nLoc(x)),"\n\n")
  cat(paste0("[populations]=",nPop(x)),"\n\n")

  for (i in 1:nPop(x)){
    cat(paste0("[pop]=",i),"\n")
    popi <- mat[mat$popn==mat$popn[i],]
    for (j in 1:length(popi$popn)){
      cat(j,(2*popi$nobs[j]),2,popi$sum[j],(2*popi$nobs[j]-popi$sum[j]),"\n")
    }
    cat("\n")
  }
  
  sink()

  if (verbose >= 3) {cat(paste("Records written to",outfilespec,"\n"))}

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(NULL)

}


