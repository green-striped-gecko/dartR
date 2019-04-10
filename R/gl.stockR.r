#' Generate a SNP matrix from a genlight \{adegenet\} object for subsequent use in R package
#' stockR.
#'
#' The script extracts the SNP data from the gl object as a matrix and transposes it
#' to comply with the format expected by stockR.
#' 
#' The script returns the transposed SNP matrix
#' 
#' @param x -- name of the genlight object containing SNP genotypes [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A matrix with the SNP scores in the form expected by stockR
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- gl.stockR(testset.gl)

# Last amended 4-Apr-19

gl.stockR <- function (x,verbose = 2) {
  
  # TIDY UP FILE SPECS
  
  funname <- match.call()[[1]]
  
  # FLAG SCRIPT START
  
  if (verbose < 0 | verbose > 5) {
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose > 0) {
    cat("Starting", funname, "\n")
  }
  
  # STANDARD ERROR CHECKING
  
  if (class(x) != "genlight") {
    cat("  Fatal Error: genlight object required!\n")
    stop("Execution terminated\n")
  }
  
  tmp <- gl.filter.monomorphs(x, verbose = 0)
  if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
    cat("  Warning: genlight object contains monomorphic loci\n")
  }
  
  # DO THE JOB
  
  if (verbose >=2) {
    cat("Converting SNP from genlight to stockR format\n")
  }

  snp.matrix <- t(as.matrix(x))
  
# FINISH UP  
  
  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return <- snp.matrix
  
}