#' Filters loci for which the SNP has been trimmed from the sequence tag along with the adaptor 
#'
#' This function checks the position of the SNP within the trimmed sequence tag and identifies those for which the SNP position is outside
#' the trimmed sequence tag. This can happen, rarely, when the sequence containing the SNP resembles the adaptor.
#' 
#' The SNP genotype can still be used in most analyses, but functions like gl2fasta() will present challenges if the SNP has been trimmed from
#' the sequence tag.
#' 
#' @param x -- name of the genlight object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A new genlight object with the recalcitrant loci deleted
#' @importFrom 
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl <- testset.gl
#' nLoc(gl)
#' gl@other$loc.metrics$SnpPosition[10] <- 100
#' gl@other$loc.metrics$SnpPosition[20] <- 100
#' gl@other$loc.metrics$SnpPosition[30] <- 100
#' gl <- gl.filter.overshoot(gl)
#' nLoc(gl)

# Last amended 17-Sep-19

gl.filter.overshoot <- function(x, verbose=2) {

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
  
  if(class(x)=="genlight"){
    if(verbose >= 2){cat("  Genlight object detected\n")}
  } else {
    cat("  Fatal Error: genlight object or distance matrix required!\n"); stop("Execution terminated\n")
  }
  
# SCRIPT SPECIFIC ERROR CHECKING
  
  if(length(x@other$loc.metrics$TrimmedSequence) != nLoc(x)) {
    stop("Fatal Error: Data must include Trimmed Sequences for each loci in a column called 'TrimmedSequence' in the @other$loc.metrics slot.\n")
  }
  if(length(x@other$loc.metrics$SnpPosition) != nLoc(x)) {
    stop("Fatal Error: Data must include position information for each loci.\n")
  }
  
# DO THE JOB

  if (verbose >=2) {cat("  Identifying loci for which the SNP has been trimmed with the adaptor\n")}

  trimmed <- as.character(x@other$loc.metrics$TrimmedSequence)
  snpos <- x@other$loc.metrics$SnpPosition
  # Shift the index for snppos to start from 1 not zero
  snpos <- snpos +1
  # Pull those loci for which the SNP position is greater than the tag length
  xx <- x[,snpos > nchar(trimmed)]
  # Report the number of such loci
  if (verbose >=3){
    cat("  No. of loci with SNP falling outside the trimmed sequence:",nLoc(xx),"\n")
    if(nLoc(xx) > 0){cat("\n",paste(locNames(xx),"\n"))}
  }
  if (verbose >= 2){
    cat("  Deleting loci with SNP falling outside the trimmed sequence\n")
  }
  xx <- x[,snpos <= nchar(trimmed)]
 
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("\nCompleted:",funname,"\n")
  }
    
    return(xx)
}

