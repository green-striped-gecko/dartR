#' Filters loci for which the SNP has been trimmed from the sequence tag along with the adaptor 
#'
#' This function checks the position of the SNP within the trimmed sequence tag and identifies those for which the SNP position is outside
#' the trimmed sequence tag. This can happen, rarely, when the sequence containing the SNP resembles the adaptor.
#' 
#' The SNP genotype can still be used in most analyses, but functions like gl2fasta() will present challenges if the SNP has been trimmed from
#' the sequence tag.
#' 
#' Not fatal, but should apply this filter before gl.filter.secondaries, for obvious reasons.
#' 
#' @param x -- name of the genlight object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2, unless specified using gl.set.verbosity]
#' @return A new genlight object with the recalcitrant loci deleted
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' result <- gl.filter.overshoot(testset.gl, verbose=3)

gl.filter.overshoot <- function(x, verbose=NULL) {

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
    stop("Fatal Error: genlight object required!")
  }
  
  if (all(x@ploidy == 1)){
    stop("  Detected Presence/Absence (SilicoDArT) data. Please supply a SNP dataset\n")
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
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
    cat("  Deleting loci with SNPs falling outside the trimmed sequence\n")
  }
  # extracting indexes of loci to keep
  index <- which((snpos <= nchar(trimmed)) == TRUE)
  # loci to keep
  xx <- x[, index]
  # updating loc.metrics
  xx@other$loc.metrics <- x@other$loc.metrics[index, ]
 
# ADD TO HISTORY
  nh <- length(xx@other$history)
  xx@other$history[[nh + 1]] <- match.call()
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
  return(xx)
}

