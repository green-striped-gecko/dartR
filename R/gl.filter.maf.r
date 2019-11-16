#' Filter loci on the basis of minor allele frequency (MAF) in a genlight {adegenet} object
#'
#' This script calculates the minor allele frequency for each locus and updates the locus
#' metadata for FreqHomRef, FreqHomSnp, FreqHets and MAF (if it exists). It then uses 
#' the updated metadata for MAF to filter loci.
#' 
#' Note the this filter applies to MAF calculated across all individuals, without regard
#' to population structure. It is a means of removing overall rare alleles. To apply this to
#' single populations, use sepPop and lapply.
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param threshold -- threshold MAF -- loci with a MAF less than the threshold will be removed [default 0.01]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2, unless specified using gl.set.verbosity]
#' @return The reduced genlight dataset
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' result <- gl.filter.monomorphs(testset.gl)
#' result <- gl.filter.maf(result, threshold=0.05, verbose=3)

gl.filter.maf <- function(x, threshold=0.01, verbose=NULL) {
  
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
      cat("Starting",funname,"[Build =",build,"\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
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

# FUNCTION SPECIFIC ERROR CHECKING

  if (threshold > 0.5 | threshold <= 0) {
    cat("  Warning: threshold must be in the range (0,0.5], but usually small, set to 0.05\n")
    threshold <- 0.05
  }

# DO THE JOB

# Recalculate the relevant loc.metrics
  
  if (verbose >= 3) {cat("  Removing monomorphic loci and recalculating FreqHoms and FreqHets\n")}

  x <- utils.recalc.maf(x,verbose=0)
  
  # Remove loci with NA count <= 1-threshold
  index <- x@other$loc.metrics$maf >= threshold
  x2 <- x[ ,index]
  x2@other$loc.metrics <- x@other$loc.metrics[index,]
 
  if (verbose > 2) {
    cat("    Initial number of loci:", nLoc(x), "\n")
    cat("      Number of loci deleted:", nLoc(x) - nLoc(x2), "\n")
    cat("    Final number of loci:", nLoc(x2), "\n")
  }

# ADD TO HISTORY  
  nh <- length(x2@other$history)
  x2@other$history[[nh + 1]] <- match.call()
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(x2)
}  
