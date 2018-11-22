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
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The reduced genlight dataset
#' @export
#' @author Arthur Georges (Post to https://groups.google.com/d/forum/dartr)
#' @examples
#' f <- gl.filter.maf(testset.gl, threshold=0.05)


gl.filter.maf <- function(x, threshold=0.01, v=2) {
  
# ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.drop.pop.r!\n"); stop("Execution terminated\n")
  }
  
  if (v < 0 | v > 5){
    cat("Warning: Verbosity must take on an integer value between 0 and 5, set to 2\n")
    v <- 2
  }
  
  if (v > 0) {
    cat("Starting gl.report.maf: Minimum Allele Frequency\n")
  }
  
  if (threshold > 0.5 | threshold <= 0) {
    cat("Warning: threshold must be in the range (0,0.5], but usually small, set to 0.05\n")
    threshold <- 0.05
  }

# Recalculate the relevant loc.metrics
  
  if (v >= 3) {cat("  Removing monomorphic loci and recalculating FreqHoms and FreqHets\n")}

  x <- utils.recalc.maf(x,v=v)
  
  # Remove loci with NA count <= 1-threshold
  index <- x@other$loc.metrics$maf >= threshold
  x2 <- x[ ,index]
  x2@other$loc.metrics <- x@other$loc.metrics[index,]
 
  if (v > 2) {
    cat("  Initial number of loci:", nLoc(x), "\n")
    cat("    Number of loci deleted:", nLoc(x) - nLoc(x2), "\n")
    cat("  Final number of loci:", nLoc(x2), "\n")    
  }
  
  if (v > 0) {
        cat("Completed gl.filter.maf\n\n")
  }
  
  return(x2)
}  
