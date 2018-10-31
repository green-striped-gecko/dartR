#' A utility script to recalculate the the minor allele frequency by locus, typically after some populations have been deleted
#'
#' The locus metadata supplied by DArT does not have MAF included, so it is calculated and
#' added to the locus.metadata by this script. The minimum allele frequency will change when
#' some individuals are removed from the dataset. This script recalculates the MAF and 
#' places these recalculated values in the appropriate place in the genlight object.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The modified genlight dataset
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' f <- utils.recalc.maf(testset.gl)


utils.recalc.maf <- function(x, v=2) {
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.drop.pop.r!\n"); stop("Execution terminated\n")
  }
  
  if (v > 0) {
    cat("Starting gl.report.maf: Minimum Allele Frequency\n")
  }
  
  if (v < 0 | v > 5){
    cat("Warning: Verbosity must take on an integer value between 0 and 5, set to 3\n")
    v <- 3
  }
  
  # Recalculate the relevant loc.metrics
  
  if (v >= 3) {cat("  Removing monomorphic loci and recalculating FreqHoms and FreqHets\n")}
  
  x <- gl.filter.monomorphs(x, v = v)
  x <- dartR:::utils.recalc.freqhets(x,v=v)
  x <- dartR:::utils.recalc.freqhomref(x,v=v)
  x <- dartR:::utils.recalc.freqhomsnp(x,v=v)
  
  # Calculate and plot overall MAF
  
  if (v >= 3) {cat("Calculating MAF\n")}
  if (is.null(x@other$loc.metrics$maf)) {
    if (v >= 3){
      cat("  Locus metric maf does not exist, creating slot @other$loc.metrics$maf\n")
    }
    x@other$loc.metrics$maf <- array(NA,nLoc(x))
  } else {
    if (v >= 3){cat("  Recalculating  minor allele frequency\n")}
  }
  
  homref <- x@other$loc.metrics$FreqHomRef
  homalt <- x@other$loc.metrics$FreqHomSnp
  het <- x@other$loc.metrics$FreqHets
  
  for (i in 1:nLoc(x)){
    x@other$loc.metrics$maf[i] <- min((homref[i]*2 + het[i]), (homalt[i]*2 + het[i]))/2
  }
  
  if (v > 0) {
    cat("Completed gl.filter.maf\n\n")
  }
  
  return(x)
}  
