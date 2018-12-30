#' A utility script to recalculate the the frequency of the homozygous alternate SNP by locus after some populations have been deleted
#'
#' The locus metadata supplied by DArT has FreqHomSnp included,
#' but the frequency of the homozygous alternate will change when some individuals are removed from the dataset. 
#' This script recalculates the FreqHomSnp and places these recalculated values in the appropriate place in the genlight object.
#' Note that the frequency of the homozygote alternate SNPS is calculated from the individuals that could be scored.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The modified genlight object
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' #result <- dartR:::utils.recalc.freqhomsnp(testset.gl)

utils.recalc.freqhomsnp <- function(x, v=2) {

  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  # Work around a bug in adegenet if genlight object is created by subsetting
  x@other$loc.metrics <- x@other$loc.metrics[1:nLoc(x),]
  
  if (v > 0) {
    cat("Starting utils.recalc.freqhomref: Recalculating frequency of homozygotes, alternate allele\n")
  }
  if (is.null(x@other$loc.metrics$FreqHomSnp)) {
    x@other$loc.metrics$FreqHomSnp <- array(NA,nLoc(x))
    if (v >= 3){
      cat("  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n")
    }
  }

  # Do the deed
     t <- as.matrix(x)
     for (i in 1:nLoc(x)) {
       x@other$loc.metrics$FreqHomSnp[i] <- length(which(t[,i] == 2))/(nInd(x)-length(which(is.na(t[,i]))))
     }
     
     if (v > 0) {
       cat("Completed: utils.recalc.freqhomsnp\n\n")
     }
   
   return(x)
}