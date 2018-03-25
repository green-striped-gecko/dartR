#' A utility script to recalculate the the frequency of the heterozygous SNPs by locus after some populations have been deleted
#'
#' The locus metadata supplied by DArT has FreqHets included,
#' but the frequency of the heterozygotes will change when some individuals are removed from the dataset. 
#' This script recalculates the FreqHets and places these recalculated values in the appropriate place in the genlight object.
#' Note that the frequency of the homozygote reference SNPS is calculated from the individuals that could be scored.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The modified genlight object
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' result <- utils.recalc.freqhets(testset.gl)

utils.recalc.freqhets <- function(x, v=2) {

  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for utils.recalc.freqhets!\n"); stop("Execution terminated\n")
  }
  if (v > 0) {
    cat("Starting utils.recalc.freqhets: Recalculating frequency of heterozygotes\n")
  }

  # Do the deed
     t <- as.matrix(x)
     for (i in 1:nLoc(x)) {
       x@other$loc.metrics$FreqHets[i] <- length(which(t[,i] == 1))/(nInd(x)-length(which(is.na(t[,i]))))
     }

     if (v > 0) {
       cat("Completed utils.recalc.freqhets\n\n")
     }
     
   return(x)
}