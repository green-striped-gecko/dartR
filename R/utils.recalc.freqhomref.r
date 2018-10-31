#' A utility script to recalculate the the frequency of the homozygous reference SNP by locus after some populations have been deleted
#'
#' The locus metadata supplied by DArT has FreqHomRef included,
#' but the frequency of the homozygous reference will change when some individuals are removed from the dataset. 
#' This script recalculates the FreqHomRef and places these recalculated values in the appropriate place in the genlight object.
#' Note that the frequency of the homozygote reference SNPS is calculated from the individuals that could be scored.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The modified genlight object
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' #result <- dartR:::utils.recalc.freqhomref(testset.gl)

utils.recalc.freqhomref <- function(x, v=2) {

  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for utils.recalc.freqhomref!\n"); stop("Execution terminated\n")
  }
  if (v > 0) {
    cat("Starting utils.recalc.freqhomref: Recalculating frequency of homozygotes, reference allele\n")
  }
  if (is.null(x@other$loc.metrics$FreqHomRef)) {
    x@other$loc.metrics$FreqHomRef <- array(NA,nLoc(x))
    if (v >= 3){
      cat("  Locus metric FreqHomRef does not exist, creating slot @other$loc.metrics$FreqHomRef\n")
    }
  }  

  # Do the deed
     t <- as.matrix(x)
     for (i in 1:nLoc(x)) {
       x@other$loc.metrics$FreqHomRef[i] <- length(which(t[,i] == 0))/(nInd(x)-length(which(is.na(t[,i]))))
     }

     if (v > 0) {
       cat("Completed utils.recalc.freqhomref\n\n")
     }
   
   return(x)
}