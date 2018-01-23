#' A utility script to recalculate the the frequency of the homozygous reference SNP by locus after some populations have been deleted
#'
#' The locus metadata supplied by DArT has FreqHomRef included,
#' but the frequency of the homozygous reference will change when some individuals are removed from the dataset. 
#' This script recalculates the FreqHomRef and places these recalculated values in the appropriate place in the genlight object.
#' Note that the frequency of the homozygote reference SNPS is calculated from the individuals that could be scored.
#'
#' @param gl -- name of the genlight object containing the SNP data [required]
#' @param v -- v=0, silent; v=1, low verbosity; v=2, high verbosity [default 1]
#' @return The modified genlight object
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' result <- dartR:::utils.recalc.freqhomref(testset.gl)

utils.recalc.freqhomref <- function(gl, v=1) {
 x <- gl
   
  if(class(x) == "genlight") {
     #cat("Reporting for a genlight object\n")
   } else {
     cat("Fatal Error: Specify a genlight object\n")
     stop()
  }

  # Do the deed
     t <- as.matrix(x)
     for (i in 1:nLoc(x)) {
       x@other$loc.metrics$FreqHomRef[i] <- length(which(t[,i] == 0))/(nInd(x)-length(which(is.na(t[,i]))))
     }

   if (v>0) {cat("FreqHomRef recalculated\n")}
   
   return(x)
}