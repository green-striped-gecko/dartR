#' A utility script to recalculate the OneRatioRef, OneRatioSnp, PICRef,	PICSnp, and	AvgPIC by locus after some 
#' populations have been deleted.
#'
#' The locus metadata supplied by DArT has OneRatioRef, OneRatioSnp, PICRef,	PICSnp, and	AvgPIC included,
#' but the allelec composition will change when some individuals are removed from the dataset and so the initial statistics will
#' no longer apply. This script recalculates these statistics and places the recalculated values in the appropriate place 
#' in the genlight object.
#'
#' @param gl -- name of the genlight object containing the SNP data [required]
#' @param v -- v=0, silent; v=1, low verbosity; v=2, high verbosity [default 1]
#' @return The modified genlight object
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' result <- utils.recalc.avgpic(testset.gl)

utils.recalc.avgpic <- function(gl, v=1) {

  if(class(gl) == "genlight") {
     #cat("Reporting for a genlight object\n")
   } else {
     cat("Fatal Error: Specify a genlight object\n")
     stop()
  }

  # Do the deed
     t <- as.matrix(gl)
     for (i in 1:nLoc(gl)) {
       c0 <- length(t[t[,i]==0 & !is.na(t[,i]),i])
       c1 <- length(t[t[,i]==1 & !is.na(t[,i]),i])
       c2 <- length(t[t[,i]==2 & !is.na(t[,i]),i])
       c <- (c0+c1+c2)
       gl@other$loc.metrics$OneRatioRef[i] <- (c0+c1)/c
       gl@other$loc.metrics$OneRatioSnp[i] <- (c1+c2)/c
       OneRatioRef <- gl@other$loc.metrics$OneRatioRef[i]
       OneRatioSnp <- gl@other$loc.metrics$OneRatioSnp[i]
       ZeroRatioRef <- 1 - OneRatioRef
       ZeroRatioSnp <- 1 - OneRatioSnp
       gl@other$loc.metrics$PICRef[i] <- 1 - ((OneRatioRef*OneRatioRef) + (ZeroRatioRef*ZeroRatioRef))
       gl@other$loc.metrics$PICSnp[i] <- 1 - ((OneRatioSnp*OneRatioSnp) + (ZeroRatioSnp*ZeroRatioSnp))
       gl@other$loc.metrics$avgPIC[i] <- (gl@other$loc.metrics$PICRef[i] + gl@other$loc.metrics$PICSnp[i])/2
     }

   if (v>0) {cat("OneRatioRef, OneRatioSnp, PICRef, PICSnp, and AvgPIC recalculated\n")}
   
   return(gl)
}