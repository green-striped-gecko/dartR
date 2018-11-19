#' Calculates the pairwise Hamming distance between DArT trimmed DNA sequences
#'
#' Hamming distance is calculated as the number of base differences between two 
#' sequences which can be expressed as a count or a proportion. Typically, it is
#' calculated between two sequences of equal length. In the context of DArT
#' trimmed sequences, which differ in length but which are anchored to the left
#' by the restriction enzyme recognition sequence, it is sensible to compare the
#' two trimmed sequences starting from immediately after the common recognition
#' sequence and terminating at the last base of the shorter sequence. 
#' 
#' Hamming distance can be computed 
#' by exploiting the fact that the dot product of two binary vectors x and (1-y)
#' counts the corresponding elements that are different between x and y.
#' This approach can also be used for vectors that contain more than two possible 
#' values at each position (e.g. A, C, T or G).
#'
#' If a pair of DNA sequences are of differing length, the longer is truncated.
#'
#' The algorithm is that of Johann de Jong 
#' \url{https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/}
#' as implimented in utils.hamming.r
#'
#' @param gl -- genlight object [required]
#' @param rs -- number of bases in the restriction enzyme recognition sequence [default = 4]
#' @return Histogram of Hamming distance for the gl object
#' @importFrom graphics hist
#' @importFrom stats sd
#' @export
#' @author Arthur Georges (bugs? Post to https://groups.google.com/d/forum/dartr)
#' @examples
#' gl.report.hamming(testset.gl)

gl.report.hamming <- function(gl, rs=5) {
  
  x <- gl
  
  if(class(x) == "genlight") {
    cat("Analysing a genlight object\n")
  } else {
    cat("Fatal Error: Specify a genlight object\n")
    stop()
  }
  if(length(x@other$loc.metrics$TrimmedSequence) == 0) {
    cat("Fatal Error: Data must include Trimmed Sequences\n"); stop()
  }
  
  s <- as.character(x@other$loc.metrics$TrimmedSequence)

  if (nLoc(x)==1) {
    cat("Only one locus. No distances calculated.\n")
  } else {
    cat("Hamming distance ranges from zero (sequence identity) to 1 (no bases shared at any position)\n")
    cat("Calculating pairwise Hamming distances between trimmed reference sequence tags\n")
  count=0
  nL <- nLoc(x)
  
  # Calculate the number of iterations in loops below to set dimensions of d
  # niter = sum i=1 to nL-1 of (nL - i)
  # niter = sum i=1 to nL-1 of (nL) - sum i=1 to nL-1 of (i)
  # niter = nL(nL-1) - (nL-1)nL/2 [triangle number]
  # niter = nL(nL-1)/2 which seem intuitive
  d <- rep(NA,(((nL-1)*nL)/2))
  
  pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  getTxtProgressBar(pb)
  for (i in 1:(nL-1)){
    for (j in ((i+1):nL)){
      count <- count + 1
      d[count] <- utils.hamming(s[i],s[j],r=rs)
    }
    setTxtProgressBar(pb, i/(nL-1))
  }
  }
   mn <- round(mean(d),2)
   sdev <- round(sd(d),2)
   mi <- round(min(d),2)
   mx <- round(max(d),2)
   cat(paste0("\n\n  Mean Hamming Distance ",mn,"+/-",sdev,"SD (",mi," - ",mx,")\n"))
   hist(d, main="Hamming D", c="red")
   result <- hist(d, plot=FALSE)
   return <- result
}