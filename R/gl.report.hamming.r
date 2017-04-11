#' Calculates the pairwise Hamming distance between DArT trimmed DNA sequences
#'
#' Hamming distance can be computed 
#' by exploiting the fact that the dot product of two binary vectors x and (1 – y) 
#' counts the corresponding elements that are different between x and y.
#' This approach can also be used for vectors that contain more than two possible 
#' values at each position (e.g. A, C, T or G).
#'
#' If a pair of DNA sequences are of differing length, the longer is truncated.
#'
#' The algorithm is that of Johann de Jong 
#' (https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/)
#' as implimented in utils.hamming.r
#'
#' @param gl -- genlight object [required]
#' @return Histogram of Hamming distance for the gl object
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl.report.hamming(testset.gl)

gl.report.hamming <- function(gl) {
  
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
  d <- rep(NA,(nLoc(x)-1))
  pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  getTxtProgressBar(pb)
  for (i in 1:(nLoc(x)-1)){
    for (j in ((i+1):nLoc(x))){
      count <- count + 1
      d[count] <- utils.hamming(s[i],s[j])
    }
    setTxtProgressBar(pb, i/nLoc(x))
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