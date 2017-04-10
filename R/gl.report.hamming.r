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

  if (nLoc(x)==1) {
    cat("Only one locus. No distances calculated.\n")
  } else {
  cat("Calculating pairwise Hamming Distances\n")
  count=0
  d <- rep(NA,(nLoc(x)-1))
  pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  getTxtProgressBar(pb)
  for (i in 1:(nLoc(x)-1)){
    for (j in ((i+1):nLoc(x))){
      count <- count + 1
      s1 <- as.character(x@other$loc.metrics$TrimmedSequence[i])
      s2 <- as.character(x@other$loc.metrics$TrimmedSequence[j])
      d[count] <- utils.hamming(s1,s2)
    }
    setTxtProgressBar(pb, i/nLoc(x))
  }
  }
   hist(d, main="Hamming D", c="red")
   return <- hist(d, main="Hamming D", c="red")
}