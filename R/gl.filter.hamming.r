#' Filters loci in a genlight object based on pairwise Hamming distance between sequence tags
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
#' Only one of two loci are retained if their Hamming distance is less that a specified
#' percentage. 5 base differences out of 100 bases is a 20% Hamming distance.
#'
#' @param gl -- genlight object [required]
#' @param t -- a threshold Hamming distance for filtering loci [default 0.2]
#' @return a genlight object filtered on Hamming distance.
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl <- gl.filter.hamming(testset.gl, t=0.35)

gl.filter.hamming <- function(gl=gl, t=0.2) {
  
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
  
  x@other$loc.metrics$TrimmedSequence <- as.character(x@other$loc.metrics$TrimmedSequence)
  
  cat("Calculating pairwise Hamming Distances\n")
  count=0
  flag <- rep(FALSE,(nLoc(x)-1))
  pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  getTxtProgressBar(pb)
  for (i in 1:(nLoc(x)-1)){
    for (j in ((i+1):nLoc(x))){
      count <- count + 1
      s1 <- x@other$loc.metrics$TrimmedSequence[i]
      s2 <- x@other$loc.metrics$TrimmedSequence[j]
      if(utils.hamming(s1,s2) <= (t*100)) {
        flag[i] <- TRUE
        # cat("Deleting locus\n")
      }
    }
    setTxtProgressBar(pb, i/(nLoc(x)-1))
  }
  index <- !flag
  x <- x[,(index)]
  # That pesky genlight bug
  x@other$loc.metrics <- x@other$loc.metrics[(index),]

  return <- x
}