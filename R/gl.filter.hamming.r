#' Filters loci in a genlight object based on pairwise Hamming distance between sequence tags
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
#' by exploiting the fact that the dot product of two binary vectors x and (1 – y) 
#' counts the corresponding elements that are different between x and y.
#' This approach can also be used for vectors that contain more than two possible 
#' values at each position (e.g. A, C, T or G).
#'
#' If a pair of DNA sequences are of differing length, the longer is truncated.
#'
#' The algorithm is that of Johann de Jong 
#' \code{https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/}
#' as implimented in utils.hamming.r
#' 
#' Only one of two loci are retained if their Hamming distance is less that a specified
#' percentage. 5 base differences out of 100 bases is a 20% Hamming distance.
#'
#' @param gl -- genlight object [required]
#' @param t -- a threshold Hamming distance for filtering loci [default 0.2]
#' @param rs -- number of bases in the restriction enzyme recognition sequence [default = 4]
#' @return a genlight object filtered on Hamming distance.
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl <- gl.filter.hamming(testset.gl, t=0.25)

gl.filter.hamming <- function(gl=gl, t=0.2, rs=5) {
  
  x <- gl
  n0 <- nLoc(x)
  
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
  
  cat("Hamming distance ranges from zero (sequence identity) to 1 (no bases shared at any position)\n")
  cat("Calculating pairwise Hamming distances between trimmed reference sequence tags\n")
  count=0
  #flag <- rep(FALSE,(nLoc(x)-1))
  nL <- nLoc(x)
  index <- rep(TRUE,(nL-1))
  pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  getTxtProgressBar(pb)
  for (i in 1:(nL-1)){
    s1 <- x@other$loc.metrics$TrimmedSequence[i]
    for (j in ((i+1):nL)){
      count <- count + 1
      s2 <- x@other$loc.metrics$TrimmedSequence[j]
      if(utils.hamming(s1,s2,r=rs) <= t) {
        index[i] <- FALSE
        break
      }
    }
    setTxtProgressBar(pb, i/(nL-1))
  }
  #index <- flag
  x <- x[,(index)]
  # That pesky genlight bug
  x@other$loc.metrics <- x@other$loc.metrics[(index),]
  
  # REPORT A SUMMARY
  cat("\n\nSummary of filtered dataset\n")
  cat(paste("  Initial No. of loci:",n0,"\n"))
  cat(paste("  Hamming d >",t,"\n"))
  cat(paste("  Loci deleted",(n0-nLoc(x)),"\n"))
  cat(paste("  Final No. of loci:",nLoc(x),"\n"))
  cat(paste("  No. of individuals:", nInd(x),"\n"))
  cat(paste("  No. of populations: ", length(levels(factor(pop(x)))),"\n"))

  return <- x
}