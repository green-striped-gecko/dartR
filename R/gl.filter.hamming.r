#' Filters loci based on pairwise Hamming distance between sequence tags
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
#' Only one of two loci are retained if their Hamming distance is less that a specified
#' percentage. 5 base differences out of 100 bases is a 20% Hamming distance.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param threshold -- a threshold Hamming distance for filtering loci [default threshold <= 0.2]
#' @param rs -- number of bases in the restriction enzyme recognition sequence [default = 4]
#' @param pb -- switch to output progress bar [default FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2, unless specified using gl.set.verbosity]
#' @return a genlight object filtered on Hamming distance.
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' # SNP data
#'   result <- gl.filter.hamming(testset.gl, threshold=0.25, verbose=3)


gl.filter.hamming <- function(x, threshold=0.2, rs=5, pb=FALSE, verbose=NULL) {
  
  n0 <- nLoc(x)
  
# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
  # STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    cat("  Processing Presence/Absence (SilicoDArT) data\n")
  } else if (all(x@ploidy == 2)){
    cat("  Processing a SNP dataset\n")
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)!")
  }

# FUNCTION SPECIFIC ERROR CHECKING

  if(length(x@other$loc.metrics$TrimmedSequence) == 0) {
    stop("Fatal Error: Data must include Trimmed Sequences\n")
  }
  if(threshold < 0 || threshold > 1){
    cat("  Warning: Parameter 'threshold' must be an integer between 0 and 1, set to 0.2\n")
    threshold = 0.2
  }

# DO THE JOB
  
  if (verbose >= 3) {
    cat("  Note: Hamming distance ranges from zero (sequence identity) to 1 (no bases shared at any position)\n")
    cat("  Note: Calculating pairwise Hamming distances between trimmed reference sequence tags\n")
  }
  
  x@other$loc.metrics$TrimmedSequence <- as.character(x@other$loc.metrics$TrimmedSequence)

  count=0
  nL <- nLoc(x)
  index <- rep(TRUE,(nL-1))
  if (pb) {
    pbar <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(pbar)
  }
  if (verbose >= 2) {
    cat("  Calculating Hamming distances between sequence tags\n")
    cat("  Filtering loci with Hamming Distance is less than",threshold,"\n")
  }
  for (i in 1:(nL-1)){
    s1 <- x@other$loc.metrics$TrimmedSequence[i]
    for (j in ((i+1):nL)){
      count <- count + 1
      s2 <- x@other$loc.metrics$TrimmedSequence[j]
      if(utils.hamming(s1,s2,r=rs) <= threshold) {
        index[i] <- FALSE
        if (verbose >= 3){cat(" Deleting:",locNames(x)[i],locNames(x)[j],"\n")}
        break
      }
    }
    if (pb)  setTxtProgressBar(pbar, i/(nL-1))
  }
  
  x <- x[,(index)]
  # That pesky genlight bug
  x@other$loc.metrics <- x@other$loc.metrics[(index),]
  
  
  # REPORT A SUMMARY
  if (verbose >= 3){
    cat("\n  Summary of filtered dataset\n")
    cat(paste("    Initial No. of loci:",n0,"\n"))
    cat(paste("    Hamming d >",threshold,"\n"))
    cat(paste("    Loci deleted",(n0-nLoc(x)),"\n"))
    cat(paste("    Final No. of loci:",nLoc(x),"\n"))
    cat(paste("    No. of individuals:", nInd(x),"\n"))
    cat(paste("    No. of populations: ", length(levels(factor(pop(x)))),"\n\n"))
  }
  
# ADD TO HISTORY
  
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()      
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(x)
  
}
