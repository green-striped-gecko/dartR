#' Calculates the Hamming distance between two DArT trimmed DNA sequences
#'
#' Hamming distance is calculated as the number of base differences between two 
#' sequences which can be expressed as a count or a proportion. Typically, it is
#' calculated between two sequences of equal length. In the context of DArT
#' trimmed sequences, which differ in length but which are anchored to the left
#' by the restriction enzyme recognition sequence, it is sensible to compare the
#' two trimmed sequences starting from immediately after the common recognition
#' sequence and terminating at the last base of the shorter sequence. 
#' 
#' The Hamming distance between the rows of a matrix can be computed quickly 
#' by exploiting the fact that the dot product of two binary vectors x and (1-y)
#' counts the corresponding elements that are different between x and y.
#' This matrix multiplication can also be used for matrices with more than two possible 
#' values, and different types of elements, such as DNA sequences.
#'
#' The function calculates the Hamming distance between all columns of a 
#' matrix X, or two matrices X and Y. Again matrix multiplication is used, this
#' time for counting, between two columns x and y, the number of cases in which 
#' corresponding elements have the same value (e.g. A, C, G or T). This counting 
#' is done for each of the possible values individually, while iteratively adding 
#' the results. The end result of the iterative adding is the sum of all 
#' corresponding elements that are the same, i.e. the inverse of the Hamming 
#' distance. Therefore, the last step is to subtract this end result H from the 
#' maximum possible distance, which is the number of rows of matrix X.
#' 
#' If the two DNA sequences are of differing length, the longer is truncated. The 
#' initial common restriction enzyme recognition sequence is ignored.
#'
#' The algorithm is that of Johann de Jong 
#' \url{https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/}
#'
#' @param str1 -- string containing the first sequence [required]
#' @param str2 -- string containing the second sequence [required]
#' @param r -- number of bases in the restriction enzyme recognition sequence [default = 4]
#' @return Hamming distance between the two strings
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})

utils.hamming <- function(str1, str2, r=4) {
  
  # Build = Jacob
  
  # Make the strings the same length and remove the recognition sequence
  strmin <- min(nchar(str1),nchar(str2))
  str1 <- substr(str1,r,strmin)
  str2 <- substr(str2,r,strmin)
  # Split the strings into characters and hold each in a vector
  x <- strsplit(str1,"")[[1]]
  y <- strsplit(str2,"")[[1]]
  # Apply the de Jong algorithm
  uniqs <- union(x, y)
  H <- t(x == uniqs[1]) %*% (y == uniqs[1])
  for ( uniq in uniqs[-1] ) {
      H <- H + t(x == uniq) %*% (y == uniq)
  }
  # Calculate the hamming distance as a proportion
  result <- as.numeric(length(x) - H)/strmin
  
  return(result)
}

