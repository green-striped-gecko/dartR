#' Calculates the Hamming distance between two DArT trimmed DNA sequences
#'
#' The Hamming distance between the rows of a matrix can be computed fast, 
#' by exploiting the fact that the dot product of two binary vectors x and (1 – y) 
#' counts the corresponding elements that are different between x and y.
#' Matrix multiplication can also be used for matrices with more than two possible 
#' values, and different types of elements.
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
#' If the two DNA sequences are of differing length, the longer is truncated.
#'
#' The algorithm is that of Johann de Jong 
#' (https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/)
#'
#' @param str1 -- string containing the first sequence [required]
#' @param str2 -- string containing the second sequence [required]
#' @return Hamming distance between the two strings
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' str1 <- "AcTGGCTAGC"
#' str2 <- "AaTGGCTAG"
#' utils.hamming(str1,str2)

utils.hamming <- function(str1, str2) {
  # Make the strings the same length
  strmin <- min(nchar(str1),nchar(str2))
  str1 <- substr(str1,1,strmin)
  str2 <- substr(str2,1,strmin)
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


