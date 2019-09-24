#' Calculate a similarity(distance) matrix for individuals on the proportion of shared alleles 
#'
#' This script calculates a individual based distance matrix. It uses an C++ implementation, so package Rcpp needs to be installed.
#' 
#' @param x -- name of the genlight containing the SNP genotypes [required]
#' @importFrom Rcpp cppFunction
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' res <- gl.propShared(testset.gl)
#' res[1:5,1:7] #show only a small part of the matrix


gl.propShared <- function(x) {
  xx <- as.matrix(x)
  cppFunction('NumericMatrix glpropSharedC(NumericMatrix x) {
  int nrow = x.nrow();
  NumericMatrix out(nrow,nrow);
  for (int i=0; i<(nrow-1); i++) {
     for (int j=(i+1); j<nrow; j++) {
     out(j,i) = 1-mean( abs( na_omit(x(i,_)- x(j,_) ))/2);
     }
  }
  return out;
}')
  res <- glpropSharedC(xx)
  res <- as.matrix(as.dist(res))
  diag(res)<- 1
  res
}