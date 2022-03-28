#' Calculates a similarity (distance) matrix for individuals on the proportion of
#'  shared alleles
#'
#' This script calculates an individual based distance matrix. It uses an C++
#'  implementation, so package Rcpp needs to be installed and it is therefore
#'   really fast (once it has compiled the function after the first run).
#'
#' @param x Name of the genlight containing the SNP genotypes [required].
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' #takes some time at the first run of the function...
#' \dontrun{
#' res <- gl.propShared(bandicoot.gl)
#' res[1:5,1:7] #show only a small part of the matrix
#' }

gl.propShared <- function(x) {
  # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "Rcpp"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it."
    ))
  }
  
  xx <- as.matrix(x)
  glpropSharedC <- function() {
    
  }  #to hack package checking...
  Rcpp::cppFunction(
    "NumericMatrix glpropSharedC(NumericMatrix x) {
  int nrow = x.nrow();
  NumericMatrix out(nrow,nrow);
  for (int i=0; i<(nrow-1); i++) {
     for (int j=(i+1); j<nrow; j++) {
     out(j,i) = 1-mean( abs( na_omit(x(i,_)- x(j,_) ))/2);
     }
  }
  return out;
}"
  )
  res <- glpropSharedC(xx)
  res <- as.matrix(as.dist(res))
  diag(res) <- 1
  colnames(res) <- rownames(res) <- indNames(x)
  res
}
