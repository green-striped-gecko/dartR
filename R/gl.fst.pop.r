#' Calculate a pairwise fst values for populations in a genlight object
#'
#' This script calculates pairwise fst values based on the implementation in the StAMPP package (?stamppFst). It allows to run bootstrap to estimate probability of fst values to be different from zero. For detailed information please check the help pages (?stamppFst).
#' 
#' @param x -- name of the genlight containing the SNP genotypes [required]
#' @param nboots -- number of bootstraps to perform across loci to generate confidence intervals and p-values
#' @param percent -- the percentile to calculate the confidence interval around [defalut = 95]
#' @param nclusters -- the number of proccesor threads or cores to use during calculations.
#' @return A matrix of distances between populations (class dist), if nboots =1, otherwise a list with Fsts (in a matrix), Pvalues (a matrix of pvalues), Bootstraps results (data frame of all runs). Hint: Use \code{as.matrix(as.dist(fsts))} if you want to have a squared matrix with symmetric entries returned, instead of a dist object.
#' @importFrom StAMPP stamppFst
#' @export
#' @author Bernd Gruber (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl.fst.pop(possums.gl, nboots=1)

gl.fst.pop <- function(x, nboots=100, percent=95, nclusters =1) {

  fsts <- stamppFst(x, nboots=nboots, percent=percent, nclusters = nclusters)
  return (fsts)
}
