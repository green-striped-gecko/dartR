#' Filter loci that contain private (and fixed alleles) between two populations. 
#'
#' This script is meant to be used prior to \code{gl.nhybrids} to maximise the information content of the snps used to identify hybrids (currently newhybrids does allow only 200 SNPs). The idea is to use first all loci that have fixed alleles between the potential source populations and then "fill up" to 200 loci using loci that have private alleles between those. The functions filters for those loci (if invers is set to TRUE, the opposite is returned (all loci that are not fixed and have no private alleles - not sure why yet, but maybe useful.)
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param pop1 -- name of the first parental population (in quotes)
#' @param pop2 -- name of the second parental population (in quotes)
#' @param invers -- switch to filter for all loci that have no private alleles and are not fixed
#' @return The reduced genlight dataset, containing now only fixed and private alleles
#' @export
#' @author Bernd Gruber & Ella Kelly (University of Melbourne) (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' f <- gl.filter.pa.pop(testset.gl, pop1=)




gl.filter.pa.pop<-function(x, pop1, pop2, invers=FALSE){
  pops <- seppop(x)
  p1 <- as.matrix(pops[[pop1]])
  p2 <- as.matrix(pops[[pop2]])
  p1alf <- colMeans(p1, na.rm = T)/2
  p2alf <- colMeans(p2, na.rm = T)/2
  priv1 <- c(names(p1alf)[p2alf == 0 & p1alf != 0], names(p1alf)[p2alf == 1 & p1alf != 1]) # private alleles for pop 1
  priv2 <-  c(names(p2alf)[p1alf == 0 & p2alf != 0], names(p2alf)[p1alf == 1 & p2alf != 1]) # private alleles for pop 2
  pfLoci<-unique(c(priv1, priv2)) # put all together
  index <- locNames(x) %in% pfLoci
  if (invers) index <- !index
  x <- x[, index]
  x@other$loc.metrics <- x@other$loc.metrics[index, ]
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()
  return(x)
}
