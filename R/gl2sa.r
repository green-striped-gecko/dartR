#' Convert genlight objects to the format used in the SNPassoc package
#'
#' This function exports a genlight object into a SNPassoc object. See \link[SNPassoc]{setupSNP}
#' @param x -- genlight object
#' @export
#' @return Returns an object of class 'snp' to be used with \pkg{SNPassoc}
#' @author Bernd Guber (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl2plink(testset.gl)

gl2sa <- function(x)
{
  pop <- gl2gi(x)
  xxx <- pegas::as.loci(pop)[,-1]
  sa <- SNPassoc::setupSNP(data.frame(xxx), 1:ncol(xxx))
  return(sa)
}
