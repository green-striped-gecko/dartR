#' Converts a vcf file into a genlight object
#'
#' @param vcffile -- a vcf file (works only for diploid data)
#' @param verbose set to 2
#' @return A genlight object, with most slots filled.
#' @importFrom vcfR read.vcfR vcfR2genlight
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})


gl.read.vcf <- function(vcffile, verbose = 1) {
  vcf <- read.vcfR(file = vcffile, verbose = verbose)
  x <- vcfR2genlight(vcf)
  
  #add history
  x@other$history <- list(match.call())
  x <- gl.recalc.metrics(x)
  
  if (verbose) {
    cat(
      "genlight object does not have individual metrics. You need to add them 'manually' to the @other$ind.metrics slot.\n"
    )
  }
  return(x)
  
}
