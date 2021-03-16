#' Converts a vcf file into a genlight object
#'
#' This function needs package vcfR, please install it. The converted genlight object does not have individual metrics. You need to add them 'manually' to the other$ind.metrics slot.
#' @param vcffile -- a vcf file (works only for diploid data)
#' @param verbose set to 2
#' @return A genlight object.
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' test <- gl.read.vcf("your_data.vcf")
#' }

gl.read.vcf <- function(vcffile, verbose = 2) {

  x=NULL

  if (!(requireNamespace("vcfR", quietly = TRUE))) {
    stop("To use this function you need to install package: vcfR.")
  } else {
  vcf <- vcfR::read.vcfR(file = vcffile, verbose = verbose)
  x <- vcfR::vcfR2genlight(vcf)
  ploidy(x) <- 2
  x <- gl.compliance.check(x)

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
  
}
