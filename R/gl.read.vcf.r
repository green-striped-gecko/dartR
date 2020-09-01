#' Converts a vcf file into a genlight object
#'
#' This function needs package vcfR. At the time of writing (August 2020) the package was no longer available from CRAN. To install the package check their github repository. \url{https://github.com/knausb/vcfR}. 
#' Once installed run with installed=TRUE
#' @param vcffile -- a vcf file (works only for diploid data)
#' @param verbose set to 2
#' @param installed switch to run the function once vcfR package is installed
#' @return A genlight object, with most slots filled.
# @importFrom vcfR read.vcfR vcfR2genlight
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})


gl.read.vcf <- function(vcffile, verbose = 1, installed=FALSE) {

  x=NULL
  if (!installed) cat("This functions requires to have package vcfR installed, which is no longer on CRAN. Check ?gl.read.vcf for help.")
  
  #Change 'if (FALSE) {' below to 'if (TRUE) {' to run the function.
  if (installed) {
  if (!(requireNamespace("vcfR", quietly = TRUE))) {
    stop("To use this function you need to install package: vcfR. Please refer to the help of the function for instructions (?gl.read.vcf).")
  } else {  
  vcf <- vcfR::read.vcfR(file = vcffile, verbose = verbose)
  x <- vcfR::vcfR2genlight(vcf)
  
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
}
