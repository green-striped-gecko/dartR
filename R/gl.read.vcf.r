#' Converts a vcf file into a genlight object
#'
#' This function needs package vcfR, please install it. The converted genlight
#' object does not have individual metrics. You need to add them 'manually' to
#' the other$ind.metrics slot.
#' @param vcffile A vcf file (works only for diploid data) [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return A genlight object.
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' obj <- gl.read.vcf(system.file('extdata/test.vcf', package='dartR'))
#' }

gl.read.vcf <- function(vcffile,
                        verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jackson",
                     verbosity = verbose)
    
    x = NULL
    
    if (!(requireNamespace("vcfR", quietly = TRUE))) {
        stop(error("To use this function you need to install package: vcfR."))
    } else {
        vcf <- vcfR::read.vcfR(file = vcffile, verbose = verbose)
        x <- vcfR::vcfR2genlight(vcf)
        ploidy(x) <- 2
        x <- gl.compliance.check(x)
        
        # add history
        x@other$history <- list(match.call())
        x <- gl.recalc.metrics(x)
        
        if (verbose > 2) {
            cat(
                important(
                    "Genlight object does not have individual metrics. You need to add them 'manually' to the @other$ind.metrics slot.\n"
                )
            )
        }
        return(x)
    }
    
}
