#' A utility script to recalculate the frequency of the homozygous alternate
#' SNP by locus after some populations have been deleted
#'
#' The locus metadata supplied by DArT has FreqHomSnp included, but the
#' frequency of the homozygous alternate will change when some individuals are
#' removed from the dataset.
#'
#' This script recalculates the FreqHomSnp and places these recalculated values
#'  in the appropriate place in the genlight object.
#'
#' Note that the frequency of the homozygote alternate SNPS is calculated from
#' the individuals that could be scored.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report [default 2].
#' @return The modified genlight object.
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @seealso \code{utils.recalc.metrics} for recalculating all metrics,
#' \code{utils.recalc.callrate} for recalculating CallRate,
#' \code{utils.recalc.freqhomref} for recalculating frequency of homozygous
#' reference, \code{utils.recalc.avgpic} for recalculating AvgPIC,
#' \code{utils.recalc.freqhet} for recalculating frequency of heterozygotes,
#' \code{gl.recalc.maf} for recalculating minor allele frequency,
#' \code{gl.recalc.rdepth} for recalculating average read depth
#' @examples
#' #out <- utils.recalc.freqhomsnp(testset.gl)

utils.recalc.freqhomsnp <- function(x,
                                    verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # Check monomorphs have been removed up to date
    if (x@other$loc.metrics.flags$monomorphs == FALSE) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  Warning: Dataset contains monomorphic loci which will be included in the ",
                    funname,
                    " calculations\n"
                )
            )
        }
    }
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (is.null(x@other$loc.metrics$FreqHomSnp)) {
        x@other$loc.metrics$FreqHomSnp <- array(NA, nLoc(x))
        if (verbose >= 2) {
            cat(
                report(
                    "  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n"
                )
            )
        }
    }
    
    # DO THE JOB
    
    t <- as.matrix(x)
    if (verbose >= 2) {
        cat(report("  Recalculating locus metric freqHomSnp\n"))
    }
    x@other$loc.metrics$FreqHomSnp <- colMeans(t == 2, na.rm = T)
    x@other$loc.metrics.flags$FreqHomSnp <- TRUE
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x)
}
