#' A utility script to recalculate the minor allele frequency by locus,
#' typically after some populations have been deleted
#'
#' The locus metadata supplied by DArT does not have MAF included, so it is
#' calculated and added to the locus.metadata by this script. The minimum allele
#' frequency will change when some individuals are removed from the dataset.
#' This script recalculates the MAF and places these recalculated values in the
#'  appropriate place in the genlight object.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report [default 2].
#' @return The modified genlight dataset.
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @seealso \code{utils.recalc.metrics} for recalculating all metrics,
#' \code{utils.recalc.callrate} for recalculating CallRate,
#' \code{utils.recalc.freqhomref} for recalculating frequency of homozygous
#' reference, \code{utils.recalc.freqhomsnp} for recalculating frequency of
#' homozygous alternate, \code{utils.recalc.freqhet} for recalculating frequency
#' of heterozygotes, \code{gl.recalc.avgpic} for recalculating AvgPIC,
#' \code{gl.recalc.rdepth} for recalculating average read depth
#' @examples
#' #f <- dartR::utils.recalc.maf(testset.gl)

utils.recalc.maf <- function(x,
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
    
    if (is.null(x@other$loc.metrics$maf)) {
        if (verbose >= 3) {
            cat(
                report(
                    "  Locus metric maf does not exist, creating slot @other$loc.metrics$maf\n"
                )
            )
        }
        x@other$loc.metrics$maf <- array(NA, nLoc(x))
    }
    
    # DO THE JOB
    
    if (verbose >= 2) {
        cat(report("  Recalculating FreqHoms and FreqHets\n"))
    }
    
    x <- utils.recalc.freqhets(x, verbose = verbose)
    x <- utils.recalc.freqhomref(x, verbose = verbose)
    x <- utils.recalc.freqhomsnp(x, verbose = verbose)
    
    # Calculate and plot overall MAF
    
    if (verbose >= 2) {
        cat(report("  Recalculating Minor Allele Frequency (MAF)\n"))
    }
    
    alf <- gl.alf(x)[, 2]
    x@other$loc.metrics$maf <- ifelse(alf > 0.5, 1 - alf, alf)
    x@other$loc.metrics.flags$maf <- TRUE
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x)
}
