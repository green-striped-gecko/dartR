#' Recalculates locus metrics when individuals or populations are deleted from a
#'  genlight \{adegenet\} object
#'
#' When individuals,or populations, are deleted from a genlight object, the
#' locus metrics no longer apply. For example, the Call Rate may be different
#' considering the subset of individuals, compared with the full set. This
#'  script recalculates those affected locus metrics, namely, avgPIC, CallRate,
#' freqHets, freqHomRef, freqHomSnp, OneRatioRef, OneRatioSnp, PICRef and
#'  PICSnp.
#'  Metrics that remain unaltered are RepAvg and TrimmedSeq as they are
#'   unaffected by the removal of individuals.
#'
#' The script optionally removes resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r).
#'
#' The script returns a genlight object with the recalculated locus metadata.
#'
#' @param x Name of the genlight object containing SNP genotypes [required].
#' @param mono.rm If TRUE, removes monomorphic loci [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report
#'   [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object with the recalculated locus metadata.
#' @export
#' @author Custodian: Luis Mijangos (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'   gl <- gl.recalc.metrics(testset.gl, verbose=2)
#' @seealso \code{\link{gl.filter.monomorphs}}

gl.recalc.metrics <- function(x,
                              mono.rm = FALSE,
                              verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # DO THE JOB
    
    # Recalculate statistics
    
    if (datatype == "SNP") {
        x <- utils.recalc.avgpic(x, verbose = verbose)
        x <- utils.recalc.callrate(x, verbose = verbose)
        x <- utils.recalc.maf(x, verbose = verbose)
    }
    if (datatype == "SilicoDArT") {
        x <- utils.recalc.avgpic(x, verbose = verbose)
        x <- utils.recalc.callrate(x, verbose = verbose)
    }
    
    if (verbose >= 2) {
        cat(report("  Locus metrics recalculated\n"))
    }
    
    if (mono.rm) {
        x <- gl.filter.monomorphs(x, verbose = 0)
        if (verbose >= 2) {
            cat(report("  Monomorphic loci deleted\n"))
        }
    }
    
    # ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x)
    
}
