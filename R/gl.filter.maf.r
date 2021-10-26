#' @title gl.filter.maf
#' @title Filter loci on the basis of minor allele frequency (MAF) in a genlight
#'  {adegenet} object
#' @description
#' This script calculates the minor allele frequency for each locus and updates
#' the locus metadata for FreqHomRef, FreqHomSnp, FreqHets and MAF (if it
#' exists). It then uses the updated metadata for MAF to filter loci.
#'
#' Note the this filter applies to MAF calculated across all individuals,
#'  without regard to population structure. It is a means of removing overall
#'  rare alleles. To apply this to single populations, use seppop and lapply.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param threshold Threshold MAF -- loci with a MAF less than the threshold
#' will be removed [default 0.01].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#' @return The reduced genlight dataset
#' @export
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' result <- gl.filter.monomorphs(testset.gl)
#' result <- gl.filter.maf(result, threshold=0.05, verbose=3)

gl.filter.maf <- function(x,
                          threshold = 0.01,
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
    
    # Work around a bug in adegenet if genlight object is created by subsetting
    if (nLoc(x) != nrow(x@other$loc.metrics)) {
        stop(
            error(
                "The number of rows in the loc.metrics table does not match the number of loci in your genlight object!"
            )
        )
    }
    
    # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
        if (verbose >= 2) {
            cat(
                report(
                    "  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n"
                )
            )
        }
        pop(x) <- array("pop1", dim = nInd(x))
        pop(x) <- as.factor(pop(x))
    }
    
    # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose = 0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
        cat(warn("  Warning: genlight object contains monomorphic loci\n"))
    }
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (threshold > 0.5 | threshold <= 0) {
        cat(
            warn(
                "  Warning: threshold must be in the range (0,0.5], but usually small, set to 0.05\n"
            )
        )
        threshold <- 0.05
    }
    
    # DO THE JOB
    
    # Recalculate the relevant loc.metrics
    
    if (verbose >= 3) {
        cat(report(
            "  Removing monomorphic loci and recalculating FreqHoms and FreqHets\n"
        ))
    }
    
    x <- utils.recalc.maf(x, verbose = 0)
    
    # Remove loci with NA count <= 1-threshold
    index <- x@other$loc.metrics$maf >= threshold
    x2 <- x[, index]
    x2@other$loc.metrics <- x@other$loc.metrics[index,]
    
    if (verbose > 2) {
        cat("  Initial number of loci:", nLoc(x), "\n")
        cat("  Number of loci deleted:", nLoc(x) - nLoc(x2), "\n")
        cat("  Final number of loci:", nLoc(x2), "\n")
    }
    
    # ADD TO HISTORY
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x2)
}
