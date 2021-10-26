#' @name gl.reassign.pop
#' @title Assign a individual metric as pop in a genlight \{adegenet\} object
#' @description
#' Individuals are assigned to populations based on the
#' individual/sample/specimen metrics file (csv) used with gl.read.dart().
#'
#' One might want to define the population structure in accordance with another
#' classification, such as using an individual metric (e.g. sex, male or
#' female). This script discards the current population assignments and replaces
#' them with new population assignments defined by a specified individual
#' metric.
#'
#' The script returns a genlight object with the new population assignments.
#' Note that the original population assignments are lost.
#'
#' @param x Name of the genlight object containing SNP genotypes [required].
#' @param as.pop Specify the name of the individual metric to set as the pop
#'  variable [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object with the reassigned populations.
#' @export
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SNP data
#'    popNames(testset.gl)
#'    gl <- gl.reassign.pop(testset.gl, as.pop='sex',verbose=3)
#'    popNames(gl)
#' # Tag P/A data
#'    popNames(testset.gs)
#'    gs <- gl.reassign.pop(testset.gs, as.pop='sex',verbose=3)
#'    popNames(gs)

gl.reassign.pop <- function(x,
                            as.pop,
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
    
    pop(x) <- as.matrix(x@other$ind.metrics[as.pop])
    if (verbose >= 2) {
        cat(report(
            "  Setting population assignments to individual metric",
            as.pop,
            "\n"
        ))
    }
    
    if (verbose >= 3) {
        cat("  Summary of recoded dataset\n")
        cat(paste("    No. of loci:", nLoc(x), "\n"))
        cat(paste("    No. of individuals:", nInd(x), "\n"))
        cat(paste("    No. of populations: ", nPop(x), "\n"))
    }
    
    # ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x)
}
