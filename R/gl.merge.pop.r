#' @name gl.merge.pop
#' @title Merge two or more populations in a genlight object into one population
#' @description
#' Individuals are assigned to populations based on the specimen metadata data
#' file (csv) used with gl.read.dart().
#'
#' This script assigns individuals from two nominated populations into a new
#' single population. It can also be used to rename populations.
#'
#' The script returns a genlight object with the new population assignments.
#'
#' @param x Name of the genlight object containing SNP genotypes [required].
#' @param old A list of populations to be merged [required].
#' @param new Name of the new population [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object with the new population assignments.
#' @export
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#'    gl <- gl.merge.pop(testset.gl, old=c('EmsubRopeMata','EmvicVictJasp'), new='Outgroup')

gl.merge.pop <- function(x,
                         old = NULL,
                         new = NULL,
                         verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jackson",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    if (verbose >= 1) {
        if (length(old) == 1) {
            cat(report("  Renaming a population\n"))
        } else if (length(old) > 1) {
            cat(report("  Merging a list of populations into one\n"))
        } else {
            stop(error(
                "Fatal Error: At least one old population label must be provided\n"
            ))
        }
    }
    
    # SCRIPT SPECIFIC ERROR TESTING
    
    if (is.null(new)) {
        stop(error("Fatal Error: A new population label must be specified\n"))
    }
    if (class(x) != "genlight") {
        stop(error(
            "Fatal Error: genlight object required for gl.keep.pop.r!\n"
        ))
    }
    if (verbose >= 2) {
        if (length(old) == 1) {
            cat("  Renaming", old, "as", new, "\n")
        } else {
            cat("  Merging", old, "into", new, "\n")
        }
    }
    
    # DO THE JOB
    
    # Merge or rename
    
    for (i in 1:length(old)) {
        levels(pop(x))[levels(pop(x)) == old[i]] <- new
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
