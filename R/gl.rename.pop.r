#' @name gl.rename.pop
#' @title Renames a population in a genlight object
#' @description
#' Individuals are assigned to populations based on the specimen metadata data
#' file (csv) used with gl.read.dart().
#'
#' This script renames a nominated population.
#'
#' The script returns a genlight object with the new population name.
#'
#' @param x Name of the genlight object containing SNP genotypes [required].
#' @param old Name of population to be changed [required].
#' @param new New name for the population [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object with the new population name.
#' @export
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#'    gl <- gl.rename.pop(testset.gl, old='EmsubRopeMata', new='Outgroup')

gl.rename.pop <- function(x,
                         old = NULL,
                         new = NULL,
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
    
    # SCRIPT SPECIFIC ERROR TESTING
    
    if (is.null(old)) {
        stop(error("Fatal Error: A population to be renamed must be specified\n"))
    }
    if (is.null(new)) {
        stop(error("Fatal Error: A new population label must be specified\n"))
    }
    if (!is(x, "genlight")) {
        stop(error(
            "Fatal Error: genlight object required!\n"
        ))
    }
    if (verbose >= 2) {
            cat("  Renaming", old, "as", new, "\n")
    }
    
    # DO THE JOB

    levels(pop(x))[levels(pop(x)) == old] <- new
    
    # ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x)
}
