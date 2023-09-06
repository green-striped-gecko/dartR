#' Converts a genind object into a genlight object
#'
#' @param gi A genind object [required].
#' @param parallel Switch to deactivate parallel version. It might not be worth
#' to run it parallel most of the times [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report [default 2].
#' @return A genlight object, with all slots filled.
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @details
#' Be aware due to ambiguity which one is the reference allele a combination of
#'  gi2gl(gl2gi(gl)) does not return an identical object (but in terms of
#'  analysis this conversions are equivalent)

gi2gl <- function(gi,
                  parallel = FALSE,
                  verbose = NULL) {
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)
    
    # STANDARD ERROR CHECKING
    
    if (!is(gi, "genind")) {
        stop(error("  Fatal Error: genind object required!\n"))
    }
    
    # DO THE JOB
    
    locna <- gi@loc.n.all
    ccc <- 1
    for (i in 2:length(locna)) {
        if (locna[i - 1] == 1) {
            ccc[i] <- ccc[i - 1] + 1
        } else {
            ccc[i] <- ccc[i - 1] + 2
        }
    }
    
    gl <-
        new(
            "genlight",
            gi@tab[, ccc],
            pop = pop(gi),
            other = gi@other,
            ploidy = 2,
            loc.names = locNames(gi),
            ind.names = indNames(gi),
            parallel = parallel
        )
    
    gl <- gl.compliance.check(gl)
    
    if (is.null(gl@other$loc.metrics.flags$monomorphs)) {
        gl@other$loc.metrics.flags$monomorphs <- FALSE
    }
    # FLAG SCRIPT END
    
    # ADD TO HISTORY
    
    nh <- length(gl@other$history)
    gl@other$history[[nh + 1]] <- match.call()
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(gl)
}
