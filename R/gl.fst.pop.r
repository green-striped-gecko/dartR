#' Calculates a pairwise Fst values for populations in a genlight object
#'
#' This script calculates pairwise Fst values based on the implementation in the
#' StAMPP package (?stamppFst). It allows to run bootstrap to estimate
#' probability of Fst values to be different from zero. For detailed information
#' please check the help pages (?stamppFst).
#'
#' @param x Name of the genlight containing the SNP genotypes [required].
#' @param nboots Number of bootstraps to perform across loci to generate
#' confidence intervals and p-values [default 100].
#' @param percent Percentile to calculate the confidence interval around
#'  [default 95].
#' @param nclusters Number of processor threads or cores to use during
#' calculations [default 1].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return A matrix of distances between populations (class dist), if nboots =1,
#' otherwise a list with Fsts (in a matrix), Pvalues (a matrix of pvalues),
#' Bootstraps results (data frame of all runs). Hint: Use
#' \code{as.matrix(as.dist(fsts))} if you want to have a squared matrix with
#' symmetric entries returned, instead of a dist object.
#' @importFrom StAMPP stamppFst
#' @export
#' @author Bernd Gruber (bugs? Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' out <- gl.fst.pop(possums.gl, nboots=1)

gl.fst.pop <- function(x,
                       nboots = 100,
                       percent = 95,
                       nclusters = 1,
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
    
    fsts <-
        stamppFst(x,
                  nboots = nboots,
                  percent = percent,
                  nclusters = nclusters)
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    return(fsts)
}
