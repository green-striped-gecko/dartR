#' Subsamples n loci from a genlight object and return it as a genlight object
#'
#' This is a support script, to subsample a genlight \{adegenet\} object based
#'  on loci. Two methods are used to subsample, random and based on information
#'  content.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param n Number of loci to include in the subsample [required].
#' @param method Method: 'random', in which case the loci are sampled at random;
#' or 'pic', in which case the top n loci ranked on information content are
#' chosen. Information content is stored in AvgPIC in the case of SNP data and in
#'  PIC in the the case of presence/absence (SilicoDArT) data [default 'random'].
#' @param mono.rm Delete monomorphic loci before sampling [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object with n loci
#' @export
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SNP data
#'   gl2 <- gl.subsample.loci(testset.gl, n=200, method='pic')
#' # Tag P/A data
#'   gl2 <- gl.subsample.loci(testset.gl, n=100, method='random')

gl.subsample.loci <- function(x,
                              n,
                              method = "random",
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
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (mono.rm) {
        if (x@other$loc.metrics.flags$monomorphs == FALSE) {
            if (verbose >= 2) {
                cat(report("  Deleting monomorphic loci\n"))
            }
            x <- gl.filter.monomorphs(x, verbose = 0)
        } else {
            if (verbose >= 2) {
                cat(report("  Zero monomorphic loci, none deleted\n"))
            }
        }
    }
    
    # Check monomorphs have been removed
    if (x@other$loc.metrics.flags$monomorphs == FALSE) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  Warning: Dataset contains monomorphic loci which will be included in the",
                    funname,
                    "selections\n"
                )
            )
        }
    }
    
    if (!method == "pic" & !method == "random") {
        if (verbose >= 1) {
            cat(
                warn(
                    "  Warning: parameter method must be set to 'pic' or 'random', set to random\n"
                )
            )
        }
        method <- "random"
    }
    
    if (n <= 0 | n > nLoc(x)) {
        stop(
            error(
                "Fatal Error: subsample size must be a postive integer >= 1 or <=",
                nLoc(x),
                "\n"
            )
        )
    }
    
    # DO THE JOB
    
    if (datatype == "SilicoDArT") {
        pic <- x@other$loc.metrics$PIC
    }
    if (datatype == "SNP") {
        pic <- x@other$loc.metrics$AvgPIC
    }
    
    if (method == "random") {
        if (verbose >= 2) {
            cat(report(
                "  Subsampling at random",
                n,
                "loci from",
                class(x),
                "object",
                "\n"
            ))
        }
        randsel <- sample(1:nLoc(x), n, replace = FALSE)
        
          x.new <- x[, randsel]
          x.new@other$loc.metrics <- x@other$loc.metrics[randsel, ]
        
        if (verbose >= 3) {
            cat(report("  No. of loci retained =", ncol(x.new), "\n"))
        }
        
    } else if (method == "PIC" | method == "pic") {
        x <- x[, order(-pic)]
        
        x.new <- x[, 1:n]
        x.new@other$loc.metrics <- x@other$loc.metrics[1:n, ]
    
        if (verbose >= 3) {
            cat(report("  No. of loci retained =", ncol(x.new), "\n"))
        }
        
    } else {
        stop(error("Fatal Error: method must be 'random' or 'pic'\n"))
    }
    
    # ADD TO HISTORY
    nh <- length(x.new@other$history)
    x.new@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x.new)
    
}
