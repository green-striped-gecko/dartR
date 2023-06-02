#' @name gl.keep.loc
# Preliminaries -- Parameter specifications -------------- 
#' @title Removes all but the specified loci from a genlight object
#' @description
#' This function deletes loci that are not specified to keep, and their associated metadata. 
#'
#' The script returns a dartR genlight object with the retained loci. 
#' The script works with both genlight objects
#' containing SNP genotypes and Tag P/A data (SilicoDArT).

#' @param x Name of the genlight object [required].
#' @param loc.list A list of loci to be kept
#' [required, if loc.range not specified].
#' @param first First of a range of loci to be kept
#' [required, if loc.list not specified].
#' @param last Last of a range of loci to be kept
#' [if not specified, last locus in the dataset].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress but not results; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#'
#' @export
#' @return A genlight object with the reduced data
#'
#' @family dartR-base
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
# Examples -------------
#' @examples
#' # SNP data
#'   gl2 <- gl.keep.loc(testset.gl, loc.list=c('100051468|42-A/T', '100049816-51-A/G'))
#' # Tag P/A data
#'   gs2 <- gl.keep.loc(testset.gs, loc.list=c('20134188','19249144'))
# See also ------------
#' @seealso \code{\link{gl.drop.loc}} to drop rather than keep specified loci
#'
# End Block --------------
# Function 
gl.keep.loc <- function(x,
                        loc.list = NULL,
                        first = NULL,
                        last = NULL,
                        verbose = NULL) {
  # Preliminaries -------------
  # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
   # Function-specific error checking -----------        
    if (!is.null(loc.list) && !is.null(first)) {
        flag <- "both"
        if (verbose >= 2) {
            cat(report(
                "  Both a range of loci and a list of loci to keep has been specified\n"
            ))
        }
    } else if (!is.null(loc.list)) {
        flag <- "list"
        if (verbose >= 2) {
            cat(report("  List of loci to keep has been specified\n"))
        }
    } else if (!is.null(first)) {
        flag <- "range"
        if (verbose >= 2) {
            cat(report("  Range of loci to keep has been specified\n"))
        }
    } else {
        cat(
            warn(
                "  Warning: Need to specify either a range of loci to keep, or specific loci to keep\n"
            )
        )
    }
    
    if (flag == "both" || flag == "list") {
        for (case in loc.list) {
            if (!(case %in% locNames(x))) {
                cat(
                    warn(
                        "  Warning: Listed loci",
                        case,
                        "not present in the dataset -- ignored\n"
                    )
                )
                loc.list <- loc.list[!(loc.list == case)]
            }
        }
    }
    
    if (flag == "range") {
        if (first <= 0) {
            cat(warn(
                "  Warning: Lower limit to range of loci cannot be less than 1, set to 1\n)"
            ))
            first <- 1
        }
        if (first > nLoc(x)) {
            cat(
                warn(
                    "  Warning: Upper limit to range of loci cannot be greater than the number of loci, set to",
                    nLoc(x),
                    "\n)"
                )
            )
            last <- nLoc(x)
        }
        if (first > last) {
            cat(warn(
                "  Warning: Upper limit is smaller than lower limit, reversed\n"
            ))
            tmp <- first
            first <- last
            last <- tmp
        }
    }
    
	# DO THE JOB --------------
    # Remove loci ------
    hold <- x
    
    if (verbose >= 2) {
        cat(report("  Deleting all but the specified loci\n"))
    }
    
    # Remove loci if specified
    
    if (!is.null(first) && !is.null(loc.list)) {
        list.from.range <- locNames(x)[first:last]
        loc.list <- unique(c(loc.list, list.from.range))
    } else if (!is.null(first)) {
        loc.list <- locNames(x)[first:last]
    }
    if (length(loc.list) == 0) {
        cat(warn(
            "  Warning: no loci listed to keep! Genlight object returned unchanged\n"
        ))
        x2 <- x
    } else {
      
        # Remove loci flagged for deletion
        x2 <- x[, x$loc.names %in% loc.list]
        x2@other$loc.metrics <- x@other$loc.metrics[x$loc.names %in% loc.list, ]

    }
    
# REPORT A SUMMARY -------------
    # Summary of outcomes --------------       
    if (verbose >= 3) {
        cat("  Summary of recoded dataset\n")
        cat(paste("    Original No. of loci:", nLoc(hold), "\n"))
        cat(paste("    No. of loci deleted:", nLoc(hold) - nLoc(x2), "\n"))
        cat(paste("    No. of loci retained:", nLoc(x2), "\n"))
        # cat(paste(' No. of individuals:', nInd(x2),'\n')) cat(paste(' No. of populations: ', nPop(x2),'\n'))
    }
    
    # ADD TO HISTORY ---------------
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END ----------------
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    # End block --------------
    return(x2)
}
