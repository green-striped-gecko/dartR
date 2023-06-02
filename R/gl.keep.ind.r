#' @name gl.keep.ind
# Preliminaries -- Parameter specifications -------------- 
#' @title Removes all but the specified individuals from a dartR genlight object
#' @description
#' This script deletes all individuals apart from those listed (ind.list).
#' Monomorphic loci and loci that are scored all NA are optionally deleted (mono.rm=TRUE). 
#' The script also optionally recalculates locus metatdata statistics to accommodate
#' the deletion of individuals from the dataset (recalc=TRUE).
#'
#' The script returns a dartR genlight object with the retained individuals 
#' and the recalculated locus metadata. The script works with both genlight objects
#' containing SNP genotypes and Tag P/A data (SilicoDArT).
#'
#' @param x Name of the genlight object [required].
#' @param ind.list A list of individuals to be retained [required].
#' @param recalc If TRUE, recalculate the locus metadata statistics [default FALSE].
#' @param mono.rm If TRUE, remove monomorphic and all NA loci [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].

#' @export
#' @return A reduced dartR genlight object
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#Examples -------------
#' @examples
#'   # SNP data
#'     gl2 <- gl.keep.ind(testset.gl, ind.list=c('AA019073','AA004859'))
#'   # Tag P/A data
#'    gs2 <- gl.keep.ind(testset.gs, ind.list=c('AA020656','AA19077','AA004859'))
# See also ------------
#' @seealso \code{\link{gl.drop.pop}} to drop rather than keep specified populations
# --------------
# Function 
gl.keep.ind <- function(x,
                        ind.list,
                        recalc = FALSE,
                        mono.rm = FALSE,
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
    for (case in ind.list) {
        if (!(case %in% indNames(x))) {
            cat(
                warn(
                    "  Warning: Listed individual",
                    case,
                    "not present in the dataset -- ignored\n"
                )
            )
            ind.list <- ind.list[!(ind.list == case)]
        }
    }
    if (length(ind.list) == 0) {
        stop(error("Fatal Error: no individuals listed to keep!\n"))
    }
    
# DO THE JOB -------------
    
    hold <- x
    
    if (verbose >= 2) {
        cat(report(
            "  Deleting all but the listed individuals",
            paste(ind.list, collapse = ", "),
            "\n"
        ))
    }
    
    # Delete all but the listed individuals, recalculate relevant locus metadata and remove monomorphic loci
    
    # Remove rows flagged for deletion
    ind.to.keep <- which(x$ind.names %in% ind.list)
    x <- x[ind.to.keep,]

    # Monomorphic loci may have been created -------------------
    x@other$loc.metrics.flags$monomorphs <- FALSE
    
    # Remove monomorphic loci
    if (mono.rm) {
        if (verbose >= 2) {
            cat(report("  Deleting monomorphic loc\n"))
        }
        x <- gl.filter.monomorphs(x, verbose = 0)
    }
    # Check monomorphs have been removed
    if (x@other$loc.metrics.flags$monomorphs == FALSE) {
        if (verbose >= 2) {
            cat(warn(
                "  Warning: Resultant dataset may contain monomorphic loci\n"
            ))
        }
    }
    
    # Recalculate statistics ------------------
    if (recalc) {
        x <- gl.recalc.metrics(x, verbose = 0)
        if (verbose >= 2) {
            cat(report("  Recalculating locus metrics\n"))
        }
    } else {
        if (verbose >= 2) {
            cat(report("  Locus metrics not recalculated\n"))
            x <- utils.reset.flags(x, verbose = 0)
        }
    }
    
# REPORT A SUMMARY -------------
    
    if (verbose >= 3) {
        cat("Summary of recoded dataset\n")
        cat(paste("  No. of loci:", nLoc(x), "\n"))
        cat(paste("  Original No. of individuals:", nInd(hold), "\n"))
        cat(paste("  No. of individuals:", nInd(x), "\n"))
        cat(paste("  Original No. of populations:", nPop(hold), "\n"))
        cat(paste("  No. of populations: ", nPop(x), "\n"))
    }
    
# ADD TO HISTORY -----------------
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
    
# FLAG SCRIPT END ----------------
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    # End Block
    return(x)
}
