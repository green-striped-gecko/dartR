#' @name gl.drop.ind
#' @title Removes specified individuals from a genlight object
#' @description
#' The script, having deleted individuals, optionally identifies resultant
#' monomorphic loci or loci with all values missing and deletes them
#' (using gl.filter.monomorphs.r). The script also optionally recalculates
#' statistics made redundant by the deletion of individuals from the dataset.
#'
#' The script returns a genlight object with the individuals deleted and,
#' optionally, the recalculated locus metadata.
#' 
#' #' See more about data manipulation in the [tutorial](http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial4dartrdatamanipulation22-dec-21-3.pdf).
#'
#' @param x Name of the genlight object containing SNP genotypes [required].
#' @param ind.list A list of individuals to be removed [required].
#' @param recalc Recalculate the locus metadata statistics [default FALSE].
#' @param mono.rm Remove monomorphic loci [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object with the reduced data
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#'  # SNP data
#'    gl2 <- gl.drop.ind(testset.gl,
#'    ind.list=c('AA019073','AA004859'))
#'    gl2 <- gl.drop.ind(testset.gl,
#'    ind.list=c('AA019073','AA004859'))
#'  # Tag P/A data
#'    gs2 <- gl.drop.ind(testset.gs,
#'    ind.list=c('AA020656','AA19077','AA004859'))
#'    gs2 <- gl.drop.ind(testset.gs, ind.list=c('AA020656'
#'    ,'AA19077','AA004859'),mono.rm=TRUE, recalc=TRUE)
#'
#' @seealso \code{\link{gl.keep.ind}} to keep rather than drop specified
#' individuals
#' @export

gl.drop.ind <- function(x,
                        ind.list,
                        recalc = FALSE,
                        mono.rm = FALSE,
                        verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.1",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
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
        stop(error("Fatal Error: no individuals to drop!\n"))
    }
    
    # DO THE JOB
    
    hold <- x
    
    # Remove individuals
    
    if (verbose >= 2) {
        cat(report("  Deleting specified individuals\n"))
    } 
    if(verbose >= 3){
        cat("  ",paste(ind.list, collapse = ", "),"\n")
    }
    
    # Delete listed individuals, recalculate relevant locus metadata and remove monomorphic loci
    
    # Remove rows flagged for deletion
    inds_to_drop <- which(!x$ind.names %in% ind.list)
    x <- x[inds_to_drop,]

    # Monomorphic loci may have been created
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
    
    # Recalculate statistics
    if (recalc) {
        x <- gl.recalc.metrics(x, verbose = 0)
        if (verbose >= 2) {
            cat(report("  Recalculating locus metrics\n"))
        }
    } else {
        if (verbose >= 2) {
            cat(warn("  Locus metrics not recalculated\n"))
            x <- utils.reset.flags(x, verbose = 0)
        }
    }
    
    # REPORT A SUMMARY
    
    if (verbose >= 3) {
        cat("Summary of recoded dataset\n")
        cat(paste("  No. of loci:", nLoc(x), "\n"))
        cat(paste("  Original No. of individuals:", nInd(hold), "\n"))
        cat(paste("  No. of individuals:", nInd(x), "\n"))
        cat(paste("  Original No. of populations:", nPop(hold), "\n"))
        cat(paste("  No. of populations: ", nPop(x), "\n"))
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
