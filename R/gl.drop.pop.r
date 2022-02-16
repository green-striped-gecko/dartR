#' @name gl.drop.pop
#' @title Removes specified populations from a genlight object
#' @description
#' Individuals are assigned to populations based on the specimen metadata
#' file (csv) used with \code{\link{gl.read.dart}}.
#' The script, having deleted populations, optionally identifies resultant
#' monomorphic loci or loci with all values missing and deletes them
#' (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from
#' the dataset.
#'
#' The script returns a genlight object with the new population assignments and
#' the recalculated locus metadata.
#'
#' @param x Name of the genlight object containing SNP genotypes or Tag P/A data
#' (SilicoDArT) [required].
#' @param pop.list A list of populations to be removed [required].
#' @param as.pop Temporarily assign another metric to represent population for
#' the purposes of deletions [default NULL].
#' @param recalc Recalculate the locus metadata statistics [default FALSE].
#' @param mono.rm Remove monomorphic loci [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object with the reduced data
#'
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#'  # SNP data
#'    gl2 <- gl.drop.pop(testset.gl,
#'    pop.list=c('EmsubRopeMata','EmvicVictJasp'),verbose=3)
#'    gl2 <- gl.drop.pop(testset.gl, pop.list=c('EmsubRopeMata','EmvicVictJasp'),
#'    mono.rm=TRUE,recalc=TRUE)
#'    gl2 <- gl.drop.pop(testset.gl, pop.list=c('Male','Unknown'),as.pop='sex')
#'  # Tag P/A data
#'    gs2 <- gl.drop.pop(testset.gs, pop.list=c('EmsubRopeMata','EmvicVictJasp'))
#'
#' @seealso \code{\link{gl.keep.pop}} to keep rather than drop specified
#'  populations
#' @export

gl.drop.pop <-  function(x,
                         pop.list,
                         as.pop = NULL,
                         recalc = FALSE,
                         mono.rm = FALSE,
                         verbose = NULL) {
    hold <- x
    
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
    
    # Population labels assigned?
    if (is.null(as.pop)) {
        if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
            if (verbose >= 2) {
                cat(
                    warn(
                        "  Warning: Population assignments not detected, running compliance check\n"
                    )
                )
            }
            x <- gl.compliance.check(x, verbose = 0)
        }
    }
    
    # Assign the new population list if as.pop is specified
    pop.hold <- pop(x)
    if (!is.null(as.pop)) {
        if (as.pop %in% names(x@other$ind.metrics)) {
            pop(x) <- unname(unlist(x@other$ind.metrics[as.pop]))
            if (verbose >= 2) {
                cat(report(
                    "  Temporarily assigning",
                    as.pop,
                    "as population\n"
                ))
            }
        } else {
            stop(
                error(
                    "Fatal Error: individual metric assigned to 'pop' does not exist. Check names(gl@other$loc.metrics) and select again\n"
                )
            )
        }
    }
    
    if (verbose >= 2) {
        cat(report("  Checking for presence of nominated populations\n"))
    }
    for (case in pop.list) {
        if (!(case %in% popNames(x))) {
            if (verbose >= 2) {
                cat(
                    warn(
                        "  Warning: Listed population",
                        case,
                        "not present in the dataset -- ignored\n"
                    )
                )
            }
            pop.list <- pop.list[!(pop.list == case)]
        }
    }
    if (length(pop.list) == 0) {
        stop(error("Fatal Error: no populations listed to drop!\n"))
    }
    
    # DO THE JOB
    
    # Remove populations
    
    if (verbose >= 2) {
        cat("  Deleting populations",
            paste(pop.list, collapse = ", "),
            "\n")
    }
    
    # Delete listed populations, recalculate relevant locus metadata and remove monomorphic loci
    # Remove rows flagged for deletion
    pops_to_drop <- which(!x$pop %in% pop.list)
    x <- x[pops_to_drop,]
    pop.hold <- pop.hold[pops_to_drop]

    # Monomorphic loci may have been created
    x@other$loc.metrics.flags$monomorphs == FALSE
    
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
        if (!is.null(as.pop)) {
            cat("  Summary of recoded dataset\n")
            cat(paste("    No. of loci:", nLoc(x), "\n"))
            cat(paste("    No. of individuals:", nInd(x), "\n"))
            cat(paste(
                "    No. of levels of",
                as.pop,
                "remaining: ",
                nPop(x),
                "\n"
            ))
            cat(paste("    Original no. of populations", nPop(hold), "\n"))
            cat(paste(
                "    No. of populations remaining: ",
                length(unique((
                    pop.hold
                ))),
                "\n"
            ))
        } else {
            cat("  Summary of recoded dataset\n")
            cat(paste("    No. of loci:", nLoc(x), "\n"))
            cat(paste("    No. of individuals:", nInd(x), "\n"))
            cat(paste("    Original no. of populations", nPop(hold), "\n"))
            cat(paste("    No. of populations remaining: ", nPop(x), "\n"))
        }
    }
    
    # Reassign the initial population list if as.pop is specified
    
    if (!is.null(as.pop)) {
        pop(x) <- pop.hold
        if (verbose >= 2) {
            cat(report(
                "  Restoring population assignments to initial state\n"
            ))
        }
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
