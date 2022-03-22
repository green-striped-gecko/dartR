#' @name gl.recode.pop
#' @title Recodes population assignments in a genlight object
#' @description
#' This script recodes population assignments and/or deletes populations from a
#' DaRT genlight SNP file based on information provided in a csv population
#' recode file.
#' @details
#' Individuals are assigned to populations based on the specimen metadata data
#' file (csv) used with gl.read.dart(). Recoding can be used to amalgamate
#' populations or to selectively delete or retain populations.
#'
#' The population recode file contains a list of populations in the genlight
#'  object as the first column of the csv file, and the new population
#'  assignments in the second column of the csv file. The keyword Delete used as
#'  a new population assignment will result in the associated specimen being
#'   dropped from the dataset.
#'
#' The script, having deleted populations, optionally identifies resultant
#' monomorphic loci or loci with all values missing and deletes them
#' (using gl.filter.monomorphs.r). The script also optionally recalculates the
#' locus metadata as appropriate.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param pop.recode Name of the csv file containing the population
#' reassignments [required].
#' @param recalc Recalculate the locus metadata statistics if any individuals
#' are deleted in the filtering [default FALSE].
#' @param mono.rm Remove monomorphic loci [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object with the recoded and reduced data.
#' @export
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \donttest{
#'   mfile <- system.file('extdata', 'testset_pop_recode.csv', package='dartR')
#'   nPop(testset.gl)
#'   gl <- gl.recode.pop(testset.gl, pop.recode=mfile, verbose=3)
#'  }
#' @seealso \code{\link{gl.filter.monomorphs}}
#' @seealso \code{\link{gl.recode.pop}}

gl.recode.pop <- function(x,
                          pop.recode,
                          recalc = FALSE,
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
    
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
        stop(error("Fatal Error: Population names not detected\n"))
    }
    
    recode.table <-
        read.csv(pop.recode,
                 header = FALSE,
                 stringsAsFactors = FALSE)
    
    v1 <- unique(pop(x))
    v2 <- unique(recode.table[, 1])
    v1_v2 <- v1[!(v1 %in% v2)]
    # v2_v1 <- v2[! v2 %in% v1]
    l1 <- length(v1)
    l2 <- length(v2)
    
    if (l1 != l2) {
        stop(
            error(
                "Fatal Error: Population names do not agree in number with those in the recode table\n"
            )
        )
    }
    if (!(length(v1_v2) == 0)) {
        stop(
            error(
                "Fatal Error: Some population names have no reassignment specified in the recode table:",
                v1_v2,
                "\n"
            )
        )
    }
    
    # DO THE JOB
    
    if (verbose >= 2) {
        cat(report(
            "  Reassigning entities to populations as per ",
            pop.recode,
            "\n"
        ))
    }
    
    # Store variables
    hold.nLoc <- nLoc(x)
    hold.nInd <- nInd(x)
    hold.nPop <- nPop(x)
    
    # Apply the recode to the populations
    pop.list <- as.character(pop(x))
    ntr <- length(recode.table[, 1])
    for (i in 1:nInd(x)) {
        for (j in 1:ntr) {
            if (pop.list[i] == recode.table[j, 1]) {
                pop.list[i] <- recode.table[j, 2]
            }
        }
    }
    pop(x) <- pop.list
    
    # Remove rows flagged for deletion
    
    if ("delete" %in% popNames(x) | "Delete" %in% popNames(x)) {
        if (verbose >= 2) {
            cat(
                report(
                    "  Deleting individuals in populations flagged for deletion (flagged 'Delete' or 'delete')\n"
                )
            )
        }
        deletions <-
            indNames(x)[tolower(recode.table[, 2]) == "delete"]
        if (verbose == 3) {
            cat("  Dropping\n",
                paste(deletions, collapse = ", "),
                "\n")
            cat("  A total of",
                length(deletions),
                "individuals dropped\n")
        }
        x <-
            gl.drop.pop(x,
                        pop.list = c("Delete", "delete"),
                        verbose = 0)
    }
    
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
        cat("  Summary of recoded dataset\n")
        cat(paste("  Original No. of loci:", hold.nLoc, "\n"))
        cat(paste("    New No. of loci:", nLoc(x), "\n"))
        cat(paste("  Original No. of individuals:", hold.nInd, "\n"))
        cat(paste("    New No. of individuals:", nInd(x), "\n"))
        cat(paste("  Original No. of populations:", hold.nPop, "\n"))
        cat(paste("    New No. of populations:", nPop(x), "\n"))
        if (!recalc) {
            cat(report("  Note: Locus metrics not recalculated\n"))
        }
        if (!mono.rm) {
            cat(report("  Note: Resultant monomorphic loci not deleted\n"))
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
