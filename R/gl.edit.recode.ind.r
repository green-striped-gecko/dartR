#' @name gl.edit.recode.ind
#' @title Create or edit a individual (=specimen) names, create an recode_ind
#'  file amd apply the changes to a genlight object.
#' @description
#' A script to edit individual names in a genlight object, or to create a
#' reassignment table taking the individual labels from a genlight object, or to
#' edit existing individual labels in an existing recode_ind file.
#' @details
#' Renaming individuals may be required when there have been errors in labelling
#'  arising in the process from sample to DArT files. There may be occasions
#'  where renaming individuals is required for preparation of figures. Caution
#'  needs to be exercised because of the potential for breaking the 'chain of
#'  evidence' between the samples themselves and the analyses. Recoding
#'  individuals can also be done with a recode table (csv).
#'
#' This script will input an existing recode table for editing and optionally
#' save it as a new table, or if the name of an input table is not supplied,
#' will generate a table using the individual labels in the parent genlight
#' object.
#'
#' The script, having deleted individuals, optionally identifies resultant
#' monomorphic loci or loci with all values missing and deletes them (using
#' gl.filter.monomorphs.r). The script also optionally recalculates statistics
#' made redundant by the deletion of individuals from the dataset.
#'
#' Use outpath=getwd() or outpath='.' when calling this function to direct
#' output files to your working directory.
#'
#' The script returns a genlight object with the new individual labels and the
#' recalculated locus metadata.
#'
#' @param x Name of the genlight object for which individuals are to be
#'  relabelled [required].
#' @param out.recode.file Name of the file to output the new individual labels
#'  [optional].
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN].
#' @param recalc Recalculate the locus metadata statistics [default TRUE].
#' @param mono.rm Remove monomorphic loci [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return An object of class ('genlight') with the revised individual labels.
#' @import utils
#' @export
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @seealso \code{\link{gl.recode.ind}}, \code{\link{gl.drop.ind}},
#' \code{\link{gl.keep.ind}}
#' @examples
#' \dontrun{
#' gl <- gl.edit.recode.ind(testset.gl)
#' gl <- gl.edit.recode.ind(testset.gl, out.recode.file='ind.recode.table.csv')
#' }

gl.edit.recode.ind <- function(x,
                               out.recode.file = NULL,
                               outpath = tempdir(),
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
    
    outfilespec <- file.path(outpath, out.recode.file)
    
    # DO THE JOB
    
    # Store variables
    hold <- x
    
    # Take assignments from x
    
    if (verbose >= 2) {
        cat(report(
            "  Extracting current individual labels from the x object\n"
        ))
    }
    recode.table <- cbind(indNames(x), indNames(x))
    
    # Create recode table for editting, and bring up the editor
    new <- as.matrix(edit(recode.table))
    # new <- new[,1:2]
    
    # Write out the recode table, if requested
    if (is.null(out.recode.file)) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  Warning: No output table specified, recode table not written to disk\n"
                )
            )
        }
    } else {
        if (verbose >= 2) {
            cat(report(
                paste(
                    "  Writing individual recode table to: ",
                    out.recode.file,
                    "\n"
                )
            ))
        }
        write.table(
            new,
            file = out.recode.file,
            sep = ",",
            row.names = FALSE,
            col.names = FALSE
        )
    }
    
    # Apply the new assignments
    ind.list <- as.character(indNames(x))
    
    ntr <- length(new[, 1])
    for (i in 1:nInd(x)) {
        for (j in 1:ntr) {
            if (ind.list[i] == new[j, 1]) {
                ind.list[i] <- new[j, 2]
            }
        }
    }
    # Assigning new populations to x
    if (verbose >= 2) {
        cat(report("  Assigning new individual (=specimen) names\n"))
    }
    indNames(x) <- ind.list
    
    # If there are populations to be deleted, then recalculate relevant locus metadata and remove monomorphic loci
    
    if ("delete" %in% indNames(x) | "Delete" %in% indNames(x)) {
        # Remove populations flagged for deletion
        if (verbose >= 2) {
            cat(
                report(
                    "  Deleting individuals/samples flagged for deletion (Flagged 'Delete' or 'delete')\n"
                )
            )
        }
        x <-
            gl.drop.ind(x,
                        ind.list = c("Delete", "delete"),
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
    
    if (verbose >= 2) {
        cat("  Summary of recoded dataset\n")
        # cat(paste(' Original No. of loci:',hold.nLoc,'\n')) cat(paste(' New No. of loci:',nLoc(x),'\n'))
        cat(paste("  Original No. of individuals:", nInd(hold), "\n"))
        cat(paste("    New No. of individuals:", nInd(x), "\n"))
        # cat(paste(' Original No. of populations:', hold.nPop,'\n')) cat(paste(' New No. of populations:', nPop(x),'\n'))
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
