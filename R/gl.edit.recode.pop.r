#' @name gl.edit.recode.pop
# Preliminaries -- Parameter specifications -------------- 
#' @title Creates or edits and applies a population re-assignment table
#' @description
#' A function to edit population assignments in a dartR genlight object, or to
#' create a reassignment table taking the population assignments
#' from a genlight object, or to edit existing population assignments in
#' a pop.recode.table. The amended recode table is then applied to the genlight
#' object.
#' @details
#' Genlight objects assign specimens to populations based on information in the
#' ind.metadata file provided when the genlight object is first generated.
#' Often one wishes to subset the data by deleting populations or to amalgamate
#' populations. This can be done with a pop.recode table with two columns. The
#' first column is the population assignment in the genlight object, the second
#' column provides the new assignment.
#'
#' This function will input an existing reassignment table for editing and
#' optionally save it as a new table, or if the name of an input table is not
#' supplied, will generate a table using the population assignments in the
#' parent genlight object. It will then apply the recodings to the genlight 
#' object.
#' 
#' When caution needs to be exercised because of the potential for breaking the
#' 'chain of evidence' associated with the samples, recoding individuals using
#' a recode table (csv) can provide a durable record of the changes.
#'
#' For SNP genotype data, the function, having deleted populations, optionally 
#' identifies resultant monomorphic loci or loci with all values missing 
#' and deletes them. The script also optionally recalculates the
#' locus metadata as appropriate. The optional deletion of monomorphic loci
#' and the optional recalculation of locus statistics is not available for
#' Tag P/A data (SilicoDArT).
#'
#' Use outpath=getwd() when calling this function to direct
#' output files to your working directory.
#'
#' The function returns a dartR genlight object with the new population assignments  
#' and the recalculated locus metadata. 
#'
#' @param x Name of the genlight object  [required].
#' @param pop.recode Path to recode file [default NULL].
#' @param out.recode.file Name of the file to output the new individual labels
#' [default NULL].
#' @param outpath Path where to save the output file [default tempdir(), mandated by CRAN].
#' @param recalc If TRUE, recalculate the locus metadata statistics
#' [default TRUE].
#' @param mono.rm If TRUE, remove monomorphic loci [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress but not results; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#'  
#' @import utils
#' @export
#' @return A genlight object with the revised population assignments
#' 
#' @family dartR-base
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
# Examples --------------
#' @examples
#' \dontrun{
#' gl <- gl.edit.recode.pop(testset.gl)
#' gs <- gl.edit.recode.pop(testset.gs)
#' }
#' # See also -------------------
#' @seealso \code{\link{gl.recode.pop}}, \code{\link{gl.drop.pop}},
#' \code{\link{gl.keep.pop}}, \code{\link{gl.merge.pop}},
#' \code{\link{gl.reassign.pop}}

# Function -----------
gl.edit.recode.pop <-  function(x,
                                pop.recode = NULL,
                                out.recode.file = NULL,
                                outpath = tempdir(),
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
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    outfilespec <- file.path(outpath, out.recode.file)
    
    # Function-specific error checking -----------    
    
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
        stop(error("Fatal Error: Population names not detected\n"))
    }
    
    # DO THE JOB --------------

    # Store variables
    hold <- x
    
    # Take assignments from x
    
    if (verbose >= 2) {
        cat(report("  Extracting current pop assignments from the x object\n"))
    }
    recode.table <- cbind(levels(pop(x)), levels(pop(x)))
    
    # Create recode table for editting, and bring up the editor
    new <- as.matrix(edit(recode.table))
    # new <- new[,1:2]
    
    # Write out the recode table, if requested
    if (is.null(pop.recode)) {
        cat("  No output table specified, recode table not written to disk\n")
    } else {
        if (verbose >= 2) {
            cat(report(
                paste(
                    "  Writing population recode table to: ",
                    pop.recode,
                    "\n"
                )
            ))
        }
        write.table(
            new,
            file = pop.recode,
            sep = ",",
            row.names = FALSE,
            col.names = FALSE
        )
    }
    
    # Apply the new assignments
    pop.list <- as.character(pop(x))
    ntr <- length(new[, 1])
    for (i in 1:nInd(x)) {
        for (j in 1:ntr) {
            if (pop.list[i] == new[j, 1]) {
                pop.list[i] <- new[j, 2]
            }
        }
    }
    # Assigning new populations to x
    if (verbose >= 2) {
        cat(report("  Assigning new population names\n"))
    }
    pop(x) <- pop.list
    
    # If there are populations to be deleted, then recalculate relevant locus 
    # metadata and remove monomorphic loci
    
    if ("delete" %in% x$pop | "Delete" %in% x$pop) {
        # Remove populations flagged for deletion
        if (verbose >= 2) {
            cat(
                report(
                    "  Deleting individuals/samples flagged for deletion (Flagged 'Delete' or 'delete')\n"
                )
            )
        }
        x <-
            gl.drop.pop(x,
                        pop.list = c("Delete", "delete"),
                        verbose = 0)
    }
    
    # Remove monomorphic loci -----------------
    if(datatype=="SNP"){
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
    }
    
    # Recalculate statistics ---------------
    if(datatype=="SNP"){
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
    }
    
    # REPORT A SUMMARY ---------------
    
    if (verbose >= 2) {
        cat("  Summary of recoded dataset\n")
        # cat(paste(' Original No. of loci:',nLoc(hold),'\n')) cat(paste(' New No. of loci:',nLoc(x),'\n')) cat(paste(' Original No. of
        # individuals:', nInd(hold),'\n')) cat(paste(' New No. of individuals:', nInd(x),'\n'))
        cat(paste("  Original No. of populations:", nPop(hold), "\n"))
        cat(paste("    New No. of populations:", nPop(x), "\n"))
        if (!recalc) {
            cat(report("  Note: Locus metrics not recalculated\n"))
        }
        if (!mono.rm) {
            cat(report("  Note: Resultant monomorphic loci not deleted\n"))
        }
    }
    
    # ADD TO HISTORY --------------
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END ---------------
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    # End Block --------------
    
    return(x)
}
