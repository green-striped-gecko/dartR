#' @name gl.make.recode.pop
# Preliminaries -- Parameter specifications -------------- 
#' @title Creates a proforma recode_pop_table file for reassigning population
#'  names
#' @description
#' Renaming populations may be required when there have been errors in
#' assignment arising in the process from sample to sequence files or when one
#' wishes to amalgamate populations, or delete populations. Recoding populations
#' can also be done with a recode table (csv).
#' @details
#' This function facilitates the construction of a recode table by producing a
#' proforma file with current population names in two identical columns. Edit
#' the second column to reassign populations. Use keyword 'Delete' to delete a
#' population.
#' 
#' When caution needs to be exercised because of the potential for breaking the
#'  'chain of evidence' associated with the samples, recoding individuals using
#'  a recode table (csv) can provide a clear record of the changes.
#'  
#' Use outpath=getwd() or when calling this function to direct output files 
#' to your working directory.
#' 
#' The function works with both genlight objects
#' containing SNP genotypes and Tag P/A data (SilicoDArT).
#'
#' Apply the recoding using gl.recode.pop(). 
#'
#' @param x Name of the genlight object [required].
#' @param out.recode.file File name of the output file (including extension)
#'  [default recode_pop_table.csv].
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN]. 
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A vector containing the new population names.
#' @export
#' @author Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' result <- gl.make.recode.pop(testset.gl,out.recode.file='test.csv',outpath=tempdir(),verbose=2)
#' 
# Function ---------------
gl.make.recode.pop <- function(x,
                               out.recode.file = "recode_pop_table.csv",
                               outpath = tempdir(),
                               verbose = NULL) {
  # Preliminaries --------------------
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    outfilespec <- file.path(outpath, out.recode.file)
    
    # Function specific error checking ---------------------
    
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
        stop(error("Fatal Error: Population names not detected\n"))
    }
    
    # DO THE JOB --------------------
    
    # if (verbose >= 2) {cat(' Creating draft lookup table\n')}
    mat <- cbind(levels(pop(x)), levels(pop(x)))
    if (verbose >= 2) {
        cat(report(
            "  Writing draft lookup table to",
            outfilespec,
            ". Edit before use\n"
        ))
    }
    write.table(
        mat,
        file = outfilespec,
        sep = ",",
        row.names = FALSE,
        col.names = FALSE
    )
    
    # FLAG SCRIPT END ---------------
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    # Block end ----------------------
    
    return(NULL)
}
