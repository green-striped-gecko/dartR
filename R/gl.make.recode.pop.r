#' Create a proforma recode_pop_table file for reassigning population names
#'
#' Renaming populations may be required when there have been errors in assignment arising
#' in the process from sample to DArT files or when one wishes to amalgamate populations, or delete populations.
#' Recoding populations can also be done with a recode table (csv).
#'
#' This script facilitates 
#' the construction of a recode table by producing a proforma file with
#' current population names in two identical columns. Edit the second
#' column to reassign populations. Use keyword Delete to delete a population.
#' 
#' Apply the recoding using gl.recode.pop(). 
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default recode_pop_table.csv]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A vector containing the new population names
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \donttest{
#' result <- gl.make.recode.pop(testset.gl, outfile="Emmac_recode_pop.csv")
#' }

# Last amended 3-Feb-19

 gl.make.recode.pop <- function(x, outfile="recode_pop_table.csv", outpath=tempdir(), verbose=2) {

# TIDY UP FILE SPECS

  outfilespec <- file.path(outpath, outfile)
  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(!is(x, "genlight")) {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }



  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

# DO THE JOB

 if (verbose >= 2) {cat("  Creating draft lookup table\n")}
   mat <- cbind(levels(pop(x)),levels(pop(x)))
 if (verbose >= 2) {cat("  Writing draft lookup table to",outfilespec,". Edit before use\n")}
   write.table(mat, file=outfilespec, sep=",", row.names=FALSE, col.names=FALSE)

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
 
 return(levels(pop(x)))
 
 }
