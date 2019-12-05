#' Create a proforma recode_ind file for reassigning individual (=specimen) names
#'
#' Renaming individuals may be required when there have been errors in labelling arising
#' in the process from sample to DArT files. There may be occasions where renaming
#' individuals is required for preparation of figures. Caution needs to be exercised
#' because of the potential for breaking the "chain of evidence" between the samples themselves
#' and the analyses. Recoding individuals can be done with a recode table (csv).
#'
#' This script facilitates 
#' the construction of a recode table by producing a proforma file with
#' current individual (=specimen) names in two identical columns. Edit the second
#' column to reassign individual names. Use keyword Delete to delete an individual.
#' 
#' Apply the recoding using gl.recode.ind(). Deleting individuals
#' can potentially generate monomorphic loci or loci with all
#' values missing. Clean this up with gl.filter.monomorphic().
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default default_recode_ind.csv]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A vector containing the new individual names
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' result <- gl.make.recode.ind(testset.gl, outfile="Emmac_recode_ind.csv")

# Last amended 3-Feb-19

 gl.make.recode.ind <- function(x, outfile="default_recode_ind.csv", outpath=tempdir(), verbose=2) {

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

# DO THE JOB

 if (verbose >= 2) {cat("  Creating draft lookup table\n")}
   mat <- cbind(indNames(x),indNames(x))
 if (verbose >= 2) {cat("  Writing draft lookup table to",outfilespec,". Edit before use\n")}
   write.table(mat, file=outfilespec, sep="," , row.names=FALSE, col.names=FALSE)

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
   
 
 return(indNames(x))
 
 }
