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
#' @param out.recode.file -- file name of the output file (including extension) [default default_recode_ind.csv]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A vector containing the new individual names
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' result <- gl.make.recode.ind(testset.gl, out.recode.file="Emmac_recode_ind.csv")

 gl.make.recode.ind <- function(x, out.recode.file="default_recode_ind.csv", outpath=tempdir(), verbose=2) {

   # TIDY UP FILE SPECS
   
   funname <- match.call()[[1]]
   build <- "Jacob"
   outfilespec <- file.path(outpath, out.recode.file)
   
   # FLAG SCRIPT START
   
   if (verbose < 0 | verbose > 5){
     cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
     verbose <- 2
   }
   
   if (verbose >= 1){
     cat("Starting",funname,"[ Build =",build,"]\n")
   }
   
   # STANDARD ERROR CHECKING
   
   if(class(x)!="genlight") {
     stop("Fatal Error: genlight object required!\n")
   }
   
   if (all(x@ploidy == 1)){
     if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
     data.type <- "SilicoDArT"
   } else if (all(x@ploidy == 2)){
     if (verbose >= 2){cat("  Processing a SNP dataset\n")}
     data.type <- "SNP"
   } else {
     stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
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
   
 
 return(NULL)
 
 }
 
 # # Test scripts
 # gl.make.recode.ind(gl,out.recode.file="test.csv",outpath=getwd(),verbose=3)
 # gl.make.recode.ind(gs,out.recode.file="test.csv",outpath=getwd(),verbose=3)
