#' Converts a genlight object to PLINK file format
#'
#' This function exports a genlight object into PLINK format and save it into a file
#'
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007).
#' PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics 81:551-575.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default plink.csv]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() when calling this function or set.tempdir <- getwd() elsewhere in your script
#' to direct output files to your working directory.
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @export
#' @author Bernd Guber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2plink(testset.gl)

gl2plink <- function(x, outfile="plink.csv", outpath=tempdir(), verbose=2) {

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
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  
  if (verbose >= 2){
    if (all(x@ploidy == 1)){
      stop("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }
  
# DO THE JOB

  tt <- as.matrix(x)
  if (verbose >= 2){ cat("  Writing data to output file\n")}
  write.csv(tt, file=outfilespec, row.names = TRUE, na = "-9")

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(NULL)

}












