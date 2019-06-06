#' Convert a genlight object to format suitable for input to genalex
#'
#' The output csv file contains the snp data and other relevant lines suitable for genalex. This script is a wrapper for genind2genalex {poppr}
#' 
#' Reference: Peakall, R. and Smouse P.E. (2012) GenAlEx 6.5: genetic analysis in Excel. Population genetic software for teaching and research-an update. Bioinformatics 28, 2537-2539.
#' http://bioinformatics.oxfordjournals.org/content/28/19/2537
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default 'genalex.csv']
#' @param outpath -- path where to save the output file [default tempdir()]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @importFrom poppr genind2genalex
#' @return NULL
#' @export
#' @author Katrin Hohwieler, wrapper Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' gl2genalex(testset.gl, outfile="testset.csv")
#' }


gl2genalex <- function(x, outfile="genalex.csv", outpath=tempdir(), verbose=2) {
  
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

  # Work around a bug in adegenet if genlight object is created by subsetting
      if (nLoc(x)!=nrow(x@other$loc.metrics)) { stop("The number of rows in the loc.metrics table does not match the number of loci in your genlight object!")  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nLoc(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# DO THE JOB
  
  gind <- gl2gi(x, verbose=0)
  poppr::genind2genalex(gind, filename = outfilespec, sequence = TRUE, overwrite = FALSE)

  if (verbose > 2) {cat(paste("    Records written to",outfile,":",nInd(x),"\n"))}
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(NULL)

}


