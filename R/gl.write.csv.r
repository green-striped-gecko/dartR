#' Write out data from a gl object \link{adegenet} to csv file
#'
#' This script writes to file the SNP genotypes with specimens as entities (columns) and
#' loci as attributes (rows). Each row has associated locus metadata. Each column, with header
#' of specimen id, has population in the first row. 
#' 
#' The data coding differs from the DArT 1row format in that 0 = reference homozygous, 2 =
#' alternate homozygous, 1 = heterozyous, and NA = missing SNP assignment. 
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default outfile.csv]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return saves a glenlight object to csv, returns NULL
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.write.csv(testset.gl, outfile="SNP_1row.csv")

gl.write.csv <- function(x, outfile="outfile.csv", outpath=tempdir(), verbose=2) {

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

# Add individual names and population names to rows
  x1 <- cbind.data.frame(pop(x),as.matrix(x))
# Transpose to have id as columns, loci as rows  
  x1 <- t(x1)
# Create two filler rows to bring number of data rows and number of locus.metadata rows together
  filler1 <- rep("*",length(x@other$loc.metrics[1,]))
# Bind the filler rows to the locus metadata
  x2 <- rbind(filler1, as.matrix(x@other$loc.metrics))
# Bind the locus metadata to the data  
  x3 <- cbind.data.frame(x2,x1)
# Output
  if (verbose >= 2){cat("  Writing records to",outfilespec,"\n")}
  write.table(x3, file=outfilespec,sep=",",row.names=FALSE)

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(NULL)
}
