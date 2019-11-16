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
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return saves a glenlight object to csv, returns NULL
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' # SNP data
#'   gl.write.csv(testset.gl, outfile="SNP_1row.csv")
#' # Tag P/A data
#'   gl.write.csv(testset.gs, outfile="PA_1row.csv")

gl.write.csv <- function(x, outfile="outfile.csv", outpath=tempdir(), verbose=NULL) {

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  outfilespec <- file.path(outpath, outfile)
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
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
