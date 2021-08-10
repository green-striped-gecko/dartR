#' Convert a genlight object to format suitable for input to Bayescan
#'
#' The output text file contains the snp data and relevant BAyescan command lines to guide input.
#' 

#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default bayescan.txt]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return NULL
#' @export
#' @author Custodian: Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2bayescan(testset.gl)

#' @references Foll M and OE Gaggiotti (2008) A genome scan method to identify selected loci appropriate for both dominant and codominant markers: A Bayesian perspective. Genetics 180: 977-993.

gl2bayescan <- function(x, outfile="bayescan.txt", outpath=tempdir(), verbose=NULL) {

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

  if (verbose >= 2) {cat(paste("Extacting SNP data and creating records for each individual\n"))}
  
# Prepare the data  
  mat <- gl.percent.freq(x, verbose=verbose)
  mat <- mat[order(mat$popn),]
  
# convert to character so it can be used in the for loop
  mat$popn <- as.character(mat$popn)

# Create the bayescan input file  
  if (verbose >= 2) {cat(paste("Writing text input file for Bayescan",outfilespec,"\n"))}
  sink(outfilespec)
  
  cat(paste0("[loci]=",nLoc(x)),"\n\n")
  cat(paste0("[populations]=",nPop(x)),"\n\n")

# creating a counter to be used in the for loop
  pop_id <- 0
# adding  one parenthesis that was missing and using pop names in the for loop 
  # for (i in 1:nPop(x)) {
    for (i in as.character(unique(pop(x)))){
        # counter
         pop_id <- pop_id + 1
        # change the variable used in the loop
        # cat(paste0("[pop]=", i), "\n")
        cat(paste0("[pop]=", pop_id), "\n")
        # change the variable used in the loop
        # popi <- mat[mat$popn == mat$popn[i], ]
        popi <- mat[mat$popn == i, ]
        for (j in 1:length(popi$popn)) {
            cat(j, (2 * popi$nobs[j]), 2, popi$sum[j], (2 * popi$nobs[j] - 
                popi$sum[j]), "\n")
        }
        cat("\n")
    }
  
  sink()

  if (verbose >= 3) {cat(paste("Records written to",outfilespec,"\n"))}

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(NULL)

}


