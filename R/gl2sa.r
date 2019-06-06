#' Convert genlight objects to the format used in the SNPassoc package
#'
#' This function exports a genlight object into a SNPassoc object. See \link[SNPassoc]{setupSNP}
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @export
#' @return Returns an object of class 'snp' to be used with \pkg{SNPassoc}
#' @references Gonz?lez, J.R., Armengol, L., Sol?, X., Guin?, E., Mercader, J.M., Estivill, X. and Moreno, V. (2017). SNPassoc: an R package to perform whole genome association studies. Bioinformatics 23:654-655.
#' @author Bernd Guber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2plink(testset.gl)

# Last amended 3-Feb-19

gl2sa <- function(x, verbose=2){

# TIDY UP FILE SPECS

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
      if (nLoc(x)!=nrow(x@other$loc.metrics)) { stop("The number of rows in the @other$loc.metrics table does not match the number of loci in your genlight object!! Most likely you subset your dataset using the '[ , ]' function of adegenet. This function does not subset the number of loci [you need to subset the loci metrics by hand if you are using this approach].")  }

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

  if (verbose >= 2){ cat("  Writing data to SNPassoc object\n")}
  pop <- gl2gi(x)
  xxx <- pegas::as.loci(pop)[,-1]
  sa <- SNPassoc::setupSNP(data.frame(xxx), 1:ncol(xxx))

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(sa)
}
