#' Filters individuals with average heterozygosity greater than a specified threshold 
#' 
#' Calculates the observed and expected heterozygisities for each population (method="pop") 
#' or the observed heterozyosity for each individual (method="ind") in a genlight object.
#' 
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param threshold -- filter individuals above the threshold [default 0.5]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return the filtered genlight object
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom plyr join
#' @importFrom pegas heterozygosity
#' @examples
#' tmp <- gl.filter.heterozygosity(testset.gl,threshold=0.06,verbose=3)

# Last amended 3-Feb-19

gl.filter.heterozygosity <- function(x, 
                                     threshold=0.7,
                                     verbose=2) {
  
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
    x@other$loc.metrics <- x@other$loc.metrics[1:nLoc(x),]

  # Set a population if none is specified (such as if the genlight object has been 
  # generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, 
                             individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nLoc(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
      cat("  Warning: genlight object contains monomorphic loci\n")
    }
    
  # DO THE JOB

    # Convert to matrix
    m <- as.matrix(x)
    
    # For each individual determine counts of hets, homs and NAs
    c.na <- array(NA, nInd(x))
    c.hets <- array(NA, nInd(x))
    c.hom0 <- array(NA, nInd(x))
    c.hom2 <- array(NA, nInd(x))
    for (i in 1:nInd(x)){
      c.na[i] <- sum(is.na(m[i,]))/nInd(x)
      c.hets[i] <- sum(m[i,]==1,na.rm=TRUE)/(nInd(x)-c.na[i])
      #c.hom0[i] <- sum(m[i,]==0,na.rm=TRUE)/(nInd(x)-c.na[i])
      #c.hom2[i] <- sum(m[i,]==2,na.rm=TRUE)/(nInd(x)-c.na[i])
    }
    
    x.kept <- x[c.hets <= threshold,]
    x.discarded <- x[c.hets > threshold,]

# REPORT THE RESULTS
    if(verbose >= 2){
      cat("Initial number of individuals:",nInd(x),"\n")
      cat("  Number of outlier individuals:",nInd(x.discarded),"with heterozygosity >=",
          threshold, "\n")
      cat(paste0("    Deleted: ",indNames(x.discarded),"[",as.character(pop(x.discarded)),
                 "],"))
      cat("\n  Number of individuals retained:",nInd(x.kept),"\n\n")
    }
    
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  # Return the result
  return(x.kept) 
}

