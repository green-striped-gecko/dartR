#' Filters individuals with average heterozygosity greater than a specified upper threshold or less than a specified lower threshold.
#' 
#' Calculates the observed heterozyosity for each individual in a genlight object and filters individuals based on specified threshold values.
#' Use gl.report.heterozygosity to determine the appropriate thresholds.
#' 
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param t.upper -- filter individuals > the threshold [default 0.7]
#' @param t.lower -- filter individuals < the threshold [default 0]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return the filtered genlight object
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom plyr join
#' @importFrom pegas heterozygosity
#' @examples
#' tmp <- gl.filter.heterozygosity(testset.gl,t.upper=0.06,verbose=3)
#' gl.report.heterozygosity(tmp,method="ind")

# Last amended 28-Feb-19

gl.filter.heterozygosity <- function(x, t.upper=0.7, t.lower=0, verbose=2) {
  
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
  
  if (t.upper < 0 | t.upper > 1){
    cat("  Warning: Parameter 't.upper' must lie between 0 and 1, set to the default of 0.7\n")
    t.upper <- 0.7
  }
  
  if (t.lower < 0 | t.lower > 1){
    cat("  Warning: Parameter 't.lower' must lie between 0 and 1, set to the default of 0.005\n")
    t.upper <- 0.005
  }

  # Set a population if none is specified (such as if the genlight object has been 
  # generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, 
                             individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
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
    #c.hom0 <- array(NA, nInd(x))
    #c.hom2 <- array(NA, nInd(x))
    for (i in 1:nInd(x)){
      c.na[i] <- sum(is.na(m[i,]))/nInd(x)
      c.hets[i] <- sum(m[i,]==1,na.rm=TRUE)/(nInd(x)-c.na[i])
      #c.hom0[i] <- sum(m[i,]==0,na.rm=TRUE)/(nInd(x)-c.na[i])
      #c.hom2[i] <- sum(m[i,]==2,na.rm=TRUE)/(nInd(x)-c.na[i])
    }
    
    if (verbose >=2 ) {
      cat("\nRetaining individuals with heterozygosity in the range",t.lower,"to",t.upper,"\n")
    }
    x.kept <- x[c.hets >= t.lower & c.hets <= t.upper,]
    if (any(c.hets > t.upper)) {
      x.discarded.upper <- x[c.hets > t.upper,]
    }  
    if (any(c.hets < t.lower)) {
      x.discarded.lower <- x[c.hets < t.lower,]
    }  

# REPORT THE RESULTS
    if(verbose >= 3){
      cat("Initial number of individuals:",nInd(x),"\n\n")

      if (any(c.hets > t.upper)) {
        cat("  Number of outlier individuals:",nInd(x.discarded.upper),"with heterozygosity >",t.upper, "\n")
        cat("    Deleted: ")
        cat(paste0(indNames(x.discarded.upper),"[",as.character(pop(x.discarded.upper)),"],"))
        cat("\n\n")
      } else {
        if (!(t.upper == 1)) {cat("  No outlying individuals with heterozygosity >",t.upper, "\n")}
      }

      if (any(c.hets < t.lower)) {
        cat("  Number of outlier individuals:",nInd(x.discarded.lower),"with heterozygosity <",t.lower, "\n")
        cat("    Deleted: ")
        cat(paste0(indNames(x.discarded.lower),"[",as.character(pop(x.discarded.lower)),"],"))
        cat("\n\n")
      } else {
        if ((!t.lower == 0)) {cat("  No outlying individuals with heterozygosity <",t.lower, "\n")}
      }

      cat("  Number of individuals retained:",nInd(x.kept),"\n\n")
    }
    
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("\nCompleted:",funname,"\n")
  }
    
  #add to history
    nh <- length(x.kept@other$history)
    x.kept@other$history[[nh + 1]] <- match.call()

  # Return the result
  return(x.kept) 
}

