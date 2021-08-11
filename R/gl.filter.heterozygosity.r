#' @name gl.filter.heterozygosity
#' @title Filters individuals with average heterozygosity greater than a specified upper threshold or less than a specified lower threshold.
#' @description 
#' Calculates the observed heterozyosity for each individual in a genlight object and filters individuals based on specified threshold values.
#' Use gl.report.heterozygosity to determine the appropriate thresholds.
#' 
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param t.upper -- filter individuals > the threshold [default 0.7]
#' @param t.lower -- filter individuals < the threshold [default 0]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2, unless specified using gl.set.verbosity]
#' @return the filtered genlight object
#' @export
#' @author Custodian: Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom plyr join
#' @examples
#'  result <- gl.filter.heterozygosity(testset.gl,t.upper=0.06,verbose=3)
#'  tmp <- gl.report.heterozygosity(result,method="ind")

gl.filter.heterozygosity <- function(x, 
                                     t.upper=0.7, 
                                     t.lower=0, 
                                     verbose=NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(x,accept="SNP",verbose=verbose)
  
# SCRIPT SPECIFIC ERROR CHECKING
  
  if (t.upper < 0 | t.upper > 1){
    stop(error("Fatal Error:Parameter 't.upper' must lie between 0 and 1\n"))
  }
  
  if (t.lower < 0 | t.lower > 1){
    stop(error("Fatal Error:Parameter 't.upper' must lie between 0 and 1\n"))
  }
  
  # Check for monomorphic loci
  
  if (!x@other$loc.metrics.flags$monomorphs) {
    cat(warn("  Warning: genlight object contains monomorphic loci which will be factored into heterozygosity estimates\n"))
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
      cat(report("  Retaining individuals with heterozygosity in the range",t.lower,"to",t.upper,"\n"))
      cat(paste("  Minimum individual heterozygosity",round(min(c.hets),4),"\n"))
      cat(paste("  Maximum individual heterozygosity",round(max(c.hets),4),"\n"))

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
      cat("  Initial number of individuals:",nInd(x),"\n")

      if (any(c.hets > t.upper)) {
        cat("  Number of outlier individuals (heterozygosity  >",t.upper,"):",nInd(x.discarded.upper),"\n")
        cat("    Deleted: ")
        cat(paste0(indNames(x.discarded.upper),"[",as.character(pop(x.discarded.upper)),"],"))
        #cat("\n")
      } else {
        if (!(t.upper == 1)) {cat("  Zero outlier individuals with heterozygosity >",t.upper, "\n")}
      }

      if (any(c.hets < t.lower)) {
        cat("  Number of outlier individuals:",nInd(x.discarded.lower),"with heterozygosity <",t.lower, "\n")
        cat("    Deleted: ")
        cat(paste0(indNames(x.discarded.lower),"[",as.character(pop(x.discarded.lower)),"],"))
        cat("\n\n")
      } else {
        if ((!t.lower == 0)) {cat("  No outlying individuals with heterozygosity <",t.lower, "\n")}
      }

      cat(report("  Number of individuals retained:",nInd(x.kept),"\n"))
    }
    
# ADD ACTION TO HISTORY

    nh <- length(x.kept@other$history)
    x.kept@other$history[[nh + 1]] <- match.call()
    
# FLAG SCRIPT END

  if (verbose > 0) {
      cat(report("\nCompleted:",funname,"\n"))
  }

  return(x.kept) 
}

