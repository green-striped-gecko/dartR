#' Generate percentage allele frequencies by locus and population
#'
#' This is a support script, to take SNP data or SilocoDArT presence/absence data grouped into populations in a genelight or genind object \{adegenet\}
#' and generate a table of allele frequencies for each population and locus
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A matrix with allele frequencies (genelight) or presence/absence frequencies (genind) broken down by population and locus
#' @export
#' @importFrom reshape2 melt
#' @importFrom plyr ddply
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' m <-  gl.percent.freq(testset.gl)
#' m

# Last amended 3-Feb-19

gl.percent.freq<- function(x, verbose=2) {

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
      if (nLoc(x)!=nrow(x@other$loc.metrics)) { stop("The number of rows in the loc.metrics table does not match the number of loci in your genlight object!")  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# DO THE JOB
  
  if (verbose > 0) {
    cat("Starting gl.percent.freq: Calculating allele frequencies for populations\n")
    cat("  This may take some time -- be patient\n")
  }
  
# Calculate the required statistics, to be backward compatible 
  nmissing <- apply(as.matrix(x),2, tapply, pop(x), function(x){sum(is.na(x))})
  nobs <- apply(as.matrix(x),2, tapply, pop(x), function(x){sum(!is.na(x))})
  n <- nmissing + nobs
  sum <- apply(as.matrix(x),2, tapply, pop(x), function(x){sum(x,na.rm=TRUE)})
  f <- apply(as.matrix(x),2, tapply, pop(x), function(x){mean(x, na.rm=TRUE)/2})
  f <- f*100

# Convert the matricies to long format pop, locus, value
  nobs <- melt(nobs, na.rm=FALSE)
  nmissing <- melt(nmissing, na.rm=FALSE)
  n <- melt(n, na.rm=FALSE)
  f <- melt(f, na.rm=FALSE)
  sum <- melt(sum, na.rm=FALSE)
  
  if(nPop(x) == 1) {
    m <- cbind(levels(pop(x)),rownames(sum),sum,nobs,nmissing,f,n)
  } else {
    m <- cbind(sum,nobs[,3],nmissing[,3],f[,3],n[,3])
  }  
  colnames(m) <- c("popn","locus","sum","nobs","nmissing","frequency","n")
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(m)
  
}
