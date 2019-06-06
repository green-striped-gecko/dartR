#' Assign an individual of unknown provenance to population
#'
#' This script assigns an individual of unknown provenance to one or more target populations based on first, an analysis
#' of private alleles, and then, if the assignment remains ambigous, on the basis of a weighted likelihood index.
#' 
#' The algorithm first identifies those target populations for which the individual has no private alleles. If no single 
#' population emerges from this analysis, or if a higher threshold than 0 is chosen for the number of tollerable private
#' alleles, then the following process is followed.
#' (a) The space defined by the loci is ordinated to yield a series of orthogonal axes (independent), a necessary condition
#' for combining likelihoods calculated from each axis.
#' (b) A workable subset of dimensions is chosen, normally equal to the number of target populations or the number of dimensions
#' with substantive eigenvalues, whichever is the smaller. 
#' (c) The log-likelihood of the value for the unknown on each axis is calculated, weighted by the eigenvalue for that axis, 
#' and summed over all dimensions as an assignment index. The assignment index is calculated for a point on the boundary of
#' the 95% (or as specified) confidence envelope.
#' 
#' There are three considerations to the assignment. First, consider only those populations for which the unknown has no
#' private alleles. Private alleles are an indication that the unknown does not belong to a target population (provided that
#' the sample size is adequate, say >=10). 
#' 
#' Second, consider the PCoA plot for populations where no private alleles have been
#' detected and the position of the unknown in relation to the confidence ellipses. Note, this is considering only the 
#' top two dimensions of the ordination, and so an unknown lying outside the confidence ellipse can be interpreted as it lying 
#' outside the confidence envelope. However, if the unknown lies inside the confidence ellipse in two dimensions, then it may still lie outside
#' the confidence envelope. This is good for eliminating populations from consideration, but does not provide confidence in
#' assignment. 
#' 
#' Third, consider the assignment probabilities. This approach calculates the squared Generalised Linear Distance (Mahalanobis
#' distance) of the unknown from the centroid for each population, and calculates the probability associated with its quantile
#' under the zero truncated normal distribution. This index takes into account position of the unknown in relation to the 
#' confidence envelope in all selected dimensions of the ordination.
#' 
#' Each of these approaches provides evidence, none are 100% definitive. They need to be interpreted cautiously.
#'
#' @param x -- name of the input genlight object [required]
#' @param id -- identity label of the focal individual whose provenance is unknown [required]
#' @param nmin -- minimum sample size for a target population to be included in the analysis [default 10]
#' @param dim -- number of dimensions to retain in the dimension reduction [default k, number of populations]
#' @param alpha -- probability level for bounding ellipses in the PCoA plot [default 0.05]
#' @param threshold -- populations to retain for consideration; those for which the focal individual has less than or equal to threshold loci with private alleles [default 0]
#' @param verbose -- verbosity: 0, silent or errors only; 1, begin and end; 2, progress log; 3, progress and results; 5, full report [default 3]
#' @return A genlight object containing the focal individual (assigned to population "unknown") and #' populations for which the focal individual is not distinctive (number of loci with private alleles less than or equal to thresold t.
#' @importFrom stats dnorm qnorm
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' # Test run with a focal individual from the Macleay River (EmmacMaclGeor)
#' x <- gl.assign(testset.gl, id="UC_00146", nmin=10, alpha=0.05, threshold=1)

# Last amended 3-Feb-19

gl.assign <- function (x, id, nmin=10, dim=NULL, alpha= 0.05, threshold=0, verbose=3) {

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
      pop(x) <- array("pop1",dim = nLoc(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# FUNCTION SPECIFIC ERROR CHECKING

  alpha <- 1-alpha

  if (!(id %in% indNames(x))) {
    cat("Fatal Error: Unknown must be listed among the individuals in the genlight object!\n"); stop()
  }
  if (alpha > 1 || alpha < 0) {
    cat("Warning: Value of alpha must be between 0 and 1, set to 0.05\n");
    alpha <- 0.05
    alpha <- 1- alpha
  }
  if (!is.null(dim)) {
    if(dim < 1) {
      cat("Warning: Value of dim must be greater than zero, set to NULL\n");
      dim <- NULL
    }  
  }
  if (threshold < 0) {
    cat("Warning: Threshold value cannot be negative, set to 0\n")
    threshold <- 0
  }

# DO THE JOB

# Identify populations that can be eliminated on the basis of private alleles
# Retain the remainder for analysis
  x2 <- gl.report.pa(x, id=id, nmin=nmin, threshold=threshold, verbose=verbose)

# Check that there is more than one population to assign (excluding 'unknown')
  if (nPop(x2)==1) {
    cat("There are no populations retained for assignment. Conclude that the unknown does not belong to one of the target populations or rerun with a higher threshold\n")
    stop("Terminating execution")
  }
  if (nPop(x2)==2) {
    return(cat("There are no further populations to compare for assignment.",levels(pop(x2))[1],"is the best assignment\n"))
    #stop("Terminating execution")
  }
  cat("\n\nCOMPUTING ASSIGNMENT BASED ON CONFIDENCE ENVELOPES\n\n")
# Ordinate a reduced space of K = nPop(x2) dimensions
  pcoa <- gl.pcoa(x2, nfactors=nPop(x2),verbose=FALSE)

#  gl.pcoa.plot(pcoa,x2, xaxis=3, yaxis=4, ellipse=TRUE)
  
# Determine the number of dimensions for confidence envelope (the ordination and dimension reduction)
  # From the eigenvalue distribution
    s <- sum(pcoa$eig)
    e <- round(pcoa$eig*100/s,1)
    e <- e[e>mean(e)]
    first.est <- length(e)
  # From the number of populations, including the unknown
    sec.est <- nPop(x2)
  # From a hard set maximum
    third.est <- 8
    cat("  Number of populations, including the unknown:",sec.est,"\n")
    cat("  Number of dimensions with substantial eigenvalues:",first.est,"\n")
    cat("  Hard coded upper limit to dimensions:",third.est,"\n")
    if (is.null(dim)){
      cat("  User specified dimensions to retain: Not specified\n")
    } else {
      cat("  User specified dimensions to retain:",dim,"\n")
    }  
    dim <- min(first.est, sec.est, third.est, dim)
    cat("    Dimension of confidence envelope set at",dim,"\n")

# Re-run the PCoA with calculated dimension
    pcoa <- gl.pcoa(x2, nfactors=dim)

# Plot the PCoA
    suppressMessages(print(gl.pcoa.plot(pcoa,x2, xaxis=1, yaxis=2, ellipse=TRUE, p=alpha)))
# Add population names to the scores   
  c <- cbind(pcoa$scores,as.character(pop(x2)))
  colnames(c)[dim+1]<- "pop"

# Create a set of data without the unknown
  clouds <- c[c[,"pop"]!="unknown",]
  unknown <- c[c[,"pop"]=="unknown",]
  unknown <- as.numeric(unknown[1:dim])
  

  cat("\nLikelihood Index for assignment of unknown",id,"to listed populations\n\n") 
# For each population
  p <- as.factor(clouds[,"pop"])
  for (i in 1:length(levels(p))) {
    # Pull out population i
      m <- clouds[clouds[,"pop"]==levels(p)[i],]
    # Discard the population labels
      m <- m[,1:dim]
    # Convert to numeric, and reformat as a matrix
      n <- as.numeric(m)
      n <- matrix(data=n, ncol=ncol(m), nrow=nrow(m))
    # Calculate statistics for each dimension
      means <- colMeans(n)
      N <- nrow(n)
      sd <- sqrt(N/(N-1) * (colMeans(n*n)-colMeans(n)^2))
    # Standardise the unknown
      std.unknown <- (unknown-means)/sd
      std.boundary <- rep(qnorm((1-(1-alpha)/2), 0, 1),dim)
    # Calculate log likelihoods
      li <- dnorm(std.unknown,0,1)
      li <- log(li)
      li.boundary <- dnorm(std.boundary,0,1)
      li.boundary <- log(li.boundary)
    # Weight by eigenvalues
      li <- li*pcoa$eig[1:dim]
      li.boundary <- li.boundary*pcoa$eig[1:dim]
    # Calculate an assignment index
      index <- round(sum(li),4)
      index.boundary <- round(sum(li.boundary),4)
    # Print out the result for population i
      assign <- "no"
      if (abs(index) <= abs(index.boundary)) {assign <- "yes"}
      tmp <- data.frame(levels(p)[i],index,index.boundary, assign)
      if (i==1) {
        df <- tmp
      } else {
        df <- rbind(df,tmp)
      }
  }  
  colnames(df) <- c("Population","Index","CE","Assign")
  df <- df[rev(order(df$Index)),]
  print(df)
  best <- as.character(df$Population[df$Assign=="yes"][1])
  cat("  Index is a weighted log-likelihood\n")
  cat("  CE is the value of the Index on the boundary of the",alpha*100,"% confidence envelope\n")
  cat("  Best assignment is the population with the largest value of the Index, in this case",best,"\n\n")
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(df)
}
