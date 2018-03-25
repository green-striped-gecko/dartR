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
#' @param t -- populations to retain for consideration; those for which the focal individual has less than or equal to t loci with private alleles [default 0]
#' @param v -- verbosity: 0, silent or errors only; 1, begin and end; 2, progress log; 3, progress and results; 5, full report [default 2]
#' @return A genlight object containing the focal individual (assigned to population "unknown") and #' populations for which the focal individual is not distinctive (number of loci with private alleles less than or equal to thresold t.
#' @importFrom stats dnorm qnorm
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' x <- testset.gl
#' # Test run with a focal individual from the Macleay River (EmmacMaclGeor)
#' x <- gl.assign(testset.gl, id="UC_00146", nmin=10, alpha=0.05, t=1)
#' 

gl.assign <- function (x, id, nmin=10, dim=NULL, alpha= 0.01, t=0, v=2) {
  
  if (v >= 1) {
    cat("Starting gl.assign\n")
  }
  
  alpha <- 1-alpha

  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.recode.pop.r!\n"); stop()
  }
  if (!(id %in% indNames(x))) {
    cat("Fatal Error: Unknown must be listed among the genotypes in the genlight object!\n"); stop()
  }
  if (alpha >1 || alpha <0) {
    cat("Fatal Error: Value of alpha must be between 0 and 1, typically 0.05!\n"); stop()
  }
  
# Identify populations that can be eliminated on the basis of private alleles
# Retain the remainder for analysis
  x2 <- gl.report.pa(x, id=id, nmin=nmin, t=t)

# Check that there is more than one population to assign (excluding 'unknown')
  if (nPop(x2)==1) {
    cat("There are no populations retained for assignment. Conclude that the unknown does not belong to one of the target populations or rerun with a higher threshold\n")
    stop("Terminating execution")
  }
  if (nPop(x2)==2) {
    cat("There are no further populations to compare for assignment.",levels(pop(x2))[1],"is the best assignment\n")
    stop("Terminating execution")
  }
  cat("\n\nCOMPUTING ASSIGNMENT BASED ON CONFIDENCE ENVELOPES\n\n")
# Ordinate a reduced space of K = nPop(x2) dimensions
  pcoa <- gl.pcoa(x2, nfactors=nPop(x2),v=FALSE)

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
    cat("  User specified dimensions to retain:",dim,"\n")
    dim <- min(first.est, sec.est, third.est, dim)
    cat("Dimension of confidence envelope set at",dim,"\n")

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
  cat("  Best assignment is the population with the largest value of the Index, in this case",best,"\n")
  
  if (v >= 1) {
    cat("Completed gl.assign\n\n")
  }
  return(df)
}
