#' Calculate likelihood of population assignment
#'
#' This script calculates the likelihood that a focal individual of unknown provenance belongs to one of several
#' target populations.
#' 
#' The algorithm first identifies those target populations for which the individual has no private alleles, then for each
#' population, ordinates the locus space for a reduced representation of the variation among individuals, including the 
#' unknown individual. The mahanobalis distance is calculated for each individual in the target population from the centroid 
#' (without the unknown) and for the focal individual.The likelihood of population membership is then calculated for the 
#' focal individual. This is repeated for all target populations and the results presented as likelihoods and 
#' posterior probabilities.
#'
#' @param gl -- name of the input genlight object [required]
#' @param id -- identity label of the focal individual whose provenance is unknown
#' @param nmin -- minimum sample size for a target population to be included in the analysis
#' @param t -- populations to retain for consideration; those for which the focal individual has less than or equal to t loci with private alleles
#' @return A genlight object containing the focal individual (assigned to population "unknown") and 
#' populations for which the focal individual is not distinctive (number of loci with private alleles less than or equal to thresold t.
#' @import msm
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' x <- testset.gl
#' # Test run with a focal individual from the Macleay River (EmmacMaclGeor)
#' x <- gl.assign(testset.gl, id="UC_00146", nmin=10, t=1)
#' 

gl.assign <- function (gl, id, nmin=10, t=0) {
x <- gl

# Identify populations that can be eliminated on the basis of private alleles
# Retain the remainder for analysis
  x2 <- gl.report.pa(x, id=id, nmin=nmin, t=t)
  x2 <- gl.filter.monomorphs(x2)
  
# Check that there is more than one population to assign (excluding 'unknown')
  if (nPop(x2)==1) {
    cat("There are no populations retained for assignment. Rerun with a higher threshold\n")
    stop("Terminating execution")
  }
  if (nPop(x2)==2) {
    cat("There are no further populations to compare for assignment.",levels(pop(x2))[1],"is the best assignment\n")
    stop("Terminating execution")
  }
  
# Ordinate a reduced space of K = nPop(x2) dimensions
  pcoa <- gl.pcoa(x2, nfactors=nPop(x2))
  print(gl.pcoa.plot(pcoa,x2, xaxis=1, yaxis=2, ellipse=TRUE))
#  gl.pcoa.plot(pcoa,x2, xaxis=3, yaxis=4, ellipse=TRUE)
  
# Add population names to the scores   
  c <- cbind(pcoa$scores,as.character(pop(x2)))
  colnames(c)[nPop(x2)+1]<- "pop"

# Create a set of data without the unknown
  clouds <- c[c[,"pop"]!="unknown",]
  unknown <- c[c[,"pop"]=="unknown",]
  unknown <- as.numeric(unknown[1:nPop(x2)])
  
  cat("\nProbability of assignment of unknown",id,"to listed populations\n\n") 

# For each population
  p <- as.factor(clouds[,"pop"])
  for (i in 1:length(levels(p))) {
    # Pull out population i
      m <- clouds[clouds[,"pop"]==levels(p)[i],]
    # Discard the population labels
      m <- m[,1:nPop(x2)]
    # Convert to numeric, and reformat as a matrix
      n <- as.numeric(m)
      n <- matrix(data=n, ncol=ncol(m), nrow=nrow(m))
    # Calculate the squared M distances for each of the known individuals in population i  
      m.dist <- mahalanobis(n, colMeans(n), cov(n))
      m.dist <- round(m.dist,4)
    # Calculate the mean and standard deviation
      mean.m <- mean(m.dist)
      mean.m <- round(mean.m,4)
      sd.m <- sqrt(var(m.dist))
      sd.m <- round(sd.m,4)
    # Calculate the squared M distance for the unknown  
      m.dist.unk <- mahalanobis(unknown, colMeans(n), cov(n))
      m.dist.unk <- round(m.dist.unk,4)
    # Determine the probability of the unknown (or one more deviant) using the truncated normal density function  
      prob <- dtnorm(m.dist.unk, mean=mean.m, sd=sd.m, lower=0)
      prob <- round(prob,4)
    # Print out the result for population i
      tmp <- data.frame(levels(p)[i],mean.m,sd.m,m.dist.unk,prob)
      if (i==1) {
        df <- tmp
      } else {
        df <- rbind(df,tmp)
      }
  }  
  colnames(df) <- c("Population","Mean M","SD M","Unknown","Prob")
  df[rev(order(df$Prob)),]
  print(df)
  return(df)
}


