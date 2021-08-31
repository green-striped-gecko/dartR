#' @name gl.assign.pca
#' @title Assign an individual of unknown provenance to population based on PCA
#' @description
#' This script assigns an individual of unknown provenance to one or more target populations based on its proximity to
#' each population defined by a confidence ellipse.
#' 
#' The following process is followed.
#' (a) The space defined by the loci is ordinated to yield a series of orthogonal axes (independent), and the top two
#' dimensions are considered. Populations for which the unknown lies outside the specified confidence limits are no longer
#' considered.
#' (b) The ordination is repeated on the remaining populations to again yield a series of orthogonal (independent) axes, 
#' a necessary condition for combining likelihoods calculated from each axis.
#' (b) A workable subset of dimensions is chosen, normally equal to the number of dimensions
#' with substantive eigenvalues. 
#' (c) The log-likelihood of the value for the unknown on each axis is calculated, weighted by the eigenvalue for that axis, 
#' and summed over all dimensions as an assignment weighting. The assignment weighting is calculated for a point on the boundary of
#' the 95% (or as specified) confidence envelope, and an index is calculated as the ratio of the two. A value less than 1
#' indicates that the unknown is within the confidence envelope, a value greater than 1 indicates that the unknown lies
#' outside the confidence envelope. A smaller value of the index indicates greater conformity of the unknown genotype with
#' the target population.
#' @details
#' There are three considerations to the assignment. First, consider only those populations for which the unknown has no
#' private alleles. Private alleles are an indication that the unknown does not belong to a target population (provided that
#' the sample size is adequate, say >=10). This can be evaluated with gl.assign.pa().
#' 
#' A next step is to consider the PCoA plot for populations where no private alleles have been
#' detected and the position of the unknown in relation to the confidence ellipses as is plotted by this script. Note, this 
#' plot is considering only the top two dimensions of the ordination, and so an unknown lying outside the confidence 
#' ellipse can be interpreted as it lying outside the confidence envelope. However, if the unknown lies inside the 
#' confidence ellipse in two dimensions, then it may still lie outside the confidence envelope. This second step is good for 
#' eliminating populations from consideration, but does not provide confidence in assignment. 
#' 
#' The third step is to consider the assignment probabilities. This approach calculates the squared Generalised Linear 
#' Distance (Mahalanobis distance) of the unknown from the centroid for each population, and calculates the probability 
#' associated with its quantile under the zero truncated normal distribution. This index takes into account position of 
#' the unknown in relation to the confidence envelope in all selected dimensions of the ordination. An index is calculated
#' based on the distance of the unknown from the population centroid, scaled for the distance of the standardized
#' confidence ellipse. An index of less than 1 indicates that the unknown lies within the confidence ellipse for the focal
#' population; an index of greater than 1 indicates that the unknown lies outside the confidence ellipse for the focal
#' population. The smaller the index, the greater the confidence in the assignment.
#' 
#'   Wt is a weighted log-likelihood
#'   CE is the value of the Index on the boundary of the 95 % confidence envelope
#'   Index is the position of the unknown relative to the boundary of the 95 % confidence envelope
#'      An index value less than 1 indicates the unknown resides inside the confidence ellipse for the focal population
#'      An index value greater than 1 indicates the unknown resides outside the confidence ellipse for the focal population
#' 
#' Each of these approaches provides evidence, none are 100% definitive. They need to be interpreted cautiously.
#'
#' @param x -- name of the input genlight object [required]
#' @param unknown -- identity label of the focal individual whose provenance is unknown [required]
#' @param alpha -- probability level for bounding ellipses in the PCoA plot [default 0.05]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' 
#' @return A genlight object containing the focal individual (assigned to population "unknown") and #' populations for which the focal individual is not distinctive (number of loci with private alleles less than or equal to thresold t.
#' 
#' @importFrom stats dnorm qnorm
#' @export
#' 
#' @author Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' # Test run with a focal individual from the Macleay River (EmmacMaclGeor)
#'   x <- gl.assign.pa(testset.gl, unknown="UC_00146", nmin=10, threshold=1,verbose=3)
#'   x <- gl.assign.pca(x, unknown="UC_00146", alpha=0.05,verbose=3)


gl.assign.pca <- function (x, unknown, alpha= 0.05, verbose=3) {

# SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
# FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
# CHECK PACKAGES
  pkg <- "SIBER"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error("Package",pkg," needed for this function to work. Please install it."))
  }
  
# CHECK DATATYPE 
  datatype <- utils.check.datatype(x,verbose=0)
  if(nPop(x) < 2){
    stop(error("Fatal Error: Only one population, including the unknown, no putative source"))
  }

# FUNCTION SPECIFIC ERROR CHECKING

  alpha <- 1-alpha

  if (!(unknown %in% indNames(x))) {
    stop(error("Fatal Error: Unknown must be listed among the individuals in the genlight object!\n"))
  }
  if (alpha > 1 || alpha < 0) {
    cat(warn("  Warning: Value of alpha must be between 0 and 1, set to 0.05\n"));
    alpha <- 0.05
    alpha <- 1- alpha
  }
  
  if(nLoc(x) < nPop(x)){
    stop(error("Fatal Error: Number of loci less than number of populations"))
  }

# DO THE JOB
  vec <- as.vector(pop(x))
  vec[indNames(x)==unknown] <- "unknown"
  pop(x) <- as.factor(vec)

# Ordinate a reduced space of 2 dimensions
  pcoa <- gl.pcoa(x, nfactors=2,verbose=0)
# Plot  
  suppressWarnings(suppressMessages(gl.pcoa.plot(pcoa,x,ellipse=TRUE,plevel=alpha))) # Because the unknown pop throws an ellipse error
# Combine Pop names and pca scores
  df <- data.frame(pcoa$scores)
  df <- cbind(as.character(pop(x)),df)
  names(df) <- c("pop","x","y")
# Determine if the unknown lies within the confidence ellipses specified by alpha
  result <- data.frame() 
  count <- 0
  for (i in popNames(x)){
    if(i=="unknown" | nInd(x[pop=i])<=1){next}
    count <- count + 1
    A <- pcoa$scores[df$pop==i,]
    mu <- colMeans(A)
    sigma <- stats::cov(A)
    testset <- rbind(pcoa$scores[df$pop=="unknown",],A)
    transform <- SIBER::pointsToEllipsoid(testset, sigma, mu)
    inside.or.out <- SIBER::ellipseInOut(transform, p = alpha)
    result[count,1] <- i
    result[count,2] <- inside.or.out[1]
  }
  names(result) <- c("pop","hit")
  nhits <- length(result$pop[result$hit])
  nohits <- length(result$pop[!result$hit])
  if(verbose==3){
    if(nhits > 0){
      cat("  Putative source populations:",paste(result[result$hit==TRUE,"pop"],collapse=", "),"\n")
    } else {
      cat("  No putative source populations identified\n")
    } 
    if(all(result$hit ==TRUE)){
      cat("  No populations elimated from consideration\n")
    } else {
      cat("  Populations elimated from consideration:",paste(result[result$hit==FALSE,"pop"],collapse=", "),"\n")
    }
  }
  if(all(result$hit==TRUE)){
    x2 <- x
  } else {
    x2 <- gl.drop.pop(x,pop.list=result[result$hit==FALSE,"pop"],verbose=0)
  }
  
# Re-run the pcoa on the reduced set
  hard.limit <- 8
  pcoa <- gl.pcoa(x2,nfactors=hard.limit,verbose=0)

# Determine the number of dimensions for confidence envelope (the ordination and dimension reduction)
  # From the eigenvalue distribution
    s <- sum(pcoa$eig)
    e <- round(pcoa$eig*100/s,1)
    e <- e[e>mean(e)]
    first.est <- length(e)
  # From the number of populations, including the unknown
  #  sec.est <- nPop(x2)

    #cat("  Number of populations, including the unknown:",sec.est,"\n")
    cat(report("  Number of dimensions with substantial eigenvalues:",first.est,". Hardwired limit",hard.limit,"\n"))
    cat(report("    Selecting the smallest of the two\n"))
    dim <- min(first.est, hard.limit)
    cat(report("    Dimension of confidence envelope set at",dim,"\n"))
    pcoa$scores <- pcoa$scores[,1:dim]
    
# Add population names to the scores   
    c <- data.frame(cbind(pcoa$scores,as.character(pop(x2))),stringsAsFactors = FALSE)
    colnames(c)[dim+1]<- "pop"
    
# Create a set of data without the unknown
  clouds <- c[c[,"pop"]!="unknown",]
  Unknown <- c[c[,"pop"]=="unknown",]
  #Unknown[1:dim] <- as.numeric(Unknown[1:dim])


  cat("  Likelihood Index for assignment of unknown",unknown,"to putative source populations\n") 
# For each population
  p <- as.factor(unique(clouds[,"pop"]))
  for (i in 1:length(levels(p))) {
    # Pull out population i
      m <- clouds[clouds[,"pop"]==levels(p)[i],]
    # Discard the population labels
      m <- m[,1:dim]
      row.names(m) <- NULL
      Unknown <- Unknown[1:dim]
      row.names(Unknown) <- NULL
    # Convert to numeric, and reformat as a matrix
      n <- as.matrix(sapply(m, as.numeric))
      Unknown <- as.numeric(Unknown)
    # Calculate statistics for each dimension
      means <- colMeans(n)
      N <- nrow(n)
      sd <- sqrt(N/(N-1) * (colMeans(n*n)-colMeans(n)^2))
    # Standardise the unknown
      std.unknown <- (Unknown-means)/sd
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
    # Calculate an index
      idx <- round(index/index.boundary,2)
    # Print out the result for population i
      assign <- "no"
      if (idx <= 1) {assign <- "yes"}
      tmp <- data.frame(levels(p)[i],index,index.boundary,idx,assign)
      if (i==1) {
        df <- tmp
      } else {
        df <- rbind(df,tmp)
      }
  }  
  colnames(df) <- c("Population","Wt","CE","Index","Assign")
  df <- df[order(df$Index),]
  print(df)
  best <- as.character(df$Population[df$Assign=="yes"][1])
  # cat("  Wt is a weighted log-likelihood\n")
  # cat("  CE is the value of the Index on the boundary of the",alpha*100,"% confidence envelope\n")
  # cat("  Index is the position of the unknown relative to the boundary of the",alpha*100,"% confidence envelope\n")
  # cat("    An index value less than 1 indicates the unknown resides inside the confidence ellipse for the focal population\n")
  # cat("    An index value greater than 1 indicates the unknown resides outside the confidence ellipse for the focal population\n")
  cat(report("  Best assignment is the population with the smallest value of the Index, in this case",best,"\n"))
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat(report("Completed:",funname,"\n"))
  }

  return(df)
}