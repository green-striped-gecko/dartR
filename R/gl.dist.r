#' Calculate a distance matrix for populations defined in an \{adegenet\} genlight object
#'
#' This script calculates various distances between populations based on allele frequencies. 
#' 
#' The distance measure can be one of "manhattan", "euclidean", "canberra", "bray", 
#' "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , 
#' "binomial", "chao", "cao", "mahalanobis", "maximum", "binary" or "minkowski". Refer to the documentation for
#' dist stats or vegdist vegan for definitions.
#'
#' @param gl -- name of the genlight containing the SNP genotypes [required]
#' @param method -- Specify distance measure [method=euclidean]
#' @param binary -- Perform presence/absence standardization before analysis using decostand [binary=FALSE]
#' @param diag -- Compute and display diagonal
#' @param upper -- Return also upper triangle
#' @param p -- The power of the Minkowski distance
#' @return A matrix of distances between populations (class dist)
#' @import adegenet
#' @importFrom stats dist
#' @importFrom vegan vegdist
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl.dist(testset.gl, method="gower", diag=TRUE)

gl.dist <- function(gl, method="euclidean", binary=FALSE, diag=TRUE, upper=FALSE, p=NULL) {
  x <- gl
  veganmethod <- c("bray", 
                  "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , 
                  "binomial", "chao", "cao", "mahalanobis")
  distmethod <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  m <- method
  b <- binary
  d <- diag
  u <- upper
  pr <- p
  
  # Calculate allele frequencies for each population and locus
    x <- gl.percent.freq(x)
  # Select only pop, locus, frequency columns  
    x <- x[,c(1,2,6)]
  # Convert to a pop x locus matrix
    f <- dcast(x, popn ~ locus, value.var="frequency")
  # Reassign names to the populations, and convert from percentages to proportions 
    row.names(f)=f[,1]
    f <- f[,-c(1)]
    f <- f/100
  # Calculate distance using dist {stat}
    if (method %in% distmethod) {
      d <- dist(f, method=m, diag=d, upper=u, p=pr)
      cat(paste("  Calculating distances: ",method,"\n"))
      cat("  Refer to dist {stats} documentation for algorithm\n\n")
    }
    # Calculate distance using dist {stat}
    if (method %in% veganmethod) {
      d <- vegdist(f, method=m, binary=b, diag=d, upper=u, na.rm=TRUE)
      cat(paste("  Calculating distances: ",method,"\n"))
      cat("    Refer to vegdist {vegan} documentation for algorithm\n\n")
    }
    
    return(d)
}
