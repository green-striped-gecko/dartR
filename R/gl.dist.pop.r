#' Calculate a distance matrix for populations defined in an \{adegenet\} genlight object
#'
#' This script calculates various distances between populations based on allele frequencies. The distances are
#' calculated by scripts in the {stats} or {vegan} libraries, with the exception of the pcfixed (percent fixed
#' differences) and pa (total private alleles) distances.
#' 
#' The distance measure can be one of "manhattan", "euclidean", "pcfixed", "pa", canberra", "bray", 
#' "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , 
#' "binomial", "chao", "cao", "mahalanobis", "maximum", "binary" or "minkowski". Refer to the documentation for
#' dist stats or vegdist vegan for definitions. 
#' 
#' Distance pcfixed calculates the pair-wise count of fixed allelic differences between populations. Distance pa
#' tallies the total number of private alleles possessed by one or the other population.
#'
#' @param x -- name of the genlight containing the SNP genotypes [required]
#' @param method -- Specify distance measure [method=euclidean]
#' @param binary -- Perform presence/absence standardization before analysis using decostand [binary=FALSE]
#' @param diag -- Compute and display diagonal [FALSE]
#' @param upper -- Return also upper triangle [FALSE]
#' @param p -- The power of the Minkowski distance (typically a value ranging from 0.25 to infinity) [0.5]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A matrix of distances between populations (class dist)
#' @importFrom stats dist
#' @importFrom vegan vegdist
#' @importFrom reshape2 dcast
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.dist.pop(testset.gl, method="euclidean", diag=TRUE)

# Last amended 3-Feb-19

gl.dist.pop <- function(x, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE, p=NULL, verbose=2) {

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

# DO THE JOB

  veganmethod <- c("bray", 
                   "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , 
                   "binomial", "chao", "cao", "mahalanobis")
  distmethod <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "Chebyshev")

  if (!(method%in%veganmethod || method%in%distmethod || method=="pcfixed" || method == "pa")){
    cat("Fatal Error: Specified distance method is not among those available\n")
    cat("  Specify one of",veganmethod, distmethod,"\n")
    stop("Terminating Execution")
  }
  hard.min.p <- 0.25

  m <- method
  b <- binary
  d <- diag
  u <- upper
  pr <- p
  
  # Calculate allele frequencies for each population and locus
    f <- gl.percent.freq(x,verbose=0)
  # Select only pop, locus, frequency columns  
    f <- f[,c(1,2,6)]
  # Convert to a pop x locus matrix
    f <- reshape2::dcast(f, popn ~ locus, value.var="frequency")
  # Reassign names to the populations, and convert from percentages to proportions 
    row.names(f)=f[,1]
    f <- f[,-c(1)]
    f <- f/100
  # Calculate distance using dist {stat}
    if (m %in% distmethod) {
      if (verbose >= 2) {
        cat(paste("  Calculating distances: ",m,"\n"))
        cat("  Refer to dist {stats} documentation for algorithm\n")
      }  
      if (method == "minkowski"){
        if (pr < 0.25) {
          if (verbose >= 0){cat("  Warning:",hard.min.p,"is the practical minimum for Minkowski distance, set to,",hard.min.p,"\n\n")}
          pr <- hard.min.p
        }
        if (pr == 1) {
          if (verbose >= 2) {cat("  Note: for p = 1, Minkowski distance is equivalent to Manhattan distance\n\n")}
        }
        if (pr == 2) {
          if (verbose >= 2) {cat("  Note: for p = 2, Minkowski distance is equivalent to Euclidean distance\n\n")}
        }
        if (pr >= 30) {
          if (verbose >= 2) {cat("  Note: for large p, Minkowski distance is equivalent to the Maxiumum Metric distance\n\n")}
        }
        if (pr < 1) {
          if (verbose >= 2) {cat("  Note: for p < 1, Minkowski distance is not a metric distance, and so should be considered a measure of dissimilarity\n\n")}
        }
      }
      dd <- stats::dist(f, method=m, diag=d, upper=u, p=pr)
    }
    # Calculate distance using dist {stat}
    if (m %in% veganmethod) {
      dd <- vegan::vegdist(f, method=m, binary=b, diag=d, upper=u, na.rm=TRUE)
      if (verbose >= 2) {
        cat(paste("  Calculating distances: ",m,"\n"))
        cat("    Refer to vegdist {vegan} documentation for algorithm\n")
      }
      if (method == "bray"){
        if (verbose >= 2) {cat("  Note: the Bray-Curtis distance is non-metric, and so should be considered a dissimilarity measure. A metric alternative is the Jaccard distance.\n\n")}
      }
    }
    if (m == "pcfixed"){
      dd <- gl.fixed.diff(x,verbose=0)[[3]]
      if (verbose >= 2) {
        cat("  Calculating percent fixed differences\n")
        cat("Note: this distance may be non-metric, and so should be considered a dissimilarity measure\n")
      }  
    }
    if (m == "pa"){
      dd <- gl.report.pa.pop(x)
      dd <- dd[,c(3,4,10)]
      tmp <- dd
      tmp2 <- dd$pop1
      dd$pop1 <- dd$pop2
      dd$pop2 <- tmp2
      dd <- rbind(dd,tmp)
      dd <- dcast(dd, pop1 ~ pop2, value.var="totalpriv")
      row.names(dd) <- dd[,1]
      dd <- dd[,2:length(dd[1,])]
      dd <- dd[order(row.names(dd)),]
      dd <- dd[,order(colnames(dd))]
      if (verbose >= 2) {
        cat("  Calculating total private alleles\n")
        cat("  Note: this distance may be non-metric, and so should be considered a dissimilarity measure\n")
      }
    }
    
    dd <- as.dist(dd) 
     
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
    return(dd)
}
