#' Calculate a distance matrix for populations defined in an \{adegenet\} genlight object
#'
#' This script calculates various distances between populations based on allele frequencies. The distances are
#' calculated by scripts in the {stats} or {vegan} libraries, with the exception of the pcfixed (percent fixed
#' differences) distance.
#' 
#' The distance measure can be one of "manhattan", "euclidean", "pcfixed", "pa", canberra", "bray", 
#' "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , 
#' "binomial", "chao", "cao", "mahalanobis", "maximum", "binary" or "minkowski". Refer to the documentation for
#' dist stats or vegdist vegan for definitions. 
#' 
#' Distance pcfixed calculates the pair-wise count of fixed allelic differences between populations.
#'
#' @param x -- name of the genlight containing the SNP genotypes [required]
#' @param method -- Specify distance measure [euclidean]
#' @param plot -- if TRUE, display a histogram of the genetic distances, and a whisker plot [TRUE]
#' @param boxplot -- if 'standard', plots a standard box and whisker plot; if 'adjusted',
#' plots a boxplot adjusted for skewed distributions ['standard']
#' @param range -- specifies the range for delimiting outliers [1.5 interquartile ranges]
#' @param binary -- Perform presence/absence standardization before analysis using decostand [FALSE]
#' @param p -- The power of the Minkowski distance (typically a value ranging from 0.25 to infinity) [0.5]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [2]
#' @return An object of class 'dist' giving distances between populations
#' @importFrom stats dist
#' @importFrom vegan vegdist
#' @importFrom reshape2 dcast
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.dist.pop(testset.gl, method="euclidean")

gl.dist.pop <- function(x, method="euclidean", plot=TRUE, boxplot="standard", range=1.5, binary=FALSE, p=NULL, verbose=NULL) {

# TIDY UP FILE SPECS
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# FLAG SCRIPT START
  # set verbosity
  if (is.null(verbose) & !is.null(x@other$verbose)) verbose=x@other$verbose
  if (is.null(verbose)) verbose=2
 
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose >= 1){
    cat("Starting",funname,"\n")
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
  
# FUNCTION SPECIFIC ERROR CHECKING

# DO THE JOB

  veganmethod <- c("bray", 
                   "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , 
                   "binomial", "chao", "cao", "mahalanobis")
  distmethod <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "Chebyshev")

  if (!(method%in%veganmethod || method%in%distmethod || method=="pcfixed")){
    stop("Fatal Error: Specified distance method is not among those available\n  Specify one of",veganmethod, distmethod,"or 'pcfixed'\n")
  }
  hard.min.p <- 0.25

  m <- method
  b <- binary
  pr <- p
  u <- TRUE
  d <- TRUE
  
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
    
    # Calculate distance using vegdist {vegan}
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
    # if (m == "pa"){
    #   dd <- gl.report.pa.pop(x)
    #   dd <- dd[,c(3,4,10)]
    #   tmp <- dd
    #   tmp2 <- dd$pop1
    #   dd$pop1 <- dd$pop2
    #   dd$pop2 <- tmp2
    #   dd <- rbind(dd,tmp)
    #   dd <- dcast(dd, pop1 ~ pop2, value.var="totalpriv")
    #   row.names(dd) <- dd[,1]
    #   dd <- dd[,2:length(dd[1,])]
    #   dd <- dd[order(row.names(dd)),]
    #   dd <- dd[,order(colnames(dd))]
    #   if (verbose >= 2) {
    #     cat("  Calculating total private alleles\n")
    #     cat("  Note: this distance may be non-metric, and so should be considered a dissimilarity measure\n")
    #   }
    # }
    
    dd <- as.dist(dd) 
    
  # Revert to original order  
    ord <- rank(popNames(x))
    mat <- as.matrix(dd)[ord, ord]
    dd <- as.dist(mat)
    
# PLOT
  if (plot){
    
    # Save the prior settings for mfrow, oma, mai and pty, and reassign
    op <- par(mfrow = c(2, 1), oma=c(1,1,1,1), mai=c(0.5,0.5,0.5,0.5),pty="m")
    
    # Set margins for first plot
    par(mai=c(1,0.5,0.5,0.5))
    
    # Plot Box-Whisker plot
    if (all(x@ploidy==2)){
      title <- paste0("SNP data (DArTSeq)\nPopulation ",method," Distances")
    } else {
      title <- paste0("Tag P/A data (SilicoDArT)\nPopulation ",method," Distances")
    }  
    
    if (boxplot == "standard"){
      boxplot(dd, horizontal=TRUE, col='red', range=range, main = title)
      cat("  Standard boxplot, no adjustment for skewness\n")
    } else {
      robustbase::adjbox(dd,
                         horizontal = TRUE,
                         col='red',
                         range=range,
                         main = title)
      cat("  Boxplot adjusted to account for skewness\n")
    }  
    
    # Set margins for second plot
    par(mai=c(0.5,0.5,0,0.5))
    hist(dd, 
         main="", 
         xlab="", 
         border="blue", 
         col="red",
         xlim=c(min(dd),max(dd)),
         breaks=100)
  }  
  
# SUMMARY 
    # Print out some statistics
  if(verbose >= 2){
    cat("\n  Reporting inter-population distances\n")
    cat("  Distance measure:",method,"\n")
    cat("    No. of populations =", nPop(x), "\n")
    cat("    Average no. of individuals per population =", nInd(x)/nPop(x), "\n")
    cat("    No. of loci =", nLoc(x), "\n")
    cat("    Miniumum Distance: ",round(min(dd),2),"\n")
    cat("    Maximum Distance: ",round(max(dd),2),"\n")
    cat("    Average Distance: ",round(mean(dd),3),"\n\n")
  }  
    
# FLAG SCRIPT END

  if (plot){
    # Reset the par options    
    par(op)
  }
    
  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
    return(dd)
}