#' Identify putative parent offspring within a population
#'
#' This script examines the frequency of pedigree inconsistent loci, that is,
#' those loci that are homozygotes in the parent for the reference allele, and
#' homozygous in the offspring for the alternate allele. This condition is not
#' consistent with any pedigree, regardless of the (unkonwn) genotype of the other
#' parent. The pedigree inconsistent loci are counted as an indication of whether
#' or not it is reasonable to propose the two individuals are in a parent-offspring
#' relationship.
#' 
#' Obviously, if the two individuals are in a parent offspring relationship, the true
#' number of pedigree inconsistent loci should be zero, but SNP calling is not 
#' infallible. Some loci will be mis-called. The problem thus becomes one of determining
#' if the two focal individuals have a count of pedigree inconsistent loci less than
#' would be expected of typical unrelated individuals. There are some quite sophisticated
#' software packages available to formally apply likelihoods to the decision, but we
#' use a simple outlier comparison.
#' 
#' To reduce the frequency of mis-calls, and so emphasise the difference between true
#' parent-offspring pairs and unrelated pairs, the data can be filtered on read depth.
#' Typically minimum read depth is set to 5x, but you can examine the distribution
#' of read depths with gl.report.rdepth() and push this up with an acceptable loss
#' of loci. 12x might be a good minimum for this particular analysis. It is sensible
#' also to push the minimum reproducibility up to 1, if that does not result in an
#' unacceptable loss of loci.
#' 
#' Note that the null expectation is not well defined, and the power reduced, if the
#' population from which the putative parent-offspring pairs are drawn contains 
#' many sibs. Note also that if an individual has been genotyped twice in the dataset, 
#' the replicate pair will be assessed by this script as being in a parent-offspring 
#' relationship.
#'
#' @param x Name of the genlight object containing the SNP genotypes [required]
#' @param min.rdepth Minimum read depth to include in analysis [default = 12]
#' @param min.reproducibility Minimum reproducibility to include in analysis [default = 1]
#' @param range -- specifies the range to extend beyond the interquartile range for delimiting outliers [default = 1.5 interquartile ranges]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A set of individuals in parent-offspring relationship
#' @importFrom stats median IQR
#' @export
# @rawNamespace import(ggplot2, except = empty)

gl.report.parent.offspring <- function(x,
                                       min.rdepth=12,
                                       min.reproducibility=1,
                                       range=1.5,
                                       verbose=NULL) {
  
# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
  # FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    stop("Fatal Error: Detected Presence/Absence (SilicoDArT) data, require SNP data\n")
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)!")
  }
  
# SCRIPT SPECIFIC ERROR CHECKING
  
# DO THE JOB
  
# Generate null expectation for pedigree inconsistency, and outliers
  if (verbose >= 2){
    cat("  Generating null expectation for distribution of counts of pedigree incompatability\n")
  } 
  # Assign individuals as populations
  # x <- gl.reassign.pop(x,as.pop="id",verbose=0)
  pop(x) <- x$ind.names
  # Filter stringently on reproducibility to minimize miscalls
  x <- gl.filter.reproducibility(x,threshold=min.reproducibility,verbose=0)
  # Filter stringently on read depth, to further minimize miscalls
  x <- gl.filter.rdepth(x,lower=min.rdepth,verbose=0)
  # Filter on call rate to reduce computation time
  x <- gl.filter.callrate(x,threshold = 0.95,plot=FALSE,verbose=0)
  # Preliminaries before for loops
  nL <- nLoc(x)
  nP <- nPop(x)
  pN <- popNames(x)
  count <- array(NA,dim=c(nP,nP))
  row.names(count) <- popNames(x)
  colnames(count) <- popNames(x)
  # For pairs of individuals
  for (i in 1:(nP-1)){
  for (j in (i+1):nP){
    pair <- gl.keep.pop(x,pop.list=c(pN[i],pN[j]),verbose=0)
    mat <- as.matrix(pair)
    vect <- mat[1,]*10+mat[2,]
    homalts <- sum(vect==2 | vect==20, na.rm=T)
    count[i,j] <- homalts
  }
  }
# Remove NAs
  counts <- count[!is.na(count)]

  # Prepare for plotting
  
  if (verbose >= 2){
    cat(  "  Identifying outliers with lower than expected counts of pedigree inconsistencies\n")
  }
    title <- paste0("SNP data (DArTSeq)\nCounts of pedigree incompatable loci per pair")
  # Save the prior settings for mfrow, oma, mai and pty, and reassign
    op <- par(mfrow = c(2, 1), oma=c(1,1,1,1), mai=c(0.5,0.5,0.5,0.5),pty="m")
  # Set margins for first plot
    par(mai=c(1,0.5,0.5,0.5))

  # Plot Box-Whisker plot
    if (verbose >= 2) {cat("  Standard boxplot, no adjustment for skewness\n")} 
    whisker <- graphics::boxplot(counts, horizontal=TRUE, col='steelblue', range=range, main = title)
    lower.extremes <- whisker$out[whisker$out < stats::median(counts)]
    if (length(lower.extremes)==0){
      outliers <- NULL
    } else {
      outliers <- data.frame(Outlier=lower.extremes)
    }

# Ascertain the identity of the pairs
  if (verbose >= 2){
      cat(  "  Identifying outlying pairs\n")
  }
    if(length(lower.extremes)>0){
  tmp <- count
  tmp[lower.tri(tmp)] = t(tmp)[lower.tri(tmp)]
  for (i in 1:length(outliers$Outlier)){
    # Identify
    tmp2 <- tmp[tmp==outliers$Outlier[i]]
    outliers$ind1[i] <- popNames(x)[!is.na(tmp2)][1]
    outliers$ind2[i] <- popNames(x)[!is.na(tmp2)][2]
    # Z-scores
    zscore <- (mean(count,na.rm=TRUE)-outliers$Outlier[i])/sd(count,na.rm=TRUE)
    outliers$zscore[i] <- round(zscore,2)
    outliers$p[i] <- round(pnorm(mean=mean(count,na.rm=TRUE),sd=sd(count,na.rm=TRUE),q=outliers$zscore[i]),4)
  }
    }
# Extract the quantile threshold
  iqr <- stats::IQR(counts,na.rm = TRUE)
  qth <- quantile(counts,0.25,na.rm=TRUE)
  cutoff <- qth - iqr*range

# Set margins for second plot
  par(mai=c(0.5,0.5,0,0.5)) 
  
# Plot Histogram
   hist(counts, xlab="No. Pedigree Incompatable", col='steelblue',breaks=100, main=NULL)
   abline(v=cutoff, col="red")

# Output the outlier loci 
# if (length(whisker$out)==0){
  if(length(lower.extremes)==0){
    
  # if (verbose >= 3){
    cat("  No outliers detected\n")
  # }
} 
   # else {  
   if(length(lower.extremes)>0){  
  outliers <- outliers[order(outliers$Outlier),]
  if (verbose >= 3){
    #cat("  Outliers detected -- \n")
    print(outliers)
  }  
}

# Reset the par options    
  par(op)  
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(outliers)

}