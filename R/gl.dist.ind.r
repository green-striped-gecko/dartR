#' Calculate a distance matrix for individuals defined in an \{adegenet\} genlight object
#'
#' This script calculates various distances between individuals based on allele frequencies. The distances are
#' calculated by scripts in the {stats} or {vegan} libraries, with the exception of the pcfixed (percent fixed
#' differences) distance.
#' 
#' The distance measure for SNP data can be one of 
#' 
#'  Euclidean -- eclidean distance as computed by dist() in {stat}
#'  locus.count -- number of loci for which individuals differ, as implemented by dist.gene() in {ape}
#'  allele.count -- number of allelic differences between two individuals, as implemented by diss.dist() in {poppr}
#'  relatedness -- genetic relatedness between individuals (G matrix), as implemented by A.mat() in {rrBLUP}
#'  
#' The distance measure for Tag P/A data (binary) can be one of
#'  
#'  Simple -- simple matching, both 1 or both 0 = 0; one 1 and the other 0 = 1. Presence and absence equally weighted.
#'  Jaccard -- ignores matching 0, both 1 = 0; one 1 and the other 0 = 1. Absences could be for different reasons.
#'  Dice -- both 0 = 0; both 1 = 2; one 1 and the other 0 = 1. Absences could be for different reasons. Sometimes called the Czekanowski or Sorensen distance.
#'  Phi -- binary analogue of the Pearson Correlation coefficient.
#'  
#' Refer to the documentation in the relevant packages listed above.
#'  
#' @param x -- name of the genlight containing the SNP genotypes [required]
#' @param method -- Specify distance measure [SNP: Euclidean; P/A: Simple]
#' #' @param plot -- if TRUE, display a histogram of the genetic distances, and a whisker plot [TRUE]
#' @param boxplot -- if 'standard', plots a standard box and whisker plot; if 'adjusted',
#' plots a boxplot adjusted for skewed distributions ['standard']
#' @param range -- specifies the range for delimiting outliers [1.5 interquartile ranges]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [2]
#' @return An object of class 'dist' giving distances between individuals
#' @importFrom ape dist.gene
#' @importFrom poppr diss.dist
#' @importFrom stats dist
#' @importFrom rrBLUP A.mat
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.dist.pop(testset.gl, method="euclidean")

gl.dist.ind <- function(x, method=NULL, plot=TRUE, boxplot="standard", range=1.5, verbose=NULL) {

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
  
  if (is.null(method) && data.type=='SNP'){
    method <- 'Euclidean'
  }
  if (is.null(method) && data.type=='SilicoDArT'){
    method <- 'Simple'
  }
  method <- tolower(method)
  
  if (!(method %in% c("euclidean", "locus.count", "allele.count", "relatedness", "simple", "jaccard", "dice", "sorenson", "czekanowski", "phi"))){
    cat(" Warning: Method not in the list of options, set to euclidean for SNP data; simple matching for Tag P/A data\n")
    if (data.type == "SNP"){method <- 'euclidean'}
    if (data.type == "SilicoDArT"){method <- 'simple'}
  }

# DO THE JOB

if (data.type == "SNP"){
  
#  method.list <- c("euclidean", "locus.count", "allele.count", "relatedness")

  # Calculate euclidean distance using dist {adegenet}
    if (method == 'euclidean'){
      dd <- stats::dist(x)
      if (verbose >= 2){
        cat("  Calculating Euclidean Distances between individuals\n")
      }
    }
  # Calculate the number of loci that are different between individuals using dist.gene {ape}
  if (method == 'locus.count'){
    dd <- ape::dist.gene(as.matrix(x), method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
    if (verbose >= 2){
      cat("  Calculating number of loci for which individuals differ\n")
    }
  }
  
  # Calculate the number of allelic differences between individuals using dist.gene {poppr}
  if (method == 'allele.count'){
    dd <- poppr::diss.dist(gl2gi(x), percent = FALSE, mat = FALSE)
    if (verbose >= 2){
      cat("  Calculating number of allelic differences between individuals\n")
    }
  }
    
  # Calculate the genetic relatedness G matrix
  if (method == 'relatedness'){
    G <- A.mat(as.matrix(x)-1)
    if (verbose >= 2){
      cat("  Calculating relatedness among individuals (G matrix)\n")
    }
  }  
  
    dd <- as.dist(dd) 
    
  # Revert to original order  
    ord <- rank(popNames(x))
    mat <- as.matrix(dd)[ord, ord]
    dd <- as.dist(mat)
    
}  
  
if (data.type == "SilicoDArT"){
  
#  method.list <- c("simple", "Jaccard", "Dice", "Sorenson", "Czekanowski", "Phi")

    if (method == 'simple'){
      if (verbose >= 2){cat("  Calculating the Simple Matching Index\n")}
    }
    if (method == 'jaccard'){
      if (verbose >= 2){cat("  Calculating the Jaccard Index\n")}
    }
    if (method == 'dice'){
      if (verbose >= 2){cat("  Calculating the Dice Index (= Sorenson or Czekanowski)\n")}
    }
    if (method == 'sorenson'){
      if (verbose >= 2){cat("  Calculating the Dice Index (= Sorenson or Czekanowski\n")}
    }
    if (method == 'czekanowski'){
      if (verbose >= 2){cat("  Calculating the Dice Index (= Sorenson or Czekanowski\n")}
    }
    if (method == 'phi'){
      if (verbose >= 2){cat("  Calculating the Pearson Phi Index (= Binary correlation\n")}
    }
    dd <- utils.dist.binary(x, method=method, verbose=verbose)
}
    
# PLOT
    if (plot){
      
      # Save the prior settings for mfrow, oma, mai and pty, and reassign
      op <- par(mfrow = c(2, 1), oma=c(1,1,1,1), mai=c(0.5,0.5,0.5,0.5),pty="m")
      
      # Set margins for first plot
      par(mai=c(1,0.5,0.5,0.5))
      
      # Plot Box-Whisker plot
      if (data.type=="SNP"){
        title <- paste0("SNP data (DArTSeq)\nInter-individual ",method," distance")
      } else {
        title <- paste0("Presence/Absence data (SilicoDArT)\nInter-individual ",method," distance")
      }  
      
      if (boxplot == "standard"){
        boxplot(dd, horizontal=TRUE, col='red', range=range, main = title)
        if(verbose >= 2){cat("  Standard boxplot, no adjustment for skewness\n")}
      } else {
        robustbase::adjbox(dd,
                           horizontal = TRUE,
                           col='red',
                           range=range,
                           main = title)
        if(verbose >= 2){cat("  Boxplot adjusted to account for skewness\n")}
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
    if(verbose >= 3){
      cat("\n  Reporting inter-individual distances\n")
      cat("  Distance measure:",method,"\n")
      cat("    No. of populations =", nPop(x), "\n")
      cat("    Average no. of individuals per population =", round(nInd(x)/nPop(x),1), "\n")
      cat("    No. of loci =", nLoc(x), "\n")
      cat("    Miniumum Distance: ",round(min(dd),2),"\n")
      cat("    Maximum Distance: ",round(max(dd),2),"\n")
      cat("    Average Distance: ",round(mean(dd),3),"\n\n")
    }  
    
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
    return(dd)
}

# # # Test script
# #D <- gl.dist.ind(testset.gl,method='relatedness',verbose=3)
# D <- gl.dist.ind(testset.gs,verbose=3)
# D <- gl.dist.ind(testset.gs,method='Simple',verbose=3)
# D <- gl.dist.ind(testset.gs,method='Jaccard',verbose=3)
# D <- gl.dist.ind(testset.gs,method='Sorenson',verbose=3)
# D <- gl.dist.ind(testset.gs,method='Phi',verbose=3)
# 
# D <- gl.dist.ind(testset.gs,method='Simple',verbose=3)
