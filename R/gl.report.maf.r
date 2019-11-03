#' Report minor allele frequency (MAF) for each locus in a SNP dataset
#'
#' This script provides summary histograms of MAF for each population in the dataset as a basis for decisions on filtering.
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param maf.limit -- show histograms maf range <= maf.limit [default 0.5]
#' @param ind.limit -- show histograms only for populations of size greater than ind.limit [default 5]
#' @param loc.limit -- show histograms only for populations with more than loc.limit polymorphic loci [default 30]
#' @param silent -- if FALSE, function returns an object, otherwise NULL [default TRUE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @importFrom graphics layout hist
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.report.maf(testset.gl,silent=FALSE)
#'

gl.report.maf <- function(x, maf.limit=0.5, ind.limit=5, loc.limit=30, silent=TRUE, verbose = 2) {
  
# TIDY UP FILE SPECS
  
  build ='Jacob'
  funname <- match.call()[[1]]
  # Note does not draw upon or modify the loc.metrics.flags
  
# FLAG SCRIPT START
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  cat("Starting",funname,"[ Build =",build,"]\n")
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    cat("  Detected Presence/Absence (SilicoDArT) data\n")
    stop("Cannot calculate minor allele frequences for Tag presence/absence data. Please provide a SNP dataset.\n")
  } else if (all(x@ploidy == 2)){
    cat("  Processing a SNP dataset\n")
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (Tag P/A data) or 2 (SNP data)!\n")
  }

# FUNCTION SPECIFIC ERROR CHECKING
  
  if (maf.limit > 0.5 | maf.limit <= 0) {
    cat("Warning: maf.limit must be in the range (0,0.5], set to 0.5\n")
    maf.limit <- 0.5
  }
    
  if (ind.limit <= 0) {
    cat("Warning: ind.limit must be an integer > 0 and less than population size, set to 5\n")
    ind.limit <- 5
  }
  
  if (loc.limit <= 1) {
    cat("Warning: loc.limit must be an integer > 1 and less than the the total number of loci, set to 2\n")
    loc.limit <- 2
  }

# DO THE JOB

  layout(1,1)
  
# Recalculate the relevant loc.metrics
  
  if(x@other$loc.metrics.flags$maf==FALSE){
    cat("  Recalculating MAF\n")
    x <- utils.recalc.maf(x,verbose=1)
  }  

# Check for status -- any populations with loc > loc.limit; ind > ind.limit; and is nPop > 1
  
  count <- 0
  for (popn in popNames(x)) {
    genl <- x[pop(x)==popn]
    genl <- gl.filter.monomorphs(genl,verbose=0)
    if (nInd(genl) >= ind.limit & nLoc(genl) >= loc.limit) {
      count <- count + 1
      popn.hold <- popn
    }
  }
  if (count > 1 ) {
    flag <- 1
    title.str <- "Minor Allele Frequency\nOverall"
  }
  if (count == 0) {
    if (verbose >= 1) {cat("  No populations met minimum limits on no of individuals or loci, reporting for overall\n")}
    title.str <- "Minor Allele Frequency\nOverall"
    flag <- 0
  }
  if (count == 1) {
    if (verbose >= 3) {cat("  Only one population met minimum limits on no of individuals or loci\n")}
    title.str <- paste("Minor Allele Frequency\n",popn.hold)
    flag <- 0
  }  
  if (nPop(x) == 1) {
    flag <- 0
    if (verbose >= 1) {cat("  Only one population specified\n")}
    title.str <- paste("Minor Allele Frequency\n",pop(x)[1])
  }
    
# Calculate and plot overall MAF
  
  cat("  Calculating MAF across populations\n")
  maf <- x@other$loc.metrics$maf

    if (flag == 1){
      layout(matrix(c(1,1,1,2,3,4,5,6,7), 3, 3, byrow = TRUE))
    }
    hist(maf, 
         #breaks=seq(0,0.5,0.05), 
         col='red', 
         main=title.str, 
         xlab="Minor Allele Frequency",
         breaks=100)
  
 cat("  Calculating MAF by population\n")
  
  plot.count <- 1
  if (flag == 1){   
    for (popn in popNames(x)) {
      genl <- x[pop(x)==popn]
      genl <- gl.filter.monomorphs(genl, verbose= 0)
      if (nLoc(genl) >= loc.limit) {
      genl <- utils.recalc.maf(genl,verbose=0)
      maf <- genl@other$loc.metrics$maf
      
      if (plot.count <= 6) {
        maf <- maf[maf<maf.limit]
        hist(maf, 
            breaks=seq(0,maf.limit,len=10), 
            col=rainbow(10), 
            main=paste(popn,"\nn =",nInd(genl)), 
            xlab="Minor Allele Frequency",
            xlim=c(0,maf.limit)
        )
        plot.count <- plot.count + 1
      }  
      }
        
      if (plot.count == 7) {
        layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow = TRUE))
        cat ("  Overflow of plots across multiple pages\n") 
      }
        
      if (plot.count >= 7){
        maf <- maf[maf<maf.limit]
        hist(maf, 
           breaks=seq(0,maf.limit,len=10), 
           col=rainbow(10), 
           main=paste(popn,"\nn =",nInd(genl)), 
           xlab="Minor Allele Frequency",
           xlim=c(0,maf.limit)
        )
        plot.count <- plot.count + 1  
      }
    }   
  }
  
  if(plot.count > 6) {
        cat("Completed: gl.report.maf once plots are displayed\n  Refer to histograms which extend over multiple screens\n\n")
  } else {
        cat("Completed: gl.report.maf once plots are displayed\n  Refer to histograms\n\n")
  }
  
# FLAG SCRIPT END
  
  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  if(silent==TRUE){
    return(NULL)
  } else{
    return(NULL)
  } 

}  
