#' Report minor allele frequency (MAF) for each locus in a genlight {adegenet} object
#'
#' Summary of minor allele frequencies across loci is reported as histograms.
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param maf.limit -- show histograms maf range <= maf.limit [default 0.5]
#' @param ind.limit -- show histograms only for populations of size greater than ind.limit [default 5]
#' @param loc.limit -- show histograms only for populations with more than loc.limit polymorphic loci [default 30]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @importFrom graphics layout
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' f <- gl.report.maf(testset.gl)


gl.report.maf <- function(x, maf.limit=0.5, ind.limit=5, loc.limit=30, v=2) {
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.drop.pop.r!\n"); stop("Execution terminated\n")
  }
  if (v > 0) {
    cat("Starting gl.report.maf: Minimum Allele Frequency\n")
  }
  if (maf.limit > 0.5 | maf.limit <= 0) {
    cat("Warning: maf.limit must be in the range (0,0.5], set to 0.5\n")
    maf.limit <- 0.5
  }
  if (ind.limit <= 0) {
    cat("Warning: ind.limit must be an integer > 0 and less than population size, set to 5\n")
    ind.limit <- 5
  }
  if (loc.limit <= 1) {
    cat("Warning: ind.limit must be an integer > 1 and less than the the total number of loci, set to 2\n")
    loc.limit <- 2
  }
  if (v < 0 | v > 5){
    cat("Warning: Verbosity must take on an integer value between 0 and 5, set to 3\n")
    v <- 3
  }
  layout(1,1)
  
# Recalculate the relevant loc.metrics
  
  if (v >= 3) {cat("  Removing monomorphic loci and recalculating FreqHoms and FreqHets\n")}
  
  x <- gl.filter.monomorphs(x, v = 0)
  x <- dartR:::utils.recalc.freqhets(x,v=0)
  x <- dartR:::utils.recalc.freqhomref(x,v=0)
  x <- dartR:::utils.recalc.freqhomsnp(x,v=0)
  
# Check for status -- any populations with loc > loc.limit; ind > ind.limit; and is nPop > 1
  
  count <- 0
  for (popn in popNames(x)) {
    genl <- x[pop(x)==popn]
    genl <- gl.filter.monomorphs(genl,v=0)
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
    if (v >= 3) {cat("  No populations met minimum limits on no of individuals or loci, reporting for overall\n")}
    title.str <- "Minor Allele Frequency\nOverall"
    flag <- 0
  }
  if (count == 1) {
    if (v >= 3) {cat("  Only one population met minimum limits on no of individuals or loci\n")}
    title.str <- paste("Minor Allele Frequency\n",popn.hold)
    flag <- 0
  }  
  if (nPop(x) == 1) {
    flag <- 0
    if (v >= 3) {cat("  Only one population specified\n")}
    title.str <- paste("Minor Allele Frequency\n",pop(x)[1])
  }
    
# Calculate and plot overall MAF
  
  if (v >= 3) {cat("  Calculating MAF across populations\n")}
  
  #cat("1Pops: ",nPop(x),"\n")
  #cat("1Inds: ",nInd(x),"\n")
  #cat("1Locs: ",nLoc(x),"\n")
  
  homref <- x@other$loc.metrics$FreqHomRef
  homalt <- x@other$loc.metrics$FreqHomSnp
  het <- x@other$loc.metrics$FreqHets
  
  maf <- array(0,length(homref))
  for (i in 1:length(homref)){
    maf[i] <- min((homref[i]*2 + het[i]), (homalt[i]*2 + het[i]))/2
  }

    if (flag == 1){
      layout(matrix(c(1,1,1,2,3,4,5,6,7), 3, 3, byrow = TRUE))
    }
    hist(maf, 
         breaks=seq(0,0.5,0.05), 
         col=rainbow(10), 
         main=title.str, 
         xlab="Minor Allele Frequency"
    )
  
  if (v >= 3) {cat("  Calculating MAF by population\n")}
  
  if (flag == 1){   
    plot.count <- 1
    for (popn in popNames(x)) {
      genl <- x[pop(x)==popn]
      genl <- gl.filter.monomorphs(genl, v = 0)
      if (nLoc(genl) >= loc.limit) {
      genl <- dartR:::utils.recalc.freqhets(genl,v=0)
      genl <- dartR:::utils.recalc.freqhomref(genl,v=0)
      genl <- dartR:::utils.recalc.freqhomsnp(genl,v=0)
      homref <- genl@other$loc.metrics$FreqHomRef
      homalt <- genl@other$loc.metrics$FreqHomSnp
      het <- genl@other$loc.metrics$FreqHets
      
      #cat(popn,"Pops: ",nPop(genl),"\n")
      #cat(popn,"Inds: ",nInd(genl),"\n")
      #cat(popn,"Locs******: ",nLoc(genl),"\n")
    
      maf <- array(0,length(homref))
      for (i in 1:length(homref)){
        maf[i] <- min((homref[i]*2 + het[i]), (homalt[i]*2 + het[i]))/2
      }
        
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
  
  
    if (v > 0) {
      if(plot.count > 6) {
        cat("Completed gl.report.maf\n  Refer to histograms which extend over multiple screens\n\n")
      } else {
        cat("Completed gl.report.maf\n  Refer to histograms\n\n")
      }
    }
  
  return()
}  
