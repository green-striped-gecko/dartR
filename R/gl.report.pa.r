#' Report private alleles (and fixed allelic differences) per pair of populations
#'
#' This function reports separates the genlight object by populations and reports fixed allelic differences, private alleles 
#' and the mean absolute allele frequency differences between pairs of populations. 
#'
#' @param x -- name of the genlight object containing the SNP data broken down by population [default gl]
#' @param nmin -- minimum sample size for a population to be included in the analysis [default 10]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 3]
#' @return A dataframe giving population names, sample sizes, number of fixed allelic differences, number of private alleles held by population 1,
#' number of private alleles held by population 2, total private alleles, and mean allelic difference between the two populations.
#' @details 
#' A fixed difference at a locus occurs when the two focal populations share no alleles. This is a symmetric measure. A private allele is an allele
#' that is found in one population but not in the other. The distance defined by a count of private alleles between two populations is not symmetric.
#' 
#' The table below shows a cross tablulation of possible cases of allele frequencies between two populations 
#' (0=homozygote for Allele 1,x= both Alleles are present, 1=homozygote for Allele 2)
#' 
#' p: cases where there is a private allele in pop1 compared to pop2 (but not vice versa)
#' 
#' f: cases where there is a fixed allelic difference between pop1 and pop2.
#'
#'\tabular{ccccc}{ 
#' \tab\tab \tab \emph{pop1}\tab\cr
#' \tab\tab \strong{0} \tab   \strong{x}  \tab  \strong{1}\cr
#' \tab     \strong{0}\tab -  \tab  p \tab  p,f\cr
#'  \emph{pop2} \tab \strong{x}\tab -  \tab- \tab -\cr
#' \tab \strong{1} \tab p,f\tab p \tab   -\cr
#' }
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.report.pa(testset.gl,nmin=11)
#' 

gl.report.pa <- function(x=gl, nmin=10, verbose=3){
  
  
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
    cat("  Fatal Error: at least one genlight object is required!\n"); stop("Execution terminated\n")
  }
  
  # SCRIPT SPECIFIC ERROR CHECKING
  
  # DO THE JOB  
  
    gl <- x
  
   if (length(unique(pop(gl)))>1) {
     keep1 <- levels(pop(gl))[table(pop(gl)) >= nmin]
     if (verbose >=2) {cat("  Retaining",length(keep1),"populations with sample size greater than or equal to",nmin,":",keep1,"\n")}
     #toss1 <- levels(pop(gl))[table(pop(gl)) < nmin]
     #if (verbose >=2) {cat("  Discarding",length(toss1),"populations with sample size less than",nmin,":",toss1,"\n\n")}
     
     gl.keep <- gl.keep.pop(gl,pop.list=keep1,verbose=0)
     pops <- seppop(gl.keep)
   } else {
       stop("  Fatal Error: Only one population provided. Check pop(gl)\n ")
   }
  pc <- t(combn(length(pops),2))
  pall <- data.frame(p1=pc[,1], p2=pc[,2], pop1=names(pops)[pc[,1]], pop2=names(pops)[pc[,2]], N1=NA, N2=NA,fixed=NA, priv1=NA, priv2=NA, totalpriv=NA, mdf=NA)
  
  for (i in 1:nrow(pc)) {
    i1 =pall[i,1]
    i2 =pall[i,2]
    
    p1 <- as.matrix(pops[[i1]])
    p2 <- as.matrix(pops[[i2]])
    p1alf <- colMeans(p1, na.rm = T)/2
    p2alf <- colMeans(p2, na.rm = T)/2
    
    pall[i,5:6] <- c(nrow(p1), nrow(p2))
    pall[i,7] = sum(abs(p1alf-p2alf)==1, na.rm=T)
    
    pall[i,8] =  sum(p2alf==0 & p1alf!=0, na.rm=T) + sum(p2alf==1 & p1alf!=1, na.rm = T) 
    pall[i,9] =  sum(p1alf==0 & p2alf!=0, na.rm=T) + sum(p1alf==1 & p2alf!=1, na.rm = T)  
    pall[i,10] = pall[i,8]+pall[i,9]
    pall[i,11] = round(mean(abs(p1alf-p2alf), na.rm=T),3)
  }
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat("Completed:",funname,"\n\n")
  }

  # Return the result  
  return(pall)
}
