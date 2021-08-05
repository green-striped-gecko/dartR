#' @name gl.report.pa
#'
#' @title Report private alleles (and fixed alleles) per pair of populations
#'
#' @description 
#' This function reports private alleles in one population compared with a second population, for all populations
#' taken pairwise. It also reports a count of fixed allelic differences and the mean absolute allele frequency
#' differences between pairs of populations.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param x2 If two separate genlight objects are to be compared this can be provided here [default NULL].
#' @param verbose Verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity].
#' 
#' @details 
#' Note that the number of paired alleles between two populations is not a symmetric dissimilarity measure.
#' 
#' If no x2 is provided, the function uses the pop(gl) hierarchy to determine pairs of population, otherwise it runs a single comparison between gl1 and gl2. 
#' Hint: in case you want to run comparison between individuals you can simply redefine your pop(gl) via indNames(gl) [Assuming individual names are unique]
#'  
#'\strong{ Definition of fixed and private alleles }
#' 
#' The below is table showing the possible cases of allele frequencies between two populations 
#' (0 = homozygote for Allele 1, x = both Alleles are present, 1 = homozygote for Allele 2).
#' 
#'\itemize{
#'\item p: cases where there is a private allele in pop1 compared to pop2 (but not vice versa)
#'\item f: cases where there is a fixed allele in pop1 (and pop2, as those cases are symmetric)
#'} 
#' 
#'\tabular{ccccc}{ 
#'\tab \tab \tab \emph{pop1} \tab \cr
#'\tab \tab \strong{0} \tab   \strong{x}  \tab  \strong{1}\cr
#'\tab     \strong{0}\tab -  \tab  p \tab  p,f\cr
#'\emph{pop2} \tab \strong{x}\tab -  \tab- \tab -\cr
#'\tab \strong{1} \tab p,f\tab p \tab   -\cr
#' }
#'   
#' @return A data.frame. Each row shows for a pair of populations the number of individuals in a population, the number of loci with fixed differences (same for both populations) in pop1 (compared to pop2) and vice versa. Same for private alleles and finally the absolute mean allele frequency difference between loci (mdf).
#'
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' out <- gl.report.pa(testset.gl[1:20,])
#'
#' @seealso \code{\link{gl.list.reports}},
#'  \code{\link{gl.print.reports}}
#'  
#' @family reporting functions
#'
#' @export
#'

gl.report.pa <- function(x, 
                         x2 = NULL, 
                         verbose = NULL){
  # TRAP COMMAND
  
  funname <- match.call()[[1]]
  
  # SET VERBOSITY
  
  verbose <- gl.check.verbosity(verbose)
  
  # CHECKS DATATYPE 
  
  datatype <- utils.check.datatype(x, verbose=verbose)
  
  if(!is.null(x2)){
    datatype <- utils.check.datatype(x2)
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat(report("Starting",funname,"[ Build =",build,"]\n\n"))
    } else {
      cat(report("Starting",funname,"\n\n"))
    }
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop(error("Fatal Error: genlight object required!\n"))
  }
  if (all(x@ploidy == 1)){
    stop(error("Cannot calculate minor allele frequences for Tag presence/absence data. Please provide a SNP dataset.\n"))
  } else if (all(x@ploidy == 2)){
    if(verbose>=2){
      cat(report("  Processing a SNP dataset\n"))
      }
  } else {
    stop(error("Fatal Error: Ploidy must be universally 1 (Tag P/A data) or 2 (SNP data)!\n"))
  }
  
  if(!is.null(x2)){
    if(class(x2)!="genlight") {
      stop(error("Fatal Error: genlight object required for gl2!\n"))
    }
    if (all(x2@ploidy == 1)){
      stop(error("Fatal Error: Private alleles can only be calculated for SNP data. Please provide a SNP dataset for gl2\n"))
    } else if (all(x2@ploidy == 2)){
      if (verbose >= 2){
        cat(report("  Processing a SNP dataset",x2,"\n"))
        }
    } else {
      stop(error("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)"))
    }
  }

# FUNCTION SPECIFIC ERROR CHECKING
  
if (!is.null(x2)) {
  pops <- list(pop1=x, pop2=x2) 
  } else {
   if (length(unique(pop(x)))>1){
     pops <- seppop(x)
  } else {
    stop(error("Only one population provided. Check the @pop slot in your genlight object.\n "))
  }  
}
 
# DO THE JOB
  
  pc <- t(combn(length(pops),2))
  pall <- data.frame(p1=pc[,1], p2=pc[,2], pop1=names(pops)[pc[,1]], pop2=names(pops)[pc[,2]], N1=NA, N2=NA,fixed=NA, priv1=NA, priv2=NA, totalpriv=NA, mdf=NA)
  
  for (i in 1:nrow(pc))
  {
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
  
  # PRINTING OUTPUTS
  if(verbose >= 3){
    print(pall)
  }
  
  if(verbose >= 2){
    cat(report("  Table of private alleles and fixed differences returned\n"))
        }
  
# FLAG SCRIPT END

    if (verbose >= 1) {
    cat(report("\n\nCompleted:", funname, "\n\n"))
  }
  
  # RETURN
  
  invisible(pall)

}
