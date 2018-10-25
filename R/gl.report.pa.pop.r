#' Report private alleles (and fixed alleles) per pair of populations
#'
#' This function reports seperates the genlight object by populations and reports fixed alleles, the pairwise private alleles and the mean absolute allele frequency differences between pair of population. 
#'
#' @param gl -- name of the genlight object containing the SNP data [see Details]
#' @param gl2 -- if two seperate genlight objects are to be compared this can be provided here [see Details]
#' @return A data.frame will be returned. Each row shows for a pair of populations the number of individuals in a population, the number of loci with fixed differences (same for both populations) in pop1 (compared to pop2) and vice versa. Same for private alleles and finally the absolute mean allele frequendy difference between loci (mdf).
#' @details 
#' if no gl2 is provided, the function uses the pop(gl) hierachy to determine pairs of population, otherwise it runs a single comparison between gl and gl2. Hint: in case you want to run comparison between individuals you can simply redefine your pop(gl) via indNames(gl) [Assuming individual names are unique]
#'  
#' Definition of fixed and private allels
#' 
#' The table shows a cross table of possible cases of allele frequencies between two populations (0=homozygote for Allele 1,x= both Alleles are present, 1=homozygote for Allele 2)
#' 
#' p: cases where there is a private allele in pop1 compared to pop2 (but not vice versa)
#' 
#' f: cases where there is a fixed allele in pop1 (and pop2, as those cases are symmetric)
#'
#'\tabular{ccccc}{ 
#' \tab\tab \tab \emph{pop1}\tab\cr
#' \tab\tab \strong{0} \tab   \strong{x}  \tab  \strong{1}\cr
#' \tab     \strong{0}\tab -  \tab  p \tab  p,f\cr
#'  \emph{pop2} \tab \strong{x}\tab -  \tab- \tab -\cr
#' \tab \strong{1} \tab p,f\tab p \tab   -\cr
#' }
#' @export
#' @author Bernd Gruber (glbugs@aerg.canberra.edu.au)
#' @examples
#' gl.report.pa.pop(testset.gl[1:20,])
#' 
#' 
#' 
#private allellels per populations

gl.report.pa.pop <- function(gl, gl2=NULL)
  
{
  
  if (!is.null(gl2)) pops <- list(pop1=gl, pop2=gl2) else 
  {
   if (length(unique(pop(gl)))>1) pops <- seppop(gl) else stop("Only one population provided. Check the @pop slot in your genlight object.\n ")
  }
  
  
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
  
  return(pall)
}
