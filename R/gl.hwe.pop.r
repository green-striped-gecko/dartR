#' Filter function to facilitate analysing of dart data
#'  
#' @param gi a genlight or genind object(genlight is internally converted to genind)
#' @param pvalue the p-value for the HWE test.
#' @param plot a switch if a barplot is wanted.
#' @return This functions performs a HWE test for every population (rows) and loci (columns) and returns a true false matrix. True is reported if the p-value of an HWE-test for a particular loci and population was below the specified threshold (pvalue, default=0.05). The thinking behind this approach is that loci that are not in HWE in several populations have most likely to be treated (e.g. filtered if loci under selection are of interest). If plot=TRUE a barplot on the on the loci and the sum of deviation over all population is returned. Loci that deviate in the majority of populations can be identified via colSums on the resulting matrix.
#' @export
#' @examples
#' \dontrun{
#' gl.hwe.pop(gi=gi, pvalue = 0.05, plot = TRUE)
#' }

#########
### HWEperpop
########
gl.hwe.pop <- function(gi, pvalue=0.05, plot=TRUE)
{
# convert genlight to genind 
  if (class(gi)=="genlight") gi <- gl2gi(gi)
  
  
#  library(SNPassoc)   #package for LD and Hardy Weinberg...
#  library(pegas)
  sep.troch <- seppop(gi)  #creates a list of all populations
  #define a function
  hwepop <- function(pop)
  {
    xxx <- pegas::as.loci(pop)[,-1]
    sup <-  SNPassoc::setupSNP(data.frame(xxx), 1:ncol(xxx)) 
    sum.sup <- summary(sup)
    index.HWE <- sum.sup$HWE < pvalue 
    index.HWE <- ifelse(is.na(index.HWE), FALSE, index.HWE)      #what to do with NA?????
    index.HWE
  }
  indices.HWE.pop <- do.call(rbind,lapply(sep.troch, hwepop)) #let's do that for every population
  if (plot) barplot( colSums(indices.HWE.pop,na.rm=T))
  colnames(indices.HWE.pop)<- locNames(gi)
  indices.HWE.pop
}

