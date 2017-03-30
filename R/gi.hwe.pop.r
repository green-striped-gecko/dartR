#' Filter function to facilitate analysing of dart data
#'  
#' @param gi a genind object (often created using the gl2gi function)
#' @param pvalue the p-value for the HWE test.
#' @param plot a switch if a barplot is wanted.
#' @return tests HWE for every loci in every population (using the setupSNP function in package SNPassoc
#' @export
#' @examples
#' \dontrun{
#' HWEperpop(gi=geni, pvalue = 0.05, plot = TRUE)
#' }

#########
### HWEperpop
########
gi.hwe.pop <- function(gi, pvalue=0.05, plot=TRUE)
{
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
  indices.HWE.pop
}

