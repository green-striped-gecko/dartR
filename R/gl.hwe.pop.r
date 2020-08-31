#' Filter function to facilitate analysing of dart data
#'  
#' This function needs package SNPassoc. At the time of writing (August 2020) the package was no longer available from CRAN. To install the package check their github repository. \url{https://github.com/isglobal-brge/SNPassoc} and/or use \code{install_github("isglobal-brge/SNPassoc")} to install the function. 
#' @param gi a genlight or genind object(genlight is internally converted to genind)
#' @param pvalue the p-value for the HWE test.
#' @param plot a switch if a barplot is wanted.
#' @return This functions performs a HWE test for every population (rows) and loci (columns) and returns a true false matrix. True is reported if the p-value of an HWE-test for a particular loci and population was below the specified threshold (pvalue, default=0.05). The thinking behind this approach is that loci that are not in HWE in several populations have most likely to be treated (e.g. filtered if loci under selection are of interest). If plot=TRUE a barplot on the on the loci and the sum of deviation over all population is returned. Loci that deviate in the majority of populations can be identified via colSums on the resulting matrix.
#' @export
#' @examples
#' if (!requireNamespace("parallel", quietly = TRUE)) {
#'  stop("Package parallel needed for this function to work. Please install it.")} else {
#' library(parallel)
#' out <- gl.hwe.pop(possums.gl[,1:20], pvalue = 0.05, plot = TRUE)
#' out$indHWE
#' out$plot
#' }

#########
### HWEperpop
########
gl.hwe.pop <- function(gi, pvalue=0.05, plot=TRUE)
{
# CHECK IF PACKAGES ARE INSTALLED
  if (!(requireNamespace("parallel", quietly = TRUE))) {
    stop("To use this function you need to install package: parallel.")
  }
  pkg <- "pegas"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package",pkg," needed for this function to work. Please   install it.")  }
  
  if (!(requireNamespace("SNPassoc", quietly = TRUE))) {
    stop("To use this function you need to install package: SNPassoc. Please refer to the help of the function for instructions (?gl.hwe.pop).")
  } else {
  
  
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
  colnames(indices.HWE.pop)<- locNames(gi)
  ggp <- NA
  loci <- counts <- NA
  if (plot) {
    df <- data.frame(loci=colnames(indices.HWE.pop), counts=colSums(indices.HWE.pop, na.rm = T))
    
    ggp <-ggplot(df, aes( x=loci, y=counts))+geom_bar(stat="identity")+theme_classic()
    ggp
  }
  
  return(list(indHWE=indices.HWE.pop, plot=ggp))
  }
}
