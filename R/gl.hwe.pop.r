#' Function to run Hardy-Weinberg tests over each loci and population
#'
#' Does a Hardy-Weinberg-Test for each loci in each of the populations as defined by the pop slot in a genlight object. The function is completely rewritten as the previous function was based on package SNPassoc, which no longer is supported by CRAN. It now employs the \code{HardyWeinberg} package, which needs to be installed.  The function that is used is \code{\link[HardyWeinberg]{HWExactStats}}, but there are several other great function implemented in the package regarding HWE. Therefore you can return the data in the format the HWE package expects, via \code{HWformat=TRUE} and then use this to run other functions of the package.
#' @param x a genlight object with a population defined [pop(x) does not return NULL]
#' @param pvalue the p-value for the HWE test.
#' @param plot a switch if a barplot is returned.
#' @param HWformat switch if data should be returned in HWformat (counts of Genotypes to be used in package \code{HardyWeinberg})
#' @return This functions performs a HWE test for every population (rows) and loci (columns) and returns a true false matrix. True is reported if the p-value of an HWE-test for a particular loci and population was below the specified threshold (pvalue, default=0.05). The thinking behind this approach is that loci that are not in HWE in several populations have most likely to be treated (e.g. filtered if loci under selection are of interest). If plot=TRUE a barplot on the on the loci and the sum of deviation over all population is returned. Loci that deviate in the majority of populations can be identified via colSums on the resulting matrix. The function returns a list with up to three components: 'HWE' is the matrix over loci and populations, 'plot' is a matrixplot (ggplot) which shows the significant results for population and loci (can be amended further using ggplot syntax), and finally if 'HWEformat=TRUE' the 'HWformat' entails SNP data for each population in 'HardyWeinberg'-format to be used with other functions of the package (e.g \code{\link[HardyWeinberg]{HWPerm}} or \code{\link[HardyWeinberg]{HWExactPrevious}}).
#' @export
#' @examples 
#' out <- gl.hwe.pop(bandicoot.gl[,1:33], pvalue=0.05, plot=TRUE, HWformat=FALSE)


#########
### HWEperpop
########
gl.hwe.pop <- function(x, pvalue=0.05, plot=TRUE, HWformat=FALSE)
{
  # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "HardyWeinberg"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package",pkg," needed for this function to work. Please install it.") } 
  
  if (!is(x, "genlight")) stop("Provided object is not a genlight object! Please provide a valid input format.")
  
  if (is.null(pop(x))) stop("No population definintion is provided, therefore function terminates. If you want to have a test for a single population set population to a single population via: \npop(gl) <- rep(\"A\", nInd(gl))")

  pops <- adegenet::seppop(x)
  
  
  out <- lapply(pops, function(x) {
    pp <- as.matrix(x)
    xx <- data.frame(AA=colSums(pp==0, na.rm=T), AB=colSums(pp==1, na.rm=T), BB=colSums(pp==2, na.rm=T))
    HardyWeinberg::HWExactStats(xx)<pvalue
  } )
  #convert to matrix
  output <- matrix(unlist(out), ncol = nLoc(x), byrow = TRUE)
  rownames(output) <- names(pops)
  colnames(output) <- locNames(x)
  
  ggp <- NA
  loci <- counts <- NA
  if (plot) {
    longData<-reshape2::melt((as.matrix(output))*1)
    
    ggp <- ggplot(longData, aes(x = Var2, y = Var1,)) + 
      geom_raster(aes(fill=c("lightgrey", "red")[value+1])) + 
      scale_fill_identity() +
      labs(x="loci", y="population", title="HWE over population and loci") +
      theme_bw() + theme(axis.text.x=element_text(size=5, angle=90, vjust=0.3),
                         axis.text.y=element_text(size=9),
                         plot.title=element_text(size=11))
    
    print(ggp)
  }
  out=NULL
  if (HWformat)
  {
    out <- lapply(pops, function(x) {
      pp <- as.matrix(x)
      xx <- data.frame(AA=colSums(pp==0, na.rm=T), AB=colSums(pp==1, na.rm=T), BB=colSums(pp==2, na.rm=T))
      rownames(xx) <- locNames(x)
      xx })
    
    
  }
   
  return(list(HWE=output, plot=ggp,HWformat=out))
  
}
