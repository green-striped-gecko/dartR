#' Function to identify loci under selection per population using the outflank method of Whitlock and Lotterhos (2015)
#'  
#' @param gi a genlight of genind object, with a defined population structure 
#' @param plot a switch if a barplot is wanted.
#' @param LeftTrimFraction The proportion of loci that are trimmed from the lower end of the range of Fst before the likelihood funciton is applied.
#' @param RightTrimFraction The proportion of loci that are trimmed from the upper end of the range of Fst before the likelihood funciton is applied.
#' @param Hmin The minimum heterozygosity required before including calculations from a locus.
#' @param qthreshold The desired false discovery rate threshold for calculating q-values.
#' @param ... additional parameters (see documentation of outflank on github)
#' @return returns an index of outliers and the full outflank list
#' @details this function is a wrapper around the outflank function provided by Whitlock and Lotterhus. To be able to run this function the packages qvalue (from bioconductor) and outflank (from github) needs to be installed. To do so see example below.
#' @export
#' @importFrom stats optim pgamma quantile
#' @examples
#' \donttest{
#' gl.outflank(bandicoot.gl, plot = TRUE)
#' }
#' @references Whitlock, M.C. and Lotterhos K.J. (2015) Reliable detection of loci responsible for local adaptation: inference of a neutral model through trimming the distribution of Fst. The American Naturalist 186: 24 - 36.
#' 
#' Github repository: Whitlock & Lotterhos: \url{https://github.com/whitlock/OutFLANK} (Check the readme.pdf within the repository for an explanation. Be aware you now can run OufFLANK from a genlight object)
#' @seealso \code{\link{util.outflank}}, \code{\link{util.outflank.plotter}}, \code{\link{util.outflank.MakeDiploidFSTMat}} 




gl.outflank <- function(gi, plot=TRUE, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, qthreshold=0.05, ... )
{
# CHECK IF PACKAGES ARE INSTALLED
  pkg <- "qvalue"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package",pkg," needed for this function to work. Please install it.") }   
  

# convert genlight to genind 
if (is(gi,"genlight")) gi <- gl2gi(gi)
  
# missing value is 9!!! tempted to rewrite their model to be able to use genlight directly....
  snpmat <- as.matrix(gi)#(matrix(NA, nrow=nind, ncol=nsnp)
  snpmat <- replace(snpmat, is.na(snpmat), 9) 
  mdfm <- util.outflank.MakeDiploidFSTMat(SNPmat = snpmat, list(colnames(snpmat)), list(as.character(gi@pop)))
  #run outflank
  outf <- util.outflank(mdfm, LeftTrimFraction = LeftTrimFraction, RightTrimFraction=RightTrimFraction, Hmin = Hmin, NumberOfSamples = length(levels(gi@pop)), qthreshold = qthreshold)
  if (plot) util.outflank.plotter(outf)
  index.outflank <- !(outf$results$OutlierFlag) ## 6650 inliers and 188 outliers 
  return(list(index=index.outflank, outflank=outf))
}



