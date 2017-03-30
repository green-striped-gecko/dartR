#' Functio to identify loci under selection per population using the outflank method of Whitlock and Lotterhos (2015)
#'  
#' @param gi a genind object, with a defined population structure 
#' @param plot a switch if a barplot is wanted.
#' @return returns an index of outliers and the full outflank list
#' @details this function is a wrapper around the outflank function provided by Whitlock and Lotterhus. To be able to run this function the packages qvalue (from bioconductor) and outflank (from github) needs to be installed. To do so see example below.
#' @export
#' @examples
#' \dontrun{
#' #install qvalue from bioconductor and outflank from github 
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("qvalue") 
#' library(devtools)
#' install_github("Whitlock/OutFLANK")
#' # now we are able to run outflank using a genlight object
#' gi.outflank(gi, plot = TRUE)
#' }
#' @references Whitlock, M.C. and Lotterhos K.J. (2015) Reliable detection of loci responsible for local adaptation: inference of a neutral model through trimming the distribution of Fst. The American Naturalist 186: 24 - 36.

gi.outflank <- function(gi, plot=TRUE)
{
  #library(OutFLANK)
  #library(qvalue)
  # outflank requires (of course) a different format 
  # missing value is 9!!! tempted to rewrite their model to be able to use genlight directly....
  snpmat <- as.matrix(gi)#(matrix(NA, nrow=nind, ncol=nsnp)
  snpmat <- replace(snpmat, is.na(snpmat), 9) 
  mdfm <- OutFLANK::MakeDiploidFSTMat(SNPmat = snpmat, list(colnames(snpmat)), list(as.character(gi@pop)))
  #run outflank
  outf <- OutFLANK::OutFLANK(mdfm, LeftTrimFraction = 0.05, RightTrimFraction = 0.05, Hmin = 0.1, NumberOfSamples = length(levels(gi@pop)), qthreshold = 0.05)
  if (plot) OutFLANK::OutFLANKResultsPlotter(outf)
  #the filter if we believe in Whitlock and why wouldn't we?
  index.outflank <- !(outf$results$OutlierFlag) ## 6650 inliers and 188 outliers 
  return(list(index=index.outflank, outflank=outf))
  #sum(!index.outflank) #number of outliers 
}



