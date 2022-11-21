#' Identifies loci under selection per population using the outflank
#'  method of Whitlock and Lotterhos (2015)
#'
#' @param gi A genlight or genind object, with a defined population structure
#' [required].
#' @param plot A switch if a barplot is wanted [default TRUE].
#' @param LeftTrimFraction The proportion of loci that are trimmed from the
#' lower end of the range of Fst before the likelihood function is applied
#' [default 0.05].
#' @param RightTrimFraction The proportion of loci that are trimmed from the
#' upper end of the range of Fst before the likelihood function is applied
#' [default 0.05].
#' @param Hmin The minimum heterozygosity required before including calculations
#' from a locus [default 0.1].
#' @param qthreshold The desired false discovery rate threshold for calculating
#' q-values [default 0.05].
#' @param ... additional parameters (see documentation of outflank on github). 
#' @return Returns an index of outliers and the full outflank list
#' @details
#' This function is a wrapper around the outflank function provided by
#' Whitlock and Lotterhos. To be able to run this function the packages qvalue
#' (from bioconductor) and outflank (from github) needs to be installed. To do
#' so see example below.
#' @export
#' @importFrom stats optim pgamma quantile
#' @examples
#' \donttest{
#' gl.outflank(bandicoot.gl, plot = TRUE)
#' }
#' @references
#' Whitlock, M.C. and Lotterhos K.J. (2015) Reliable detection of loci
#' responsible for local adaptation: inference of a neutral model through
#' trimming the distribution of Fst. The American Naturalist 186: 24 - 36.
#'
#' Github repository: Whitlock & Lotterhos:
#'  \url{https://github.com/whitlock/OutFLANK} (Check the readme.pdf within the
#'   repository for an explanation. Be aware you now can run OufFLANK from a
#'   genlight object)
#' @seealso \code{\link{utils.outflank}}, \code{\link{utils.outflank.plotter}},
#'  \code{\link{utils.outflank.MakeDiploidFSTMat}}

gl.outflank <- function(gi,
                        plot = TRUE,
                        LeftTrimFraction = 0.05,
                        RightTrimFraction = 0.05,
                        Hmin = 0.1,
                        qthreshold = 0.05,
                        ...) {
    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "qvalue"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    # convert genlight to genind
    if (is(gi, "genlight")) {
        gi <- gl2gi(gi)
    }
    
    # missing value is 9!!! tempted to rewrite their model to be able to use genlight directly....
    snpmat <- as.matrix(gi)  #(matrix(NA, nrow=nind, ncol=nsnp)
    snpmat <- replace(snpmat, is.na(snpmat), 9)
    mdfm <-
        utils.outflank.MakeDiploidFSTMat(SNPmat = snpmat, list(colnames(snpmat)), list(as.character(gi@pop)))
    # run outflank
    outf <-
        utils.outflank(
            mdfm,
            LeftTrimFraction = LeftTrimFraction,
            RightTrimFraction = RightTrimFraction,
            Hmin = Hmin,
            NumberOfSamples = length(levels(gi@pop)),
            qthreshold = qthreshold
        )
    if (plot) {
        utils.outflank.plotter(outf)
    }
    index.outflank <-
        !(outf$results$OutlierFlag)  ## 6650 inliers and 188 outliers
    return(list(index = index.outflank, outflank = outf))
}
