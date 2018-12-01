#' A simulated genlight object created to run a landscape genetic example
#'
#'This a test data set to run a landscape genetics example. It contains 10 populations of 30 individuals each and each individual has 300 loci. There are no covariates for individuals or loci.
#' @name possums.gl
#' @format genlight object
#' @docType data
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr}
#' @keywords datasets
"possums.gl"

#' A genlight object created via the read.dart functions
#'
#'This a test data set to test the validity of functions within dartR and is based on a DArT SNp data set of foxes across Australia. It contains 100 individuals and 1000 SNPs.
#' @name foxes.gl
#' @format genlight object
#' @docType data
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr}
#' @keywords datasets
"foxes.gl"

#' A genlight object created via the read.dart functions
#'
#' @name testset.gl
#' @format genlight object
#' @docType data
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr}
#' @keywords datasets
"testset.gl"

#' Testfile in DArT format (as provided by DArT)
#' 
#' this test data set is provided to show a typical DArT file format. Can be used to create a genlight object using the read.dart function.
#' @name testset_SNPs_2Row
#' @format csv
#' @docType data
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr}
#' @keywords datasets
NULL

#' Recode file to be used with the function.
#' 
#' This test data set is provided to show a typical recode file format.
#' @name testset_pop_recode
#' @format csv
#' @docType data
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr}
#' @keywords datasets
NULL

#' Metadata file. Can be integrated via the dart2genlight function.
#' 
#' @name testset_metadata
#' @format csv
#' @docType data
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr}
#' @keywords datasets
NULL

#' Example data set as text file to be imported into a genlight object
#' 
#' Check ?read.genetable in pacakge PopGenReport for details on the format. 
#' @name platy
#' @format csv
#' @docType data
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr}
#' @keywords datasets
#' @examples 
#' \donttest{
#' library(PopGenReport)
#' read.csv( paste(.libPaths()[1],"/dartR/extdata/platy.csv",sep="" ))
#' platy <- read.genetable( paste(.libPaths()[1],"/dartR/extdata/platy.csv",
#' sep="" ), ind=1, pop=2, lat=3, long=4, other.min=5, other.max=6, oneColPerAll=FALSE,
#' sep="/")
#' platy.gl <- (gi2gl(platy))
#' df.loc <- data.frame(RepAvg = runif(nLoc(platy.gl)), CallRate = 1)
#' platy.gl@other$loc.metrics <- df.loc
#' gl.report.repavg(platy.gl)
#' }
NULL

