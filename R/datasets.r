#' A genlight object created via the read.dart functions
#'
#' @name testset.gl
#' @format genlight object
#' @docType data
#' @author Arthur Georges \email{(glbugs@@aerg.canberra.edu.au}
#' @keywords datasets
"testset.gl"

#' Testfile in DArT format (as provided by DArT)
#' 
#' this test data set is provided to show a typical DArT file format. Can be used to create a genlight object using the read.dart function.
#' @name testset_SNPs_2Row
#' @format csv
#' @docType data
#' @author Arthur Georges \email{(glbugs@@aerg.canberra.edu.au}
#' @keywords datasets
NULL

#' Recode file to be used with the function.
#' 
#' This test data set is provided to show a typical recode file format.
#' @name testset_pop_recode
#' @format csv
#' @docType data
#' @author Arthur Georges \email{(glbugs@@aerg.canberra.edu.au}
#' @keywords datasets
NULL

#' Metadata file. Can be integrated via the dart2genlight function.
#' 
#' @name testset_metadata
#' @format csv
#' @docType data
#' @author Arthur Georges \email{(glbugs@@aerg.canberra.edu.au}
#' @keywords datasets
NULL

#' Example data set as text file to be imported into a genlight object
#' 
#' Check ?read.genetable in pacakge PopGenReport for details on the format. 
#' @name platy
#' @format csv
#' @docType data
#' @author Bernd Gruber \email{(glbugs@@aerg.canberra.edu.au}
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

