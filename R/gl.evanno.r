#' @name gl.evanno
#'
#' @title Create an evanno plot from an sr structure run obejct
#'
#' @description 
#' This function takes a genlight object and runs a STRUCTURE analysis based on functions from \code{strataG}
#'
#' @param sr structure run object from \code{\link{gl.runstructure}} [required].
#' @param plot TRUE: all four plots are shown. FALSE: all four plots are returned by not shown.
#' @details The function is basically a convenient wrapper around the beautiful
#' strataG function \link[strataG]{evanno} (Archer et al. 2016). For a detailed
#' description please refer to this package (see references below).
#' @return an evanno plot is created and a list of all four plots is returned.

#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' #STRUCTURE needs to be installed to be able to run the example
#' #only the first 100 loci
#' #bc <- bandicoot.gl[,1:100]
#' #sr <- gl.runstructure(bc, k.range = 2:5, num.k.rep = 3, exec = "./structure.exe")
#' #gl.evanno(sr)
#' }
#' @import patchwork
###@importFrom strataG genind2gtypes structureRun
#'
#' @export
#' @seealso \code{\link{gl.runstructure}},  \link[strataG]{clumpp},
#' @references 
#' Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of population structure using multilocus genotype data. Genetics 155, 945-959.
#' 
#' Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R package for manipulating, summarizing and analysing population genetic data. Mol Ecol Resour. doi:10.1111/1755-0998.12559
#' 
#' Evanno, G., Regnaut, S., and J. Goudet. 2005. Detecting the number of clusters of individuals using the software STRUCTURE: a simulation study. Molecular Ecology 14:2611-2620.



gl.evanno <- function(sr, plot=TRUE)
{
  if (!requireNamespace("strataG", quietly = TRUE)) #not already installed?
  {
    ap <- available.packages()  #check CRAN
    
    oncran <- match("strataG", ap)
    if (is.na(oncran)) {
      warning("package strataG needs to be installed. It is currently not on CRAN, hence I try to install it manually via Github using devtools:\n  devtools::install_github('EricArcher/strataG'")
      devtools::install_github('EricArcher/strataG')
      
    }
  } else {
  
  
evno <- strataG::evanno(sr, plot=plot)
  }
}