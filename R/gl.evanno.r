#' @name gl.evanno
#' @title Create an Evanno plot from an sr structure run object
#' @description
#' This function takes a genlight object and runs a STRUCTURE analysis based on
#' functions from \code{strataG}
#' @param sr structure run object from \code{\link{gl.run.structure}} [required].
#' @param plot.out TRUE: all four plots are shown. FALSE: all four plots are
#' returned as a ggplot but not shown [default TRUE].
#' @details The function is basically a convenient wrapper around the beautiful
#' strataG function \link[strataG]{evanno} (Archer et al. 2016). For a detailed
#' description please refer to this package (see references below).
#' @return An Evanno plot is created and a list of all four plots is returned.
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' #CLUMPP and STRUCTURE need to be installed to be able to run the example
#' #bc <- bandicoot.gl[,1:100]
#' #sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3, exec = './structure.exe')
#' #ev <- gl.evanno(sr)
#' #ev
#' #qmat <- gl.plot.structure(sr, k=3, CLUMPP='d:/structure/')
#' #head(qmat)
#' #gl.map.structure(qmat, bc, scalex=1, scaley=0.5)
#' }
#' @import patchwork
### @importFrom strataG genind2gtypes structureRun
#' @export
#' @seealso \code{\link{gl.run.structure}},  \link[strataG]{clumpp},
#' @references
#' \itemize{
#' \item Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of
#' population structure using multilocus genotype data. Genetics 155, 945-959.
#' \item Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R
#' package for manipulating, summarizing and analysing population genetic data.
#' Mol Ecol Resour. doi:10.1111/1755-0998.12559
#' \item Evanno, G., Regnaut, S., and J. Goudet. 2005. Detecting the number of
#'  clusters of individuals using the software STRUCTURE: a simulation study.
#'   Molecular Ecology 14:2611-2620.
#' }

gl.evanno <- function(sr,
                      plot.out = TRUE) {
    # check for strataG (currently not on CRAN)
    if (!requireNamespace("strataG", quietly = TRUE)) {
        ap <- available.packages()  #check CRAN
        
        oncran <- match("strataG", ap)
        if (is.na(oncran)) {
            warning(
                warn(
                    "package strataG needs to be installed. It is currently not on CRAN, hence you need to install it manually via Github using devtools:\n  devtools::install_github('EricArcher/strataG')\n. Do you want to install strataG from github now? (Y/N)+Enter\n"
                )
            )
            line <- readline()
            if (substr(line, 1, 1) == "Y") {
                cat(report("\nInstalling strataG from Github...\n"))
                flush.console()
                devtools::install_github("EricArcher/strataG")
            } else{
                stop(
                    error(
                        "Package strataG is necessary to run structure. Check the help page for references ?gl.run.structure. Aborting now."
                    )
                )
            }
            
        }
    } else {
        evno <- strataG::evanno(sr, plot = plot.out)
        # for (i in 1:4) evno$plots[[i]] <- evno$plots[[i]]+theme_dartR()
        
        return(evno)
    }
}
