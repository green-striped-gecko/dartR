#' Calculates basic statistics for each loci (Hs, Ho, Fis etc.)
#' 
#' Based on function \code{\link[hierfstat]{basic.stats}}. Check ?basic.stats for help. 
#' @param gl -- a genlight object
#' @param digits -- number of digits that should be returned
#' @importFrom hierfstat genind2hierfstat basic.stats
#' @return several tables and lists with all basic stats. Check \code{\link[hierfstat]{basic.stats}} for details.
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})

gl.basic.stats <- function(gl, digits=4)
{
  out <- basic.stats(genind2hierfstat(gl2gi(gl, verbose = 0)), digits = digits)
  return(out)
}
