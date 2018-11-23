#' A very simple function to report observed Heterozygosity
#' 
#' @param gl -- a genlight object
#' @return a simple vector wiht Ho for each loci
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})

gl.Ho <- function(gl)
{
  out <- colMeans(as.matrix(gl)==1,na.rm=T)
  return(out)
}
