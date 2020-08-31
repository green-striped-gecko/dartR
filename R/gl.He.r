#' A very simple function to report expected Heterozygosity
#' 
#' @param gl -- a genlight object
#' @return a simple vector wiht Ho for each loci
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})

gl.He <- function(gl)
{
  alf <- colMeans(as.matrix(gl), na.rm = T)/2
  out <- alf*(1-alf)*2
  return(out)
}
