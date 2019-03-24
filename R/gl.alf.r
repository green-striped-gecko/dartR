#' Calculates allele frequency of the first and second allele for each loci #' 
#' A very simple function to report allele frequencies
#' @param gl -- a genlight object
#' @return a simple data.frame with alf1, alf2
#' @export
#' @rawNamespace import(adegenet, except = plot)
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})

gl.alf <- function(gl)
{
  alf <- colMeans(as.matrix(gl), na.rm = T)/2
  out <- data.frame(alf1=1-alf, alf2=alf)
  return(out)
}
