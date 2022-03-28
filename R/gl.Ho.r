#' Estimates observed Heterozygosity
#'
#' @param gl A genlight object [required]
#' @return A simple vector whit Ho for each loci
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})

gl.Ho <- function(gl) {
  out <- colMeans(as.matrix(gl) == 1, na.rm = T)
  return(out)
}
