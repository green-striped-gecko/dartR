#' Calculates allele frequency of the first and second allele for each loci
#' A very simple function to report allele frequencies
#' @param x Name of the genlight object containing the SNP data [required].
#' @return A simple data.frame with alf1, alf2.
#' @export
#' @rawNamespace import(adegenet, except = plot)
#' @author Bernd Gruber (bugs? Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' #for the first 10 loci only
#' gl.alf(possums.gl[,1:10])
#' barplot(t(as.matrix(gl.alf(possums.gl[,1:10]))))

gl.alf <- function(x) {
  alf <- colMeans(as.matrix(x), na.rm = TRUE) / 2
  out <- data.frame(alf1 = 1 - alf, alf2 = alf)
  return(out)
}
