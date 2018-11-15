#' Performs and AMOVA using genlight data. 
#'
#' This script performs an AMOVA based on the genetic distance matrix from stamppNeisD() [package StAMPP] using the amova() function from the package PEGAS for exploring within and between population variation. For detailed information use their help pages: ?pegas::amova, ?StAMPP::stamppAmova
#' 
#' @param x -- name of the genlight containing the SNP genotypes, with population information [required]
#' @param nperm -- number of permuations to perform for hypothesis testing [default 100]. Please note should be set to 1000 for analysis.]
#' @return An object of class "amova" which is a list with a table of sums of square deviations (SSD), mean square deviations (MSD), and the number of degrees of freedom, and a vector of variance components.
#' @importFrom StAMPP stamppAmova stamppNeisD
#' @export
#' @author Bernd Gruber (bugs? Post to https://groups.google.com/d/forum/dartr)

gl.amova <- function(x, nperm=100)
  {
  dd <- stamppNeisD(x, FALSE)
  amova <- stamppAmova(dist.mat = dd, geno = x, perm = nperm)
  return (amova)
}
