#' Performs and AMOVA using genlight data. 
#'
#' This script performs an AMOVA based on the genetic distance matrix from stamppNeisD() [package StAMPP] using the amova() function from the package PEGAS for exploring within and between population variation. For detailed information use their help pages: ?pegas::amova, ?StAMPP::stamppAmova
#' 
#' @param x -- name of the genlight containing the SNP genotypes, with population information [required]
#' @param permutations -- number of permuations to perform for hypothesis testing [default 100]. Please note should be set to 1000 for analysis.]
#' @return An object of class "amova" which is a list with a table of sums of square deviations (SSD), mean square deviations (MSD), and the number of degrees of freedom, and a vector of variance components.
#' @importFrom StAMPP stamppNeisD stamppAmova
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples 
#' #permuations should be higher, here set to 10 because of speed
#' gl.amova(bandicoot.gl, permutations=10)
#' 
#' 
gl.amova <- function(x, permutations=100)
{
  dd <- stamppNeisD(x, FALSE)
  amova <- StAMPP::stamppAmova(dist.mat = dd, geno = x, perm = permutations)
  return (amova)
}
