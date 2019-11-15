#' Performs and AMOVA using genlight data. 
#'
#' This script performs an AMOVA based on the genetic distance matrix from stamppNeisD() [package StAMPP] using the amova() function from the package PEGAS for exploring within and between population variation. For detailed information use their help pages: ?pegas::amova, ?StAMPP::stamppAmova. Be aware due to a conflict of the amova functions from various packages I had to "hack" StAMPP::stamppAmova to avoid a namespace conflict.
#' 
#' @param x -- name of the genlight containing the SNP genotypes, with population information [required]
#' @param permutations -- number of permuations to perform for hypothesis testing [default 100]. Please note should be set to 1000 for analysis.]
#' @return An object of class "amova" which is a list with a table of sums of square deviations (SSD), mean square deviations (MSD), and the number of degrees of freedom, and a vector of variance components.
#' @importFrom StAMPP stamppNeisD
#' @importFrom pegas amova
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples 
#' 
#' #permuations should be higher, here set to 10 because of speed
#' gl.amova(bandicoot.gl, permutations=10)
#' 
#' 
gl.amova <- function(x, permutations=100, method="pegas")
{
  dd <- StAMPP::stamppNeisD(x, FALSE)
  
  if (class(x) == "genlight") {
    geno2 <- x
    geno <- as.matrix(geno2)
    sample <- row.names(geno)
    pop.names <- pop(geno2)
    ploidy <- ploidy(geno2)
    geno = geno * (1/ploidy)
    geno[is.na(geno)] = NaN
    format <- vector(length = length(geno[, 1]))
    format[1:length(geno[, 1])] = "genlight"
    pops <- unique(pop.names)
    pop.num <- vector(length = length(geno[, 1]))
    for (i in 1:length(geno[, 1])) {
      pop.num[i] = which(pop.names[i] == pops)
    }
    genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, 
                                   ploidy, format))
    geno <- cbind(genoLHS, geno)
    geno[, 2] = as.character(pop.names)
    geno[, 4] = as.numeric(as.character(geno[, 4]))
    row.names(geno) = NULL
  }
  pop.names <- geno[, 2]
  pop.names <- factor(pop.names)
  temp <- environment(environment)
  assign("dist", dd, envir = temp)
  assign("pop.names", pop.names, envir = temp)
  assign("permutations", permutations, envir = temp)
  res <- with(temp, pegas::amova(dist ~ pop.names, nperm = permutations))
  rm(pop.names, permutations, temp)
  return(res)
  
}
