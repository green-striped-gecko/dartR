#' PCoA ordination of populations
#'
#' This script takes the data on allele frequencies for populations and undertakes a Gower 
#' PCoA ordination using a nominated distance measure. It draws population information and
#' calculates gene frequencies by drawing upon
#' data in the original genlight \{adegenet\} object (entity x attribute matrix).
#' The script is essentially a wrapper for pcoa() \{ape\}.
#'
#' @param gl -- name of the genlight object containing the SNP genotypes by specimen and population [required]
#' @param c -- Correction methods for negative eigenvalues: \"lingoes\" and \"cailliez\" Refer to \{ape\} documentation. 
#' [default \"none\"]
#' @param method -- the distance measure to be used. This must be one of "euclidean", 
#' "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @return An object of class pcoa containing the eigenvalues, factor scores and factor loadings
#' @author Arthur Georges (gl.bugs@aerg.canberra.edu.au)
#' @export
#' @importFrom ape pcoa
#' @examples
#' \dontrun{
#' pcoa <- gl.pop.pcoa(gl)
#' pcoa <- gl.pop.pcoa(gl, c="cailiez", m="Minkowski", p=2)
#' }

gl.pop.pcoa <- function(gl, c="none", method="euclidean") {

  cat("PCoA on allele frequencies -- populations", 
      "as entities, SNP loci as attributes\n")
  
    mat <- gl.gene.freq(gl, method=levels(pop(gl)))
  
    D <- dist(mat, method=method, diag=FALSE, upper=FALSE, p=2)

    t <- pcoa(D, correction=c, rn=attributes(D)$Labels)

    e <- round(t$values[,2]*100,1)
    cat(paste("Ordination yielded",length(e),"dimensions being the number of populations minus 1\n"))
    cat(paste("  PCoA Axis 1 explains",e[1],"% of the total variance\n"))
    cat(paste("  PCoA Axes 1 and 2 combined explain",e[1]+e[2],"% of the total variance\n"))
    cat(paste("  PCoA Axes 1-3 combined explain",e[1]+e[2]+e[3],"% of the total variance\n"))

    
    t2 <- list(t$values[,1],t$vectors)
    names(t2) <- c("eig", "scores")
    class(t2) <- "glPca"

    return(t2)
}
