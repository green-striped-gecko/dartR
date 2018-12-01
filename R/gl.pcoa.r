#' PCoA ordination (glPca)
#'
#' This script takes the data on SNP genotypes for individuals and undertakes a Gower PCoa ordination using Euclidean distance and drawing upon
#' data in the original genlight \{adegenet\} object (entity x attribute matrix).
#' The script is essentially a wrapper for glPca() \{adegenet\} with default settings apart from setting parallel=FALSE and
#' converting the eigenvalues to percentages.
#'
#' @param gl Name of the genlight object containing the SNP genotypes by specimen and population [required]
#' @param nfactors Number of dimensions to retain in the output file [default 5]
#' @param parallel TRUE if parallel processing is required (does fail under Windows) [default FALSE]
#' @param n.cores Number of cores to use if parallel processing is requested [default 16]
#' @param v -- verbose if TRUE, silent if FALSE [default TRUE]
#' @import adegenet
#' @return An object of class glPca containing the eigenvalues, factor scores and factor loadings
#' @export
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' pcoa <- gl.pcoa(testset.gl, nfactors=3)

gl.pcoa <- function(gl, nfactors=5, parallel=FALSE, n.cores=16, v=TRUE) {
x <- gl

    if (v==TRUE) {cat("Performing a PCoA, individuals as entities, SNP loci as attributes\n")}
    pcoa <- glPca(x, nf=nfactors, parallel=parallel, n.cores=n.cores)
    e <- pcoa$eig[pcoa$eig > sum(pcoa$eig/length(pcoa$eig))]
    e <- round(e*100/sum(pcoa$eig),1)
    if (v==TRUE) {
      cat(paste("Ordination yielded",length(e),"informative dimensions from",nInd(x)-1,"original dimensions\n"))
      cat(paste("  PCoA Axis 1 explains",e[1],"% of the total variance\n"))
      cat(paste("  PCoA Axis 1 and 2 combined explain",e[1]+e[2],"% of the total variance\n"))
      cat(paste("  PCoA Axis 1-3 combined explain",e[1]+e[2]+e[3],"% of the total variance\n"))
    }
    
    return(pcoa)
}
