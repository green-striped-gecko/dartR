#' PCoA ordination (glPca)
#'
#' This script takes the data on SNP genotypes for individuals and undertakes a Gower PCoa ordination using Euclidean distance and drawing upon
#' data in the original genlight \{adegenet\} object (entity x attribute matrix).
#' The script is essentially a wrapper for glPca() \{adegenet\} with default settings apart from setting parallel=FALSE and
#' converting the eigenvalues to percentages.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param nfactors Number of dimensions to retain in the output file [default 5]
#' @param parallel TRUE if parallel processing is required (does fail under Windows) [default FALSE]
#' @param n.cores Number of cores to use if parallel processing is requested [default 16]
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return An object of class glPca containing the eigenvalues, factor scores and factor loadings
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' pcoa <- gl.pcoa(testset.gl, nfactors=3)

# Last amended 3-Feb-19

gl.pcoa <- function(x, nfactors=5, parallel=FALSE, n.cores=16, verbose=2) {

# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }


  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}


# DO THE JOB

    if (verbose >=2) {cat("  Performing a PCoA, individuals as entities, SNP loci as attributes\n")}
    pcoa <- glPca(x, nf=nfactors, parallel=parallel, n.cores=n.cores)
    e <- pcoa$eig[pcoa$eig > sum(pcoa$eig/length(pcoa$eig))]
    e <- round(e*100/sum(pcoa$eig),1)
    if (verbose >=3) {
      cat(paste("  Ordination yielded",length(e),"informative dimensions from",nInd(x)-1,"original dimensions\n"))
      cat(paste("    PCoA Axis 1 explains",e[1],"% of the total variance\n"))
      cat(paste("    PCoA Axis 1 and 2 combined explain",e[1]+e[2],"% of the total variance\n"))
      cat(paste("    PCoA Axis 1-3 combined explain",e[1]+e[2]+e[3],"% of the total variance\n"))
    }

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
    return(pcoa)
}
