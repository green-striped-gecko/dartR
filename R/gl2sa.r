#' Convert genlight objects to the format used in the SNPassoc package
#'
#' This function exports a genlight object into a SNPassoc object. See \link[SNPassoc]{setupSNP}
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @export
#' @return Returns an object of class 'snp' to be used with \pkg{SNPassoc}
#' @references Gonz?lez, J.R., Armengol, L., Sol?, X., Guin?, E., Mercader, J.M., Estivill, X. and Moreno, V. (2017). SNPassoc: an R package to perform whole genome association studies. Bioinformatics 23:654-655.
#' @author Bernd Guber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2plink(testset.gl)

gl2sa <- function(x, verbose=NULL){

# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START
  # set verbosity
  if (is.null(verbose) & !is.null(x@other$verbose)) verbose=x@other$verbose
  if (is.null(verbose)) verbose=2
 

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }

  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  
  if (verbose >= 2){
    if (all(x@ploidy == 1)){
      stop("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }

# DO THE JOB

  if (verbose >= 2){ cat("  Writing data to SNPassoc object\n")}
  pop <- gl2gi(x)
  xxx <- pegas::as.loci(pop)[,-1]
  sa <- SNPassoc::setupSNP(data.frame(xxx), 1:ncol(xxx))

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(sa)
}
