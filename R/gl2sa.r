#' Convert genlight objects to the format used in the SNPassoc package
#'
#' This function exports a genlight object into a SNPassoc object. See \link[SNPassoc]{setupSNP}
#' #' This function needs package SNPassoc. At the time of writing (August 2020) the package was no longer available from CRAN. To install the package check their github repository. \url{https://github.com/isglobal-brge/SNPassoc} and/or use \code{install_github("isglobal-brge/SNPassoc")} to install the function. 
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @export
#' @return Returns an object of class 'snp' to be used with \pkg{SNPassoc}
#' @references Gonzalez, J.R., Armengol, L., Sol?, X., Guin?, E., Mercader, J.M., Estivill, X. and Moreno, V. (2017). SNPassoc: an R package to perform whole genome association studies. Bioinformatics 23:654-655.
#' @author Bernd Guber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \donttest{
#' library(parallel)
#' sa <-gl2sa(testset.gl )
#' }

gl2sa <- function(x, verbose=NULL){
# CHECK IF PACKAGES ARE INSTALLED
  if (!(requireNamespace("parallel", quietly = TRUE))) {
    stop("Package parallel needed for this function to work. Please install it.") }
  if (!(requireNamespace("pegas", quietly = TRUE))) {
    stop("Package pegas needed for this function to work. Please install it.") }
    
  if (!(requireNamespace("SNPassoc", quietly = TRUE))) {
    stop("To use this function you need to install package: SNPassoc. Please refer to the help of the function for instructions (?gl2sa).")
  } else {
# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
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
  sa <- SNPassoc::setupSNP(data.frame(xxx), 1:ncol(xxx), )

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(sa)
  } 
}
