#' Converts genlight objects to hiphop format
#'
#' This function exports genlight objects to the format used by the parentage assignment R package hiphop. Hiphop can be used for paternity and maternity assignment and outperforms conventional methods where closely related individuals occur in the pool of possible parents. The method compares the genotypes of offspring with any combination of potentials parents and scores the number of mismatches of these individuals at bi-allelic genetic markers (e.g. Single Nucleotide Polymorphisms).
#' 
#' @references Cockburn, A., Penalba, J.V.,Jaccoud, D.,Kilian, A., Brouwer, L.,Double, M.C., Margraf, N., Osmond, H.L., van de Pol, M. and Kruuk, L.E.B.(in revision). HIPHOP: improved paternity assignment among close relatives using a simple exclusion method for bi-allelic markers. Molecular Ecology Resources, DOI to be added upon acceptance
#' @param gl -- name of the genlight object containing the SNP data
#' @param probar -- if TRUE, a progress bar will be displayed for long loops [default = FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return Dataframe containing all the genotyped individuals (offspring and potential parents) and their genotypes scored using bi-allelic markers.
#' @export
#' @author Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' @examples
#' \donttest{
#' result <- gl2hiphop(testset.gl)
#' }

gl2hiphop <- function(gl, probar = FALSE, verbose = NULL){
    
    # TRAP COMMAND, SET VERSION
      funname <- match.call()[[1]]
      build <- "Jacob"
    
    # SET VERBOSITY
    
    if (is.null(verbose) & !is.null(gl@other$verbose)) 
        verbose = gl@other$verbose
    if (is.null(verbose)) { 
        verbose = 2
    }
    if (verbose < 0 | verbose > 5) {
        cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
        verbose <- 2
    }
    if (verbose > 0) {
        cat("Start conversion....\n")
        ptm <- proc.time()[3]
        cat("Once finished, we recommend to save the object using >save(object, file=\"object.rdata\")\n")
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
  
  if(class(gl)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }
  
  if (verbose >= 2){
    if (all(gl@ploidy == 1)){
      stop("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n")
    } else if (all(gl@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }
  
      # DO THE JOB
    
    x <- as.matrix(gl[, ])
    x[x==1] <- "het"
    x[x==2] <- 1
    x[x=="het"] <- 2
    x <- as.data.frame(x)
    x <- x %>% dplyr::mutate_all(as.numeric)
   
    if (probar) {
        pb <- txtProgressBar(min = 0, max = 1, style = 3, initial = NA)
    }
        if (probar) {
            setTxtProgressBar(pb, i/nrow(x))
        }
    if (probar) {
        close(pb)
    }
    
     # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

    return(x)
}






