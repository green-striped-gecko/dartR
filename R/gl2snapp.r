#' Convert a genlight object to nexus format suitable for phylogenetic analysis by SNAPP (via BEAUti)
#'
#' The output nexus file contains the snp data and relevant PAUP command lines suitable for BEAUti.
#' 
#' @references Bryant, D., Bouckaert, R., Felsenstein, J., Rosenberg, N.A. and RoyChoudhury, A. (2012). Inferring species trees directly from biallelic genetic markers: bypassing gene trees in a full coalescent analysis. Molecular Biology and Evolution 29:1917-1932.
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension)[default snapp.nex]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2snapp(testset.gl)

gl2snapp <- function(x, outfile="snapp.nex", outpath=tempdir(), verbose=2) {
  
  # TIDY UP FILE SPECS
  
  outfilespec <- file.path(outpath, outfile)
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
  
  if(!is(x, "genlight")) {
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
  
  if (verbose >= 2) {cat(paste("  Extacting SNP data and creating records for each individual\n"))}
  
  # Extract the reference base and the alternate base for each locus (excuse the contortion)
  m <- as.matrix(x)
  m[is.na(m)] <- "?"
  colnames(m) <- NULL
  df <- data.frame(m)
  df <- cbind(indNames(x),pop(x),df)
  indlabels <- df[,1]
  poplabels <- df[,2]
  
  # Create the snapp file
  if (verbose > 1) {cat(paste("  Writing results to nexus file",outfilespec,"\n"))}
  
  sink(outfilespec)
  
  cat("#nexus\n")
  cat("BEGIN DATA;\n")
  cat(paste0("     dimensions ntax = ",nInd(x)," nchar = ",nLoc(x)," ;\n"))
  cat("     format datatype=integerdata missing=? symbols=\"012\";\n")
  cat("matrix\n")
  for (i in 1:nInd(x)){
    cat(paste0(poplabels[i],"_",indlabels[i]))
    cat("  ")
    cat(m[i,], sep="")
    cat("\n")
  }
  cat(";\n")
  cat("end;\n\n")
  
  sink()
  
  if (verbose > 2) {cat(paste("    Records written to",outfilespec,":",nInd(x),"\n"))}
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(NULL)
  
}