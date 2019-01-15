#' Convert a genlight object to format suitable for input to Bayescan
#'
#' The output text file contains the snp data and relevant BAyescan command lines to guide input.
#' 
#' @reference: Foll M and OE Gaggiotti (2008) A genome scan method to identify selected loci appropriate for both dominant and codominant markers: A Bayesian perspective. Genetics 180: 977-993.
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension).
#' @param outpath -- path where to save the output file [default tempdir()]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2bayescan(testset.gl)

gl2bayescan <- function(x, outfile="bayescan.txt", outpath=tempdir(), v=2) {
  
  if (v > 0) {cat(paste("Starting gl2bayescan: Create text input file\n\n"))}
  if (v > 1) {cat(paste("Extacting SNP data and creating records for each individual\n"))}

# Create the bayescan file
  outfile <- file.path(outpath, outfile)
  
# Prepare the data  
  mat <- gl.percent.freq(x, v=v)
  mat <- mat[order(mat$popn),]

# Create the bayescan input file  
  if (v > 1) {cat(paste("Writing text input file for Bayescan",outfile,"\n"))}
  sink(outfile)
  
  cat(paste0("[loci]=",nLoc(x)),"\n\n")
  cat(paste0("[populations]=",nPop(x)),"\n\n")

  for (i in 1:nPop(x)){
    cat(paste0("[pop]=",i),"\n")
    popi <- mat[mat$popn==mat$popn[i],]
    for (j in 1:length(popi$popn)){
      cat(j,(2*popi$nobs[j]),2,popi$sum[j],(2*popi$nobs[j]-popi$sum[j]),"\n")
    }
    cat("\n")
  }
  
  sink()
  
  # Finish up
  
  if (v >= 3) {cat(paste("Records written to",outfile,"\n"))}
  if (v > 0) {cat("gl2bayescan Completed\n")}
  
  return(NULL)

}


