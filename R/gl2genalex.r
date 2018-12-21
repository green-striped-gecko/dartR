#' Convert a genlight object to format suitable for input to genalex
#'
#' The output csv file contains the snp data and other relevant lines suitable for genalex. This script is a wrapper for genind2genalex {poppr}
#' 
#' Reference: Peakall, R. and Smouse P.E. (2012) GenAlEx 6.5: genetic analysis in Excel. Population genetic software for teaching and research-an update. Bioinformatics 28, 2537-2539.
#' http://bioinformatics.oxfordjournals.org/content/28/19/2537
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension).
#' @param outpath -- path where to save the output file [default tempdir()]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @import from poppr genind2genalex
#' @return NULL
#' @export
#' @author Katrin Hohwieler (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2genalex(testset.gl,outfile="testset.csv)

gl2genalex <- function(x, outfile="genaxex.csv", outpath=tempdir(), v=2) {
  
# ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.report.callrate!\n"); stop()
  }
 
  if (v < 0 | v > 5){
    cat("    Warning: verbosity must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    v <- 2
  }
  
  # FLAG SCRIPT START
  
  if (v > 0) {cat(paste("Starting gl2genalex: Create genalex input file\n\n"))}
  
  gind <- gl2gi(x, v=0)
  poppr::genind2genalex(gind, filename = outfile, sequence = TRUE, overwrite = FALSE)
  

  if (v > 2) {cat(paste("    Records written to",outfile,":",nInd(x),"\n"))}
  
# FLAG SCRIPT END  
  if (v > 0) {cat("\ngl2genalex Completed\n")}
  
  return(NULL)

}


