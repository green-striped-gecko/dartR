#' Convert a genlight object to a treemix input file
#' 
#' The output file contains the snp data in the format expected by treemix -- see the treemix manual. The file will be gzipped before in order to be recognised
#' by treemix. Plotting functions provided with treemix can be obtained using source("C:/workspace/R/dartR/R/gl.treemix.plotting_funcs.r").  
#'
#' Reference: Pickrell and Pritchard (2012). Inference of population splits and mixtures from genome-wide allele frequency data. PLoS Genetics https://doi.org/10.1371/journal.pgen.1002967
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including gz extension).
#' @param outpath -- path where to save the output file (set to tempdir by default)
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2treemix(testset.gl, outpath=getwd())

gl2treemix <- function(x, outfile="treemix_input.gz", outpath=tempdir(), v=2) {
  
  # ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.report.callrate!\n"); stop()
  }
  if (v < 0 | v > 5){
    cat("    Warning: verbosity must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    v <- 2
  }
  
  if (v > 0) {cat(paste("Starting gl2treemix: Create treemix input file\n\n"))}

# Calculate allele counts
  
  freq <- gl.percent.freq(x)
  freq$ref <- freq$nobs*2-freq$sum
  freq$alt <- freq$sum
  freq$sum <- NULL
  freq$nobs <- NULL
  freq$nmissing <- NULL
  freq$frequency <- NULL
  freq$n <- NULL

# Output the file
  outfile <- file.path(outpath, outfile)
  if (v > 1) {cat(paste("    Writing results to treemix input file",outfile,"\n"))}
  sink(gzfile(outfile))

  cat (unique(as.character(freq$popn)),"\n")
  k <- 1
  for (j in 1:nLoc(x)) {
  for (i in k:(k+nPop(x)-1)) {
      cat(paste0(freq$ref[i],",",freq$alt[i])," ")
  }
  cat("\n")
  k <- k + nPop(x)
  }
  
  sink()
  
  if (v > 2) {cat(paste("    Records written to",outfile,":",nInd(x),"\n"))}
  if (v > 0) {
    cat("gl2treemix Completed\n")
    cat("Output file has been gzipped for input to treemix\n")
  }
  
  return(NULL)

}


