#' A utility function to check that the object passed to the function is a genlight object.
#'
#' Most functions require access to a genlight object, and this function checks that a genlight object has been passed,
#' whether it is a SNP dataset or a SilicoDArT object, and reports back if verbosity is >=2
#'
#' @param x -- name of the genlight object containing the SNP data or tag presence/absence data (SilicoDArT) [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return datatype, "SNP" for SNP data or  "SilicoDArT" for P/A data
#' @export
#' @examples
#' utils.check.gl(testset.gl,env=environment())

utils.check.gl <- function(x=NULL,verbose=options()$dartR_verbose,env) {
  
#### SET VERBOSITY
  env$verbose <- gl.check.verbosity(verbose)

#### CHECK FOR GENLIGHT OBJECT ####
  if (class(x) != "genlight") {
    stop(error("Fatal Error: genlight object required!"))
  }
  
#### CHECK AND REPORT DATA TYPE ####
  if (all(x@ploidy == 1)){
    if(verbose>=2){cat(report("  Processing Presence/Absence (SilicoDArT) data\n"))}
    env$datatype <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if(verbose>=2){cat(report("  Processing SNP data\n"))}
    env$datatype <- "SNP"
  } else {
    stop (error("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)"))
  }
  
  # Work around a bug in adegenet if genlight object is created by subsetting
  x@other$loc.metrics <- as.data.frame(x@other$loc.metrics)
  
  
}
