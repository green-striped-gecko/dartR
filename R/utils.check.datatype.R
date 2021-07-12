#' @name utils.check.datatype
#'
#' @title Utility function to check that the object passed to the function is a genlight object.
#'
#' @description 
#' Most functions require access to a genlight object, and this function checks that a genlight object has been passed,
#' whether it is a SNP dataset or a SilicoDArT object, and reports back if verbosity is >=2
#'
#' @param x Name of the genlight object containing the SNP data or tag presence/absence data (SilicoDArT) [required]
#' @param strict If TRUE, a fatal error occurs if the passed object is not a genlight object [default TRUE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity]
#' 
#' @details 
#' This function checks whether the input object is a genlight object and sets the datatype to "SNP" or "SilicoDArT" as appropriate.  
#' 
#' @return datatype, "SNP" for SNP data or  "SilicoDArT" for P/A data, or class(x)[1]
#' @examples
#' datatype <- utils.check.datatype(testset.gl)
#' @export

utils.check.datatype <- function(x,
                                 strict=TRUE,
                                 verbose=NULL) {
  
  #### SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  #### CHECK FOR GENLIGHT OBJECT ####
  if (!is(x, "genlight")) {
    if(strict){
      stop(error("  Warning: genlight object expected! Found",class(x)[1],".\n"))
    } else {
      cat(warn("  Warning: genlight object expected! Found",class(x)[1],".\n"))
      datatype <- class(x)[1]
    }  
  } 
  
  #### CHECK AND REPORT DATA TYPE ####
  if(is(x,"genlight")){
  if (all(x@ploidy == 1)){
    if(verbose>=2){
      cat(report("  Processing Presence/Absence (SilicoDArT) data\n"))
      }
    datatype <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if(verbose>=2){
      cat(report("  Processing SNP data\n"))
      }
    datatype <- "SNP"
  } else {
    stop (error("Fatal Error: datatype must be fragment P/A or SNP data)"))
  }
  }
  
  invisible(datatype)
  
  }
