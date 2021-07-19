#' @name utils.check.datatype
#' @title Utility function to check that the object passed to the function is a genlight object.
#'
#' @description 
#' Most functions require access to a genlight object, and this function checks that a genlight object has been passed,
#' whether it is a SNP dataset or a SilicoDArT object, and reports back if verbosity is >=2
#'
#' @param x Name of the genlight object containing the SNP data or tag presence/absence data (SilicoDArT) [required]
#' @param gl.check If TRUE, accepts a genlight object [default TRUE]
#' @param fd.check If TRUE, accepts a fixed difference object, class fd [default FALSE]
#' @param dist.check If TRUE, accepts a distance matrix, class dist [default FALSE]
#' @param mat.check If TRUE, accepts a square matrix, class matrix [default FALSE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity]
#' 
#' @details 
#' This function checks the class of passed object and sets the datatype to "SNP", "SilicoDArT", "dist", "mat",
#' or class[1](x) as appropriate.
#' 
#' Note: One and only one of gl.check, fd.check, dist.check or mat.check can be TRUE.
#' 
#' @return datatype, "SNP" for SNP data, "SilicoDArT" for P/A data, "dist" for a distance matrix, 
#' "mat" for a data matrix or class(x)[1]
#' @examples
#' datatype <- utils.check.datatype(testset.gl)
#' datatype <- utils.check.datatype(as.matrix(testset.gl),mat.check=TRUE)
#' @export

utils.check.datatype <- function(x,
                                 gl.check=FALSE,
                                 fd.check=FALSE,
                                 dist.check=FALSE,
                                 mat.check=FALSE,
                                 verbose=NULL) {
  
#### SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
#### SCRIPT SPECIFIC CHECKING
  if(sum(gl.check,fd.check,dist.check,mat.check)==0){
    # if(verbose >=2){
    #   cat(warn("  Warning: One of gl.check, fd.check, dist.check or mat.check must be TRUE\n"))
    #   cat(warn("    Set to expect a genlight object\n"))
    # }
    gl.check <- TRUE # Setting the default here
  }
  if(sum(gl.check,fd.check,dist.check,mat.check)>1){
    stop(error("Fatal Error: Only one of gl.check, fd.check, dist.check or mat.check can be TRUE"))
  }
  
#### CHECK OBJECT ####
  
# Genlight object expected
  if(is(x,"genlight") & gl.check){
    if(verbose>=2){cat(report("  Processing genlight object\n"))}
    if (all(x@ploidy == 1)){
      if(verbose>=2){
        cat(report("    with Presence/Absence (SilicoDArT) data\n"))
      }
      datatype <- "SilicoDArT"
    } else if (all(x@ploidy == 2)){
      if(verbose>=2){
        cat(report("    with SNP data\n"))
      }
      datatype <- "SNP"
    } else {
      stop(error("Fatal Error -- SNP or SilicoDArT coding misspecified, run gl <- gl.compliance.check(gl)."))
    }
  }
  
# Other objects
  if (!is(x, "genlight")) {
    if(gl.check){
      stop(error("  Warning: genlight object expected! Found",class(x)[1],".\n"))
    } 
    if(fd.check & datatype=="SilicoDArT"){
      if(!exists("x$fd")){stop(error("Fatal Error: Fixed Difference object expected! Check format of object\n"))}
      if(verbose>=2){
        cat(report("  Processing a fixed difference (fd) object with Presence/Absence (SilicoDArT) data\n"))
      }
      datatype <- "fd"
    } else if(fd.check & datatype=="SNP"){
      if(!exists("x$fd")){stop(error("Fatal Error: Fixed Difference object expected! Check format of object\n"))}
      if(verbose>=2){
        cat(report("  Processing a fixed difference (fd) object with SNP data\n"))
      }
      datatype <- "fd"
    } else if(dist.check){
      if(!is(x,"dist")){stop(error("Fatal Error: Distance matrix expected! Found ",class(x)[1],"object.\n"))}
      if(verbose>=2){
        cat(report("  Processing a distance matrix\n"))
      }
      datatype <- "dist"
    } else if(mat.check){
      if(!is(x,"matrix")){stop(error("Fatal Error: Data matrix expected! Found a",class(x)[1],"object.\n"))}
      if(verbose>=2){
        cat(report("  Processing a data matrix\n"))
      }
      datatype <- "matrix"
    } else {
      cat(warn("  Warning: Found object of class",class(x)[1],"\n"))
      datatype <- class(x)[1]
    } 
  } 
  
  invisible(datatype)
}
