#' @name utils.check.datatype
#' @title Utility function to check the class of an object passed to a function.
#'
#' @description 
#' Most functions require access to a genlight object, dist matrix, data matrix or fixed difference list (fd), 
#' and this function checks that a genlight object or one of the above has been passed, whether the genlight object is a SNP dataset 
#' or a SilicoDArT object, and reports back if verbosity is >=2
#'
#' @param x Name of the genlight object, dist matrix, data matrix, glPCA, or fixed difference list (fd) [required]
#' @param accept Vector containing the classes of objects that are to be accepted [default c("genlight","SNP","SilicoDArT"]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity]
#' 
#' @details 
#' This function checks the class of passed object and sets the datatype to "SNP", "SilicoDArT", "dist", "mat",
#' or class[1](x) as appropriate.
#' 
#' Note also that this function checks to see if there are individuals or loci scored as all missing (NA)
#' and if so, issues the user with a warning.
#' 
#' Note: One and only one of gl.check, fd.check, dist.check or mat.check can be TRUE.
#' 
#' @return datatype, "SNP" for SNP data, "SilicoDArT" for P/A data, "dist" for a distance matrix, 
#' "mat" for a data matrix, "glPCA" for an ordination file, or class(x)[1]
#' 
#' @author Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' datatype <- utils.check.datatype(testset.gl)
#' datatype <- utils.check.datatype(as.matrix(testset.gl),accept="matrix")
#' fd <- gl.fixed.diff(testset.gl)
#' datatype <- utils.check.datatype(fd,accept="fd")
#' @export

utils.check.datatype <- function(x,
                                 accept=c("genlight","SNP","SilicoDArT"),
                                 verbose=NULL) {
  
  #### SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

#### CHECK THE TYPE OF OBJECT ####

  if(is(x,"genlight")){
    if(is.null(ploidy(x))){
      stop(error("Fatal Error: ploidy not set in the genlight object, run gl <- gl.compliance.check(gl)\n"))
    }
    if(verbose>=2){cat(report("  Processing genlight object"))}
    if (all(ploidy(x) == 1)){
      if(verbose>=2){
        cat(report(" with Presence/Absence (SilicoDArT) data\n"))
      }
      datatype <- "SilicoDArT"
    } else if (all(ploidy(x) == 2)){
      if(verbose>=2){
        cat(report(" with SNP data\n"))
      }
      datatype <- "SNP"
    } else {
      stop(error("Fatal Error -- SNP or SilicoDArT coding misspecified, run gl <- gl.compliance.check(gl)"))
    }
    # Check for individuals or loci scoring all missing values (NA)
    if(verbose > 1){
    tmp <- gl.filter.allna(x,verbose=0)
    if (nLoc(tmp) < nLoc(x)){
      cat(warn("  Warning: data include loci that are scored NA across all individuals. Consider filtering using gl <- gl.filter.allna(gl)\n"))
    }
    if (nInd(tmp) < nInd(x)){
      cat(warn("  Warning: data include individuals that are scored NA across all loci. Consider filtering using gl <- gl.filter.allna(gl)\n"))
    }
    }
  }
  else if(is(x,"fd")){
      if(is(x$gl,"genlight")){
        # if(is.null(ploidy(x$gl))){
        #   stop(error("Fatal Error: ploidy not set in the genlight object, run gl <- gl.compliance.check(gl)\n"))
        # }
        if(verbose>=2){cat(report("  Processing a fixed difference (fd) object"))}
        if (all(ploidy(x$gl) == 1)){
          if(verbose>=2){
            cat(report(" with Presence/Absence (SilicoDArT) data\n"))
          }
          type <- "SilicoDArT"
        } else if (all(ploidy(x$gl) == 2)){
          if(verbose>=2){
            cat(report(" with SNP data\n"))
          }
          type <- "SNP"
        }
      } else {  
        stop(error("Fatal Error: Fixed Difference object expected! Check format of object\n"))
      } 
      datatype <- "fd"
      # if(verbose>=2 & type=="SilicoDArT"){
      #   cat(report("  Processing a fixed difference (fd) object with Presence/Absence (SilicoDArT) data\n"))
      # }
      # if(verbose>=2 & type=="SNP"){
      #   cat(report("  Processing a fixed difference (fd) object with SNP data\n"))
      # }
  } else if(is(x,"dist")){
    if(verbose>=2){
      cat(report("  Processing a distance matrix\n"))
    }
    datatype <- "dist"
  } else if(is(x,"matrix")){
    if(verbose>=2){
      cat(report("  Processing a data matrix\n"))
    }
    datatype <- "matrix"
  } else if (is(x,"glPca")){
    if(verbose>=2){
      cat(report("  Processing an ordination file (glPca)\n"))
    }
    datatype <- "glPca"
  } else {
    cat(warn("  Warning: Found object of class",class(x)[1],"\n"))
    datatype <- class(x)[1]
  }
  
  #### CHECK WHETHER TO THROW AN ERROR ####
  
  if(!(datatype %in% accept)){
    stop(error("Fatal Error: inappropriate object passed to function, found",datatype,"expecting",paste(accept,collapse=" or ")))
  }
  
  invisible(datatype)
}
