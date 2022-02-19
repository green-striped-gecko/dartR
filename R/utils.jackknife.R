#' @name utils.jackknife
#' @title Conducts jackknife resampling using a genlight object
#' @description
#' Jackknife resampling is a statistical procedure where for a dataset of sample 
#' size n, subsamples of size n-1 are used to compute a statistic. The collection 
#' of the values obtained can be used to evaluate the variability around the point 
#' estimate. 
#' 
#' Note that when n is very small, jackknife resampling is not recommended.
#' 
#' @inheritParams gl.drop.ind
#' @param FUN the name of the function to be used to calculate the statistic
#' @param unit The unit to use for resampling. One of c("loc", "ind", "pop"): 
#' loci, individuals or populations
#' @param ... any additional arguments to be passed to FUN
#' @return A list of length n where each element is the output of FUN
#'
#' @author Custodian: Carlo Pacioni -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' @export

utils.jackknife <- function(x, 
                            FUN, 
                            unit="loc", 
                            recalc = FALSE, 
                            mono.rm = FALSE, 
                            n.cores = 1,
                            verbose = NULL, 
                            ...) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING check if packages are installed
  if(!is.character(unit)) stop(error("The argument 'unit' should be character vector"))
  if(length(unit == 1)) {
    if(!unit %in% c("loc", "ind", "pop")) {
      stop(error('The argument "unit" should one of the following: "loc", "ind", "pop"'))
    }
    } else {
      stop(error('The argument "unit" should be of length 1'))
  }
  
  # set variables based on unit
  if(unit == "loc") {
    subsetFUN <- "gl.drop.loc"
    subsetList <- locNames(x)
    } else {
    if(unit == "ind") {
      subsetFUN <- "gl.drop.ind"
      subsetList <- indNames(x)
      } else {
        subsetFUN <- "gl.drop.pop"
        subsetList <- unique(pop(x))
      }
    }
  argmts <- list(...)
  #subsetList.item <- subsetList[1]
  
  # PUT CHECKS HERE for correct unit and x and FUN is a char, and ... is a list
  jacknife <- function(x, jckfun, subsetFUN, subsetList.item, rec = recalc, 
                       mono = mono.rm, opt.argt=argmts) {
    xsub <- do.call(subsetFUN, 
                    args = if(subsetFUN == "gl.drop.loc") {
                      list(x, subsetList.item, verbose=0)
                      } else {
                        list(x, subsetList.item, recal=rec, mono.rm=mono, verbose=0)
                      }
                    )
    oldVerb <- gl.check.verbosity()
    gl.set.verbosity(0)
    on.exit(gl.set.verbosity(oldVerb))
    res <- do.call(jckfun, c(list(x=xsub), opt.argt))
    return(res)
  }
  
  jck <- lapply(subsetList, jacknife, x=x, jckfun=FUN, subsetFUN = subsetFUN)
  return(jck)
  
}
