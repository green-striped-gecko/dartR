#' @name gl.filter.allna
#' @title Remove loci that are all NA across individuals and/or individuals with all NA across loci
#' 
#' @description
#' This script deletes deletes loci or individuals with all calls missing (NA), from a genlight object
#'
#' A DArT dataset will not have loci for which the calls are scored all as missing (NA) for a particular individual, 
#' but such loci can arise rarely when populations or individuals are deleted. Similarly, a DArT dataset will 
#' not have individuals for which the calls are scored all as missing (NA) across all loci, 
#' but such individuals may sneak in to the dataset when loci are deleted.
#' Retaining indiviudal or loci with all NAs can cause issues for several functions. 
#' 
#' @param x Name of the input genlight object [required]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 
#' 3, progress and results summary; 5, full report [default 2, unless specified using gl.set.verbosity]
#' @return A genlight object having removed individuals that are scored NA across all loci, or loci that are scored NA across all individuals
#' 
#' @author Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SNP data
#'   result <- gl.filter.allna(testset.gl, verbose=3)
#' # Tag P/A data
#'   result <- gl.filter.allna(testset.gs, verbose=3)
#'  
#' @family filters and filter reports
#' @import utils patchwork
#' @export

gl.filter.allna <- function (x,verbose=NULL) {

# SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
# FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jody",verbosity=verbose)
  
# CHECK DATATYPE 
#  datatype <- utils.check.datatype(x,verbose=verbose) # recurrence clash
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
      stop(error("Fatal Error -- SNP or SilicoDArT coding misspecified, run gl <- gl.compliance.check(gl)."))
    }
  }
  
# DO THE JOB

  if (verbose >= 2){
    cat(report("  Identifying and removing loci and individuals scored all missing (NA)\n"))
  }  
  
  # Consider loci
  na.counter <- 0
  loc.list <- array(NA,nLoc(x))
  nL <- nLoc(x)
  matrix <- as.matrix(x)
  l.names <- locNames(x)
  for (i in 1:nL){
    row <- matrix[,i] # Row for each locus
    if (all(is.na(row))){
      loc.list[i] <- l.names[i]
      if (all(is.na(row))){
        na.counter = na.counter + 1
      }
    }
  }
  if (na.counter == 0){
    if(verbose >= 3){cat("  Zero loci that are missing (NA) across all individuals\n")}
  } else {
    loc.list <- loc.list[!is.na(loc.list)]
    if(verbose >= 3){cat("  Loci that are missing (NA) across all individuals:",paste(loc.list,collapse=", "),"\n")}
    x2 <- x[,!x$loc.names%in%loc.list]
    x2@other$loc.metrics <- x@other$loc.metrics[!x$loc.names%in%loc.list,]
    x <- x2
    if(verbose >= 2){cat("  Deleted\n")}
  }
  
  # Consider individuals
  na.counter <- 0
  ind.list <- array(NA,nInd(x))
  nI <- nInd(x)
  matrix <- as.matrix(x)
  i.names <- indNames(x)
  for (i in 1:nI){
    col <- matrix[i,] # Row for each locus
    if (all(is.na(col))){
      ind.list[i] <- i.names[i]
      if (all(is.na(col))){
        na.counter = na.counter + 1
      }
    }
  }
  if (na.counter == 0){
    if(verbose >= 3){cat("  Zero individuals that are missing (NA) across all loci\n")}
  } else {
    ind.list <- ind.list[!is.na(ind.list)]
    if(verbose >= 3){cat("  Individuals that are missing (NA) across all loci:",paste(ind.list,collapse=", "),"\n")}
    x <- x[!x$ind.names%in%ind.list]
    if(verbose >= 2){cat("  Deleted\n")}
  }

# ADD TO HISTORY
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()
  
# FLAG SCRIPT END
  if (verbose >= 1){
    cat(report("Completed:",funname,"\n"))
  }  
  
return (x)
}
