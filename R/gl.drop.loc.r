#' @name gl.drop.loc
# Preliminaries -- Parameter specifications -------------- 
#' @title Removes specified loci from a dartR genlight object
#' @description
#' This function deletes individuals and their associated metadata. 
#'
#' The script returns a dartR genlight object with the retained loci. 
#' The script works with both genlight objects
#' containing SNP genotypes and Tag P/A data (SilicoDArT).
#'
#' @param x Name of the genlight object [required].
#' @param loc.list A list of loci to be deleted
#' [required, if loc.range not specified].
#' @param first_tmp First of a range of loci to be deleted
#' [required, if loc.list not specified].
#' @param last_tmp Last of a range of loci to be deleted
#' [if not specified, last_tmp locus in the dataset].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress but not results; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#'
#' @export
#' @return A reduced dartR genlight object
#'
#' @family dartR-base
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
# Examples -------------
#' @examples
#' # SNP data
#'   gl2 <- gl.drop.loc(testset.gl, loc.list=c('100051468|42-A/T', '100049816-51-A/G'),verbose=3)
#' # Tag P/A data
#'   gs2 <- gl.drop.loc(testset.gs, loc.list=c('20134188','19249144'),verbose=3)
# See also ------------
#' @seealso \code{\link{gl.keep.loc}} to keep rather than drop specified loci
#'
# End Block --------------
# Function 
gl.drop.loc <- function(x,
                        loc.list = NULL,
                        first_tmp = NULL,
                        last_tmp = NULL,
                        verbose = NULL) {
  # Preliminaries -------------
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
   # Function-specific error checking -----------      
    if (!is.null(loc.list) && !is.null(first_tmp)) {
        flag <- "both"
        if (verbose >= 2) {
            cat(report(
                "  Both a range of loci and a list of loci to keep has been specified\n"
            ))
        }
    } else if (!is.null(loc.list)) {
        flag <- "list"
        if (verbose >= 2) {
            cat(report("  List of loci to drop has been specified\n"))
        }
    } else if (!is.null(first_tmp)) {
        flag <- "range"
        if (verbose >= 2) {
            cat(report("  Range of loci to drop has been specified\n"))
        }
    } else {
        stop(
            error(
                "Fatal Error: Need to specify either a range of loci to drop, or specific loci to drop\n"
            )
        )
    }
    
    if (flag == "both" || flag == "list") {
      
      tmp1 <- loc.list %in% locNames(x)
      tmp2 <- which(tmp1 == FALSE)
      
      if(length(tmp2)>0){
        if(verbose >= 2){
          cat(
            warn(
              "  Warning: Listed loci",
              paste(locNames(x)[tmp2],collapse = " "),
              "not present in the dataset -- ignored\n"
            ))
        }
        
        loc.list <- loc.list[-tmp2]
        
      }
    }
    
    if (flag == "range") {
        if (first_tmp <= 0) {
            cat(warn(
                "  Warning: Lower limit to range of loci cannot be less than 1, set to 1\n)"
            ))
            first_tmp <- 1
        }
        if (first_tmp > nLoc(x)) {
            cat(
                warn(
                    "  Warning: Upper limit to range of loci cannot be greater than the number of loci, set to",
                    nLoc(x),
                    "\n)"
                )
            )
            last_tmp <- nLoc(x)
        }
        if (first_tmp > last_tmp) {
            cat(warn(
                "  Warning: Upper limit is smaller than lower limit, reversed\n"
            ))
            tmp <- first_tmp
            first_tmp <- last_tmp
            last_tmp <- tmp
        }
    }
    
	# DO THE JOB --------------
    # Remove individuals ------
    
    hold <- x
    
    if (verbose >= 2) {
        cat(report("  Deleting the specified loci\n"))
    }
    
    if (!is.null(first_tmp) && !is.null(loc.list)) {
        list.from.range <- locNames(x)[first_tmp:last_tmp]
        loc.list <- unique(c(loc.list, list.from.range))
    } else if (!is.null(first_tmp)) {
        loc.list <- locNames(x)[first_tmp:last_tmp]
    }
    if (length(loc.list) == 0) {
        cat(warn(
            "  Warning: no loci listed to delete! Genlight object returned unchanged\n"
        ))
        x2 <- x
    } else {
        # Remove loci flagged for deletion
      
        x2 <- x[,which(!x$loc.names %in% loc.list)]
        x2@other$loc.metrics <- x@other$loc.metrics[!x$loc.names %in% loc.list,]
    }
    # End block -----------
	
# REPORT A SUMMARY -------------
    # Summary of outcomes --------------    
    if (verbose >= 3) {
        cat("  Summary of recoded dataset\n")
        cat(paste("    Original No. of loci:", nLoc(hold), "\n"))
        cat(paste("    No. of loci deleted:", nLoc(hold) - nLoc(x2), "\n"))
        cat(paste("    No. of loci retained:", nLoc(x2), "\n"))
        # cat(paste(' No. of individuals:', nInd(x2),'\n')) cat(paste(' No. of populations: ', nPop(x2),'\n'))
    }
    
    # ADD TO HISTORY -------------
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END --------------------
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    # End block -------------
    return(x2)
}
