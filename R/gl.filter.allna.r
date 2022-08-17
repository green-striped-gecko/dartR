#' @name gl.filter.allna
#' @title Filters loci that are all NA across individuals and/or populations 
#' with all NA across loci
#'
#' @description
#' This script deletes deletes loci or individuals with all calls missing (NA),
#'  from a genlight object
#'
#' A DArT dataset will not have loci for which the calls are scored all as
#' missing (NA) for a particular individual, but such loci can arise rarely when
#'  populations or individuals are deleted. Similarly, a DArT dataset will not
#'  have individuals for which the calls are scored all as missing (NA) across
#'  all loci, but such individuals may sneak in to the dataset when loci are
#'  deleted. Retaining individual or loci with all NAs can cause issues for
#'  several functions.
#'  
#'  Also, on occasion an analysis will require that there are some loci scored
#'  in each population. Setting by.pop=TRUE will result in removal of loci when
#'  they are all missing in any one population.
#'  
#' Note that loci that are missing for all individuals in a population are
#' not imputed with method 'frequency' or 'HW'. Consider 
#' using the function \code{\link{gl.filter.allna}} with by.pop=TRUE.
#'
#' @param x Name of the input genlight object [required].
#' @param by.pop If TRUE, loci that are all missing in any one population
#' are deleted [default FALSE]
#' @param recalc Recalculate the locus metadata statistics if any individuals
#' are deleted in the filtering [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return A genlight object having removed individuals that are scored NA
#' across all loci, or loci that are scored NA across all individuals.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SNP data
#'   result <- gl.filter.allna(testset.gl, verbose=3)
#' # Tag P/A data
#'   result <- gl.filter.allna(testset.gs, verbose=3)
#'
#' @family filter functions
#' @import utils patchwork
#' @export

gl.filter.allna <- function(x,
                            by.pop = FALSE,
                            recalc = FALSE,
                            verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Josh",
                     verbosity = verbose)
    
    # recurrence clash
    # CHECK DATATYPE datatype <- utils.check.datatype(x,verbose=verbose) 
    if (is(x, "genlight")) {
        if (is.null(ploidy(x))) {
            stop(
                error(
                    "Fatal Error: ploidy not set in the genlight object, run 
                    gl <- gl.compliance.check(gl)\n"
                )
            )
        }
        if (verbose >= 2) {
            cat(report("  Processing genlight object"))
        }
        if (all(ploidy(x) == 1)) {
            if (verbose >= 2) {
                cat(report(" with Presence/Absence (SilicoDArT) data\n"))
            }
            datatype <- "SilicoDArT"
        } else if (all(ploidy(x) == 2)) {
            if (verbose >= 2) {
                cat(report(" with SNP data\n"))
            }
            datatype <- "SNP"
        } else {
            stop(
                error(
                    "Fatal Error -- SNP or SilicoDArT coding misspecified, run 
                    gl <- gl.compliance.check(gl)."
                )
            )
        }
    }
    
    # DO THE JOB
    
    if (verbose >= 2) {
        if (by.pop==FALSE){
        cat(
            report(
                "  Identifying and removing loci and individuals scored all 
                missing (NA)\n"
            )
        )
        } else {
            cat(
                report(
                    "  Identifying and removing loci that are all missing (NA) 
                    in any one individual\n"
                )
            )
        }
    }
    
    if (by.pop==FALSE){
    # Consider loci
    if(verbose >=2){
        cat(report("  Deleting loci that are scored as all missing (NA)\n"))
    }
    na.counter <- 0
    loc.list <- array(NA, nLoc(x))
    nL <- nLoc(x)
    matrix <- as.matrix(x)
    l.names <- locNames(x)
    for (i in 1:nL) {
        row <- matrix[, i]  # Row for each locus
        if (all(is.na(row))) {
            loc.list[i] <- l.names[i]
            if (all(is.na(row))) {
                na.counter <-na.counter + 1
            }
        }
    }
    if (na.counter == 0) {
        if (verbose >= 3) {
            cat("  Zero loci that are missing (NA) across all individuals\n")
        }
    } else {
        loc.list <- loc.list[!is.na(loc.list)]
        if (verbose >= 3) {
            cat(
                "  Loci that are missing (NA) across all individuals:",
                paste(loc.list, collapse = ", "),
                "\n"
            )
        }
        x2 <- x[, !x$loc.names %in% loc.list]
    x2@other$loc.metrics <- x@other$loc.metrics[!x$loc.names %in% loc.list, ]
        x <- x2
        if (verbose >= 2) {
            cat("  Deleted\n")
        }
    }
    
    # Consider individuals
    if(verbose >=2){
    cat(report("  Deleting individuals that are scored as all missing (NA)\n"))
    }
    na.counter <- 0
    ind.list <- array(NA, nInd(x))
    nI <- nInd(x)
    matrix <- as.matrix(x)
    i.names <- indNames(x)
    for (i in 1:nI) {
        col <- matrix[i, ]  # Row for each locus
        if (all(is.na(col))) {
            ind.list[i] <- i.names[i]
            if (all(is.na(col))) {
                na.counter <-na.counter + 1
            }
        }
    }
    if (na.counter == 0) {
        if (verbose >= 3) {
            cat("  Zero individuals that are missing (NA) across all loci\n")
        }
    } else {
        ind.list <- ind.list[!is.na(ind.list)]
        if (verbose >= 3) {
            cat(
                "  Individuals that are missing (NA) across all loci:",
                paste(ind.list, collapse = ", "),
                "\n"
            )
        }
        x <- x[!x$ind.names %in% ind.list]
        if (verbose >= 2) {
            cat("  Deleted\n")
        }
    }
    }
    
    if (by.pop==TRUE){
        if(verbose >=2){
cat(report("  Deleting loci that are all missing (NA) in any one population\n"))
        }
        total <- 0
        loc.list <- NULL
        for (i in 1:nPop(x)){
            tmpop <- as.matrix(gl.keep.pop(x,popNames(x)[i],verbose=0))
            tmpsums <- apply(tmpop,2,function(x){all(is.na(x))})
            # tmpsums <-  colSums(tmpop)
            tmp.list <- locNames(x)[tmpsums==TRUE]
            # tmp.list <- locNames(x)[is.na(tmpsums)]
            count <- length(tmp.list)
            if (verbose >= 3){
               cat("    ",popNames(x)[i],": deleted",count,"loci\n")
            }  
            total <- total + count
            loc.list <- c(loc.list,tmp.list)
        }
        loc.list <- unique(loc.list)
        if (verbose >= 3){
cat("\n  Loci all NA in one or more populations:",length(loc.list),
    "deleted\n\n")
        }  
        x <- gl.drop.loc(x,loc.list=loc.list,verbose=0)
    }
    
    if (recalc) {
        # Recalculate all metrics, including Call Rate (flags reset in utils 
      #scripts)
        x <- gl.recalc.metrics(x, verbose = verbose)
    } else {
        # Reset the flags as FALSE for all metrics except all na (dealt with 
      #elsewhere)
        x@other$loc.metrics.flags$AvgPIC <- FALSE
        x@other$loc.metrics.flags$OneRatioRef <- FALSE
        x@other$loc.metrics.flags$OneRatioSnp <- FALSE
        x@other$loc.metrics.flags$PICRef <- FALSE
        x@other$loc.metrics.flags$PICSnp <- FALSE
        x@other$loc.metrics.flags$maf <- FALSE
        x@other$loc.metrics.flags$FreqHets <- FALSE
        x@other$loc.metrics.flags$FreqHomRef <- FALSE
        x@other$loc.metrics.flags$FreqHomSnp <- FALSE
        x@other$loc.metrics.flags$CallRate <- FALSE
    }
    
    # ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x)
}
