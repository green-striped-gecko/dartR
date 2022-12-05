#' @name gl.filter.ld
#' @title Filters loci based on linkage disequilibrium (LD)
#' @description
#' This function uses the statistic set in the parameter \code{stat_keep} from 
#' function \code{\link{gl.report.ld.map}} to choose the SNP to keep when two 
#' SNPs are in LD. When a SNP is selected to be filtered out in each pairwise 
#' comparison, the function stores its  name in a list. In subsequent pairwise
#'  comparisons, if the SNP is already in the list, the other SNP will be kept.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ld_report Output from function \code{\link{gl.report.ld.map}} 
#' [required].
#' @param threshold Threshold value above which loci will be removed
#' [default 0.2].
#' @param pop.limit Minimum number of populations in which LD should be more
#' than the threshold for a locus to be filtered out.
#' The default value is half of the populations [default ceiling(nPop(x)/2)].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return The reduced genlight object.
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' test <- bandicoot.gl
#' test <- gl.filter.callrate(test,threshold = 1)
#' res <- gl.report.ld.map(test)
#' res_2 <- gl.filter.ld(x=test,ld_report = res)
#' res_3 <- gl.report.ld.map(res_2)
#' }
#' if ((requireNamespace("snpStats", quietly = TRUE)) & (requireNamespace("fields", quietly = TRUE))) {
#' test <- gl.filter.callrate(platypus.gl, threshold = 1)
#' test <- gl.filter.monomorphs(test)
#' test <- test[,1:250]
#' report <- gl.report.ld.map(test)
#' }
#' @seealso \code{\link{gl.report.ld.map}}
#' @family filter functions
#' @export

gl.filter.ld <- function(x,
                         ld_report,
                         threshold = 0.2,
                         pop.limit = ceiling(nPop(x) / 2),
                         verbose = NULL) {
  
  x_hold <- x
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # Check monomorphs have been removed up to date
  if (x@other$loc.metrics.flags$monomorphs == FALSE) {
    if (verbose >= 2) {
      cat(
        warn(
          "  Warning: Data may include monomorphic loci in call rate 
          calculations for filtering\n"
        )
      )
    }
  }
  
  x <- gl.keep.pop(x,pop.list = as.character(unique(ld_report$pop)),verbose = 0)
  
  ld_tmp <- ld_report[ld_report$ld_stat >= threshold, ]
  ld_tmp$test_stat <- ld_tmp$locus_a.stat_keep >= ld_tmp$locus_b.stat_keep
  ld_tmp$pop <- as.factor(ld_tmp$pop)
  ld_tmp_pop <- split(ld_tmp, f = ld_tmp$pop)
  
  loci_list <- vector(mode = "list", length = length(ld_tmp_pop))
  
  for (i in 1:length(ld_tmp_pop)) {
    ld_pop <- ld_tmp_pop[[i]]
    for (y in 1:nrow(ld_pop)) {
      if (ld_pop[y, "test_stat"] == TRUE) {
        loci_tmp <- ld_pop[y, "locus_b.snp.name"]
      } else{
        loci_tmp <- ld_pop[y, "locus_a.snp.name"]
      }
      
      if (loci_tmp %in% loci_list[[i]]) {
        next
      } else{
        loci_list[[i]] <- c(loci_list[[i]], loci_tmp)
      }
      
    }
  }
  
  loci_list_res <- Reduce("c", loci_list)
  loci_names_tmp <- names(table(loci_list_res))
  loci_names <- loci_names_tmp[table(loci_list_res) >= pop.limit]
  
  x2 <- gl.drop.loc(x_hold, loc.list =  loci_names, verbose = 0)
  
  # REPORT A SUMMARY
  if (verbose >= 2) {
    cat("  Summary of filtered dataset\n")
    cat(paste("    LD for loci >", threshold, "\n"))
    cat(paste("    Original No. of loci :", nLoc(x), "\n"))
    cat(paste("    No. of loci retained:", nLoc(x2), "\n"))
    cat(paste("    No. of populations: ", nPop(x2), "\n"))
  }
  
  # ADD TO HISTORY
  
  nh <- length(x2@other$history)
  x2@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(invisible(x2))
  
}
