#' @name gl.filter.parent.offspring
#'
#' @title Filter putative parent offspring within a population
#'
#' @description 
#' This script removes individuals suspected of being related as parent-offspring,
#' using the output of the function \code{\link{gl.report.parent.offspring}}, which
#' examines the frequency of pedigree inconsistent loci, that is,
#' those loci that are homozygotes in the parent for the reference allele, and
#' homozygous in the offspring for the alternate allele. This condition is not
#' consistent with any pedigree, regardless of the (unknown) genotype of the other
#' parent. The pedigree inconsistent loci are counted as an indication of whether
#' or not it is reasonable to propose the two individuals are in a parent-offspring
#' relationship.
#'
#' @param x Name of the genlight object containing the SNP genotypes [required]
#' @param min.rdepth Minimum read depth to include in analysis [default = 12]
#' @param min.reproducibility Minimum reproducibility to include in analysis [default = 1]
#' @param range Specifies the range to extend beyond the interquartile range for delimiting outliers [default = 1.5 interquartile ranges]
#' @param rm.monomorphs If TRUE, remove monomorphic loci after filtering individuals [default FALSE].
#' @param plot_theme Theme for the plot. See Details for options [default theme_dartR()].
#' @param plot_colours List of two color names for the borders and fill of the
#'  plots [default two_colors].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity]
#'
#' @details 
#' If two individuals are in a parent offspring relationship, the true
#' number of pedigree inconsistent loci should be zero, but SNP calling is not 
#' infallible. Some loci will be miss-called. The problem thus becomes one of determining
#' if the two focal individuals have a count of pedigree inconsistent loci less than
#' would be expected of typical unrelated individuals. There are some quite sophisticated
#' software packages available to formally apply likelihoods to the decision, but we
#' use a simple outlier comparison.
#' 
#' To reduce the frequency of miss-calls, and so emphasize the difference between true
#' parent-offspring pairs and unrelated pairs, the data can be filtered on read depth.
#' Typically minimum read depth is set to 5x, but you can examine the distribution
#' of read depths with the function \code{\link{gl.report.rdepth}} and push this up 
#' with an acceptable loss of loci. 12x might be a good minimum for this particular
#' analysis. It is sensible also to push the minimum reproducibility up to 1, if 
#' that does not result in an unacceptable loss of loci. Reproducibility is stored 
#' in the slot \code{@other$loc.metrics$RepAvg} and is defined as the proportion 
#' of technical replicate assay pairs for which the marker score is consistent. 
#' You can examine the distribution of reproducibility with the function 
#' \code{\link{gl.report.reproducibility}}.
#' 
#' Note that the null expectation is not well defined, and the power reduced, if the
#' population from which the putative parent-offspring pairs are drawn contains 
#' many sibs. Note also that if an individual has been genotyped twice in the dataset, 
#' the replicate pair will be assessed by this script as being in a parent-offspring 
#' relationship.
#' 
#' You should run \code{\link{gl.report.parent.offspring}} before filtering. Use this report to 
#' decide min.rdepth and min.reproducibility and assess impact on your dataset.
#' 
#' Note that if your dataset does not contain RepAvg or rdepth among the locus metrics,
#' the filters for reproducibility and read depth are no used. 
#' 
#'\strong{ Function's output }
#'
#'  Plots and table are saved to the temporal directory (tempdir) and can be accessed with the function \code{\link{gl.print.reports}} and listed with the function \code{\link{gl.list.reports}}. Note that they can be accessed only in the current R session because tempdir is cleared each time that the R session is closed.
#'   
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return the filtered genlight object without A set of individuals in parent-offspring relationship. NULL if no parent-offspring relationships were found. 
#'
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' out <- gl.filter.parent.offspring(testset.gl[1:10])
#'
#' @seealso \code{\link{gl.list.reports}}, \code{\link{gl.report.rdepth}} ,
#'  \code{\link{gl.print.reports}},\code{\link{gl.report.reproducibility}},
#'  \code{\link{gl.report.parent.offspring}}
#'  
#' @family filter functions
#'
#' @importFrom stats median IQR
#' 
#' @import patchwork
#' 
#' @export
#'

gl.filter.parent.offspring <- function(x,
                                       min.rdepth=12,
                                       min.reproducibility=1,
                                       range=1.5,
                                       rm.monomorphs=FALSE,
                                       plot_theme = theme_dartR(), 
                                       plot_colours = two_colors, 
                                       verbose=NULL) {
  # TRAP COMMAND
  
  funname <- match.call()[[1]]
  
  # SET VERBOSITY
  
  verbose <- gl.check.verbosity(verbose)
  
  # CHECKS DATATYPE 
  
  datatype <- utils.check.datatype(x)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # FLAG SCRIPT START
  
  if (verbose >= 1) {
    if (verbose == 5) {
      cat(report("\n\nStarting", funname, "[ Build =", 
                 build, "]\n\n"))
    } else {
      cat(report("\n\nStarting", funname, "\n\n"))
    }
  }
  
  # DO THE JOB

  outliers <- gl.report.parent.offspring(x,
                                         min.rdepth=min.rdepth,
                                         min.reproducibility=min.reproducibility,
                                         range=range,
                                         plot_theme = plot_theme, 
                                         plot_colours = plot_colours, 
                                         verbose=verbose)
# Remove the outliers

  ind_to_remove <- unique(outliers$ind1)
  if(length(ind_to_remove)>0){
  x <- gl.drop.ind(x,ind.list = ind_to_remove,verbose=verbose)
  if (rm.monomorphs==TRUE){
    x <- gl.filter.monomorphs(x,verbose=verbose)
  }
  # REPORT THE RESULTS
  if(verbose>=2){
    cat("  \nInitial number of individuals:",nInd(x),"\n")
    cat("  Pairs of individuals in a parent offspring relationship:\n\n")
    print(outliers)
    cat("    \nIndividuals removed: ")
    cat(ind_to_remove, sep='\n')
    cat("\n")
  }
    }
  
  #case no 
  if (length(outliers)==0) {
    if(verbose>0){
      cat(important("No individuals were found to be in parent offspring relationship, therefore the genlight object is returned unchanged.\n"))
    }
  } 
  
  # ADD ACTION TO HISTORY
  
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("\n\nCompleted:", funname, "\n\n"))
  }
  
  # RETURN

  return(x)

}