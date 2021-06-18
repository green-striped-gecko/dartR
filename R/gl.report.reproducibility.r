#' @name gl.report.reproducibility
#'
#' @title Report summary of RepAvg (repeatability averaged over both alleles for 
#' each locus) or reproducibility (repeatability of the scores for fragment presence/absence)
#'
#' @description 
#' SNP datasets generated by DArT have an index, RepAvg, generated by reproducing 
#' the data independently for 30% of loci. RepAvg is the proportion of alleles that 
#' give a repeatable result, averaged over both alleles for each locus.
#' 
#' In the case of fragment presence/absence data (SilicoDArT), repeatability is
#'  the percentage of scores that are repeated in the technical replicate dataset.
#' 
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param plot_theme Theme for the plot. See Details for options [default theme_dartR()].
#' @param plot_colours List of two color names for the borders and fill of the
#'  plots [default two_colors].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity]
#'
#' @details The function \code{\link{gl.filter.reproducibility}} will filter out the
#'  loci with repeatability below a specified threshold.
#'  
#'  Quantiles are
#' partitions of a finite set of values into q subsets of (nearly) equal sizes.
#' In this function q = 20. Quantiles are useful measures because they are less
#'  susceptible to long-tailed distributions and outliers.
#'  
#'\strong{ Function's output }
#'
#'  The minimum, maximum, mean and a tabulation of repeatability quantiles against
#'  thresholds are provided. Output also includes a boxplot and a
#'  histogram to guide in the selection of a threshold for filtering on repeatability.
#'
#'  Plots and table are saved to the temporal directory (tempdir) and can be accessed with the function \code{\link{gl.print.reports}} and listed with the function \code{\link{gl.list.reports}}. Note that they can be accessed only in the current R session because tempdir is cleared each time that the R session is closed.
#'   
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return An unaltered genlight object
#'
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' # SNP data
#'   out <- gl.report.reproducibility(testset.gl)
#' # Tag P/A data
#'   out <- gl.report.reproducibility(testset.gs)
#'
#' @seealso \code{\link{gl.filter.reproducibility}}, \code{\link{gl.list.reports}},
#'  \code{\link{gl.print.reports}}
#'  
#' @family filters and filter reports
#'
#' @import patchwork
#'
#' @export
#'  

gl.report.reproducibility <- function(x,  
                                      plot_theme = theme_dartR(),  
                                      plot_colours = two_colors, 
                                      verbose = NULL) {
  
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
  
  ########### FOR METHOD BASED ON LOCUS
  
  # get title for plots
  if (all(x@ploidy==2)){
    title <- paste0("SNP data (DArTSeq)\nRepeatability by Locus")
  } else {
    title <- paste0("Fragment P/A data (SilicoDArT)\nRepeatability by Locus")
  } 
  
  if (all(x@ploidy == 2)){
    repeatability <- x@other$loc.metrics$RepAvg
  } else {
    repeatability <- x@other$loc.metrics$Reproducibility
  } 
  
  repeatability_plot <- data.frame(repeatability)
  colnames(repeatability_plot) <- "repeatability"
  
  # Boxplot
  p1 <- ggplot(repeatability_plot, aes(y = repeatability)) + 
    geom_boxplot(color = plot_colours[1], fill = plot_colours[2]) + 
    coord_flip() + 
    plot_theme + 
    ylim(c(min(repeatability), 1)) + 
    ylab(" ") + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    ggtitle(title)
  
  # Histogram
  p2 <- ggplot(repeatability_plot, aes(x = repeatability)) + 
    geom_histogram(bins = 50, color = plot_colours[1],fill = plot_colours[2]) + 
    coord_cartesian(xlim = c(min(repeatability),1)) + 
    xlab("Repeatability") + 
    ylab("Count") + 
    plot_theme
  
  # Print out some statistics
  stats <- summary(repeatability)
  cat("  Reporting Repeatability by Locus\n")
  cat("  No. of loci =", nLoc(x), "\n")
  cat("  No. of individuals =", nInd(x), "\n")
  cat("    Minimum      : ", stats[1], "\n")
  cat("    1st quantile : ", stats[2], "\n")
  cat("    Median       : ", stats[3], "\n")
  cat("    Mean         : ", stats[4], "\n")
  cat("    3r quantile  : ", stats[5], "\n")
  cat("    Maximum      : ", stats[6], "\n")
  cat("    Missing Rate Overall: ", round(sum(is.na(as.matrix(x)))/(nLoc(x) * 
                                                                      nInd(x)), 2), "\n\n")
  
  # Determine the loss of loci for a given threshold
  # using quantiles
  quantile_res <- quantile(repeatability, 
                           probs = seq(0, 1, 1/20))
  retained <- unlist(lapply(quantile_res, function(y) {
    res <- length(repeatability[repeatability >= 
                                  y])
  }))
  pc.retained <- round(retained * 100/nLoc(x), 
                       1)
  filtered <- nLoc(x) - retained
  pc.filtered <- 100 - pc.retained
  df <- data.frame(as.numeric(sub("%", "", names(quantile_res))), 
                   quantile_res, retained, pc.retained, filtered, 
                   pc.filtered)
  colnames(df) <- c("Quantile", "Threshold", 
                    "Retained", "Percent", "Filtered", "Percent")
  df <- df[order(-df$Quantile), ]
  df$Quantile <- paste0(df$Quantile, "%")
  rownames(df) <- NULL
  
  # PRINTING OUTPUTS
  # using package patchwork
  p3 <- (p1/p2) + plot_layout(heights = c(1, 4))
  print(p3)
  print(df)
  
  # SAVE INTERMEDIATES TO TEMPDIR             
  # creating temp file names
  temp_plot <- tempfile(pattern =paste0("dartR_plot",paste0(names(match.call()),"_",as.character(match.call()),collapse = "_"),"_"))
  temp_table <- tempfile(pattern = paste0("dartR_table",paste0(names(match.call()),"_",as.character(match.call()),collapse = "_"),"_"))
  
  # saving to tempdir
  saveRDS(p3, file = temp_plot)
  if(verbose>=2){
    cat(report("  Saving the plot in ggplot format to the tempfile as",temp_plot,"using saveRDS\n"))
  }
  saveRDS(df, file = temp_table)
  if(verbose>=2){
    cat(report("  Saving the report to the tempfile as",temp_table,"using saveRDS\n"))
  }
  if(verbose>=2){
    cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("\n\nCompleted:", funname, "\n\n"))
  }
  
  # RETURN
  
  invisible(x)
  
}
