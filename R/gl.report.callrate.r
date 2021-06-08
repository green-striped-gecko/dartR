#' @name gl.report.callrate
#'
#' @title Report summary of Call Rate for loci or individuals
#'
#' @description 
#' SNP datasets generated by DArT have missing values primarily
#' arising from failure to call a SNP because of a mutation at one or both of
#' the the restriction enzyme recognition sites. This function reports the
#' number of missing values for each of several quantiles. 
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param method Specify the type of report by locus (method='loc') or individual
#'  (method='ind') [default method='loc'].
#' @param plot_theme Theme for the plot. See Details for options [default theme_dartR()].
#' @param plot_colours List of two color names for the borders and fill of the
#'  plots [default two_colors].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity].
#'
#' @details The function \code{\link{gl.filter.callrate}} will filter out the
#'  loci with call rates below a specified threshold.
#'
#'  Tag Presence/Absence datasets (SilicoDArT) have missing values where it is
#'  not possible to determine reliably if the sequence tag can be called at a
#'  particular locus.
#'  
#'  Quantiles are
#' partitions of a finite set of values into q subsets of (nearly) equal sizes.
#' In this function q = 20. Quantiles are useful measures because they are less
#'  susceptible to long-tailed distributions and outliers.
#'  
#'\strong{ Function's output }
#'
#'  The minimum, maximum, mean and a tabulation of call rate quantiles against
#'  thresholds are provided. Output also includes a boxplot and a
#'  histogram to guide in the selection of a threshold for filtering on callrate.
#'
#'  Plots and table are saved to the temporal directory (tempdir) and can be accessed with the function \code{\link{gl.print.reports}} and listed with the function \code{\link{gl.list.reports}}. Note that they can be accessed only in the current R session because tempdir is cleared each time that the R session is closed.
#'   
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return Returns unaltered genlight object
#'
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' # SNP data
#'   gl.report.callrate(testset.gl)
#'   gl.report.callrate(testset.gl,method="ind")
#' # Tag P/A data
#'   gl.report.callrate(testset.gs)
#'   gl.report.callrate(testset.gs,method="ind")
#'
#' @seealso \code{\link{gl.filter.callrate}}, \code{\link{gl.list.reports}},
#'  \code{\link{gl.print.reports}}
#'  
#' @family filters and filter reports
#'
#' @import patchwork
#'
#' @export
#'  

gl.report.callrate <- function(x, 
                               method = "loc", 
                               plot_theme = theme_dartR(), 
                               plot_colours = two_colors, 
                               verbose = options()$dartR_verbose) {
  
# TRAP COMMAND
  
  funname <- match.call()[[1]]
  
# GENERAL ERROR CHECKING, SETTING VERBOSITY AND DATATYPE 
  
  datatype <- NULL
  utils.check.gl(x,env=environment())

# FUNCTION SPECIFIC ERROR CHECKING
    
  # Check that call rate is up to date and recalculate if necessary
    
    if (!x@other$loc.metrics.flags$monomorphs) {
      x <- utils.recalc.callrate(x, verbose = 0)
    }
    
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
  if (method == "loc") {
    if (datatype=="SNP") {
      title1 <- "SNP data - Call Rate by Locus"
    } else {
      title1 <- "Fragment P/A data - Call Rate by Locus"
    }
    
    callrate <- data.frame(x@other$loc.metrics$CallRate)
    colnames(callrate) <- "callrate"
    
    # Boxplot
    p1 <- ggplot(callrate, aes(y = callrate)) + 
      geom_boxplot(color = plot_colours[1], fill = plot_colours[2]) + 
      coord_flip() + 
      plot_theme + 
      xlim(range = c(-1, 1)) + 
      ylim(0,1) + ylab(" ") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      ggtitle(title1)
    
    # Histogram
        p2 <- ggplot(callrate, aes(x = callrate)) + 
      geom_histogram(bins = 50, color = plot_colours[1],fill = plot_colours[2]) + 
      coord_cartesian(xlim = c(0,1)) + 
      xlab("Call rate") + 
      ylab("Count") + 
      plot_theme
    
    # Print out some statistics
    stats <- summary(x@other$loc.metrics$CallRate)
    cat("  Reporting Call Rate by Locus\n")
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
    quantile_res <- quantile(callrate$callrate, 
                             probs = seq(0, 1, 1/20))
    retained <- unlist(lapply(quantile_res, function(y) {
      res <- length(callrate$callrate[callrate$callrate >= 
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
  }
  
########### FOR METHOD BASED ON INDIVIDUAL
  
  # get title for plots
  if (method == "ind") {
    if (all(x@ploidy == 2)) {
      title1 <- "SNP data - Call Rate by Individual"
    } else {
      title1 <- "Fragment P/A data - Call Rate by Individual"
    }
    
    # Calculate the call rate by individual
    ind.call.rate <- data.frame(1 - rowSums(is.na(as.matrix(x)))/nLoc(x))
    colnames(ind.call.rate) <- "ind.call.rate"
    
    # Boxplot
    p1 <- ggplot(ind.call.rate, aes(y = ind.call.rate)) + 
      geom_boxplot(color = plot_colours[1], fill = plot_colours[2]) + 
      coord_flip() + 
      plot_theme + 
      xlim(range = c(-1, 1)) + 
      ylim(0,1) +
      ylab(" ") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      ggtitle(title1)
    
    # Histogram
    p2 <- ggplot(ind.call.rate, aes(x = ind.call.rate)) + 
      geom_histogram(bins = 50, color = plot_colours[1], fill = plot_colours[2]) +
      coord_cartesian(xlim = c(0,1)) + 
      xlab("Call rate") + 
      ylab("Count") + 
      plot_theme
    
    # Print out some statistics
    stats <- summary(ind.call.rate$ind.call.rate)
    cat("  Reporting Call Rate by Individual\n")
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
    
    # Determine the loss of individuals for a given
    # threshold using quantiles
    quantile_res <- quantile(ind.call.rate$ind.call.rate, 
                             probs = seq(0, 1, 1/20))
    retained <- unlist(lapply(quantile_res, function(y) {
      res <- length(ind.call.rate$ind.call.rate[ind.call.rate$ind.call.rate >= 
                                                  y])
    }))
    pc.retained <- round(retained * 100/nInd(x), 
                         1)
    filtered <- nInd(x) - retained
    pc.filtered <- 100 - pc.retained
    df <- data.frame(as.numeric(sub("%", "", names(quantile_res))), 
                     quantile_res, retained, pc.retained, filtered, 
                     pc.filtered)
    colnames(df) <- c("Quantile", "Threshold", 
                      "Retained", "Percent", "Filtered", "Percent")
    df <- df[order(-df$Quantile), ]
    df$Quantile <- paste0(df$Quantile, "%")
    rownames(df) <- NULL
  }
  
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
