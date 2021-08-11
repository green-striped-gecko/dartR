#' @name gl.report.callrate
#' @title Report summary of Call Rate for loci or individuals
#' @description 
#' SNP datasets generated by DArT have missing values primarily arising from 
#' failure to call a SNP because of a mutation at one or both of the restriction
#' enzyme recognition sites. P/A datasets (SilicoDArT) have missing values 
#' because it was not possible to call whether a sequence tag was amplified or 
#' not. This function tabulates the number of missing values as quantiles. 
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param method Specify the type of report by locus (method='loc') or individual
#'  (method='ind') [default 'loc'].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot_colours Vector with two colour names for the borders and fill 
#' [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session 
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, 
#' progress log ; 3, progress and results summary; 5, full report [default NULL, 
#' unless specified using gl.set.verbosity].
#'
#' @details 
#' This function expects a genlight object, containing either SNP data or 
#' SilicoDArT (=presence/absence data).
#' 
#' Callrate is summarized by locus or by individual to allow sensible decisions 
#' on thresholds for filtering taking into consideration consequential loss of 
#' data. The summary is in the form of a tabulation and plots.
#' 
#' Plot themes can be obtained from \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' Resultant ggplots and the tabulation are saved to the session's temporary 
#' directory.
#' 
#' @return Returns unaltered genlight object
#' @author Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SNP data
#'   gl.report.callrate(testset.gl)
#'   gl.report.callrate(testset.gl,method="ind")
#' # Tag P/A data
#'   gl.report.callrate(testset.gs)
#'   gl.report.callrate(testset.gs,method="ind")
#'
#' @seealso \code{\link{gl.filter.callrate}}
#' @family filters and filter reports
#' @import patchwork
#' @export
#'  

gl.report.callrate <- function(x, 
                               method = "loc", 
                               plot.out = TRUE,
                               plot_theme = theme_dartR(), 
                               plot_colours = two_colors, 
                               save2tmp = FALSE,
                               verbose = NULL) {
  
# SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
   
# FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
 
# CHECK DATATYPE 
  datatype <- utils.check.datatype(x, verbose=verbose)

# FUNCTION SPECIFIC ERROR CHECKING
    
  # Check that call rate is up to date and recalculate if necessary
    
    if (!x@other$loc.metrics.flags$CallRate) {
      x <- utils.recalc.callrate(x, verbose = 0)
    }
    
# DO THE JOB
  
########### FOR METHOD BASED ON LOCUS
  if(method=="loc"){
  if(plot.out){  
  # get title for plots
  if (method == "loc") {
    if (datatype=="SNP") {
      title1 <- "SNP data - Call Rate by Locus"
    } else {
      title1 <- "Fragment P/A data - Call Rate by Locus"
    }
    
    callrate <- x@other$loc.metrics$CallRate

    # Calculate minimum and maximum graph cutoffs for callrate
      min <- min(callrate,na.rm=TRUE)
      min <- trunc(min*100)/100
    
    # Boxplot
    p1 <- ggplot(data.frame(callrate), aes(y = callrate)) + 
      geom_boxplot(color = plot_colours[1], fill = plot_colours[2]) + 
      coord_flip() + 
      plot_theme + 
      xlim(range = c(-1, 1)) + 
      ylim(min,1) + 
      ylab(" ") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      ggtitle(title1)
    
    # Histogram
    p2 <- ggplot(data.frame(callrate), aes(x = callrate)) + 
      geom_histogram(bins = 100, color = plot_colours[1],fill = plot_colours[2]) + 
      coord_cartesian(xlim = c(min,1)) + 
      xlab("Call rate") + 
      ylab("Count") + 
      plot_theme
  }
    
    # Print out some statistics
    stats <- summary(x@other$loc.metrics$CallRate)
    cat("  Reporting Call Rate by Locus\n")
    cat("  No. of loci =", nLoc(x), "\n")
    cat("  No. of individuals =", nInd(x), "\n")
    cat("    Minimum      : ", stats[1], "\n")
    cat("    1st quartile : ", stats[2], "\n")
    cat("    Median       : ", stats[3], "\n")
    cat("    Mean         : ", stats[4], "\n")
    cat("    3r quartile  : ", stats[5], "\n")
    cat("    Maximum      : ", stats[6], "\n")
    cat("    Missing Rate Overall: ", round(sum(is.na(as.matrix(x)))/(nLoc(x) * 
                                                                        nInd(x)), 2), "\n\n")
    
    # Determine the loss of loci for a given threshold
    quantile_res <- quantile(callrate, probs = seq(0, 1, 1/20))
    retained <- unlist(lapply(quantile_res, function(y) {
      res <- length(callrate[callrate >= y])
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
  }
  
########### FOR METHOD BASED ON INDIVIDUAL
  # Calculate the call rate by individual
  if(method=="ind"){
  ind.call.rate <- 1 - rowSums(is.na(as.matrix(x)))/nLoc(x)    
  if(plot.out){
  # get title for plots
  if (datatype=="SNP") {
      title1 <- "SNP data - Call Rate by Individual"
    } else {
      title1 <- "Fragment P/A data - Call Rate by Individual"
    }

    # Calculate minimum and maximum graph cutoffs for callrate
    min <- min(ind.call.rate)
    min <- trunc(min*100)/100 
    
    # Boxplot
    p1 <- ggplot(data.frame(ind.call.rate), aes(y = ind.call.rate)) + 
      geom_boxplot(color = plot_colours[1], fill = plot_colours[2]) + 
      coord_flip() + 
      plot_theme + 
      xlim(range = c(-1, 1)) + 
      ylim(min,1) +
      ylab(" ") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      ggtitle(title1)
    
    # Histogram
    p2 <- ggplot(data.frame(ind.call.rate), aes(x = ind.call.rate)) + 
      geom_histogram(bins = 100, color = plot_colours[1], fill = plot_colours[2]) +
      coord_cartesian(xlim = c(min,1)) + 
      xlab("Call rate") + 
      ylab("Count") + 
      plot_theme
  }
    # Print out some statistics
    stats <- summary(ind.call.rate)
    cat("  Reporting Call Rate by Individual\n")
    cat("  No. of loci =", nLoc(x), "\n")
    cat("  No. of individuals =", nInd(x), "\n")
    cat("    Minimum      : ", stats[1], "\n")
    cat("    1st quartile : ", stats[2], "\n")
    cat("    Median       : ", stats[3], "\n")
    cat("    Mean         : ", stats[4], "\n")
    cat("    3r quartile  : ", stats[5], "\n")
    cat("    Maximum      : ", stats[6], "\n")
    cat("    Missing Rate Overall: ", round(sum(is.na(as.matrix(x)))/(nLoc(x) * 
                                                                        nInd(x)), 2), "\n\n")
    
    # Determine the loss of individuals for a given
    # threshold using quantiles
    quantile_res <- quantile(ind.call.rate, 
                             probs = seq(0, 1, 1/20))
    retained <- unlist(lapply(quantile_res, function(y) {
      res <- length(ind.call.rate[ind.call.rate >= y])
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
  if(plot.out){
    # using package patchwork
    p3 <- (p1/p2) + plot_layout(heights = c(1, 4))
    print(p3)
  }
  print(df)
  
# SAVE INTERMEDIATES TO TEMPDIR             

  # creating temp file names
  if(save2tmp){
    if(plot.out){
  temp_plot <- tempfile(pattern = "Plot_")
  match_call <- paste0(names(match.call()),"_",as.character(match.call()),collapse = "_")
  # saving to tempdir
  saveRDS(list(match_call,p3), file = temp_plot)
  if(verbose>=2){
    cat(report("  Saving the ggplot to session tempfile\n"))
  }
    }
    temp_table <- tempfile(pattern = "Table_")
    saveRDS(list(match_call,df), file = temp_table)
  if(verbose>=2){
    cat(report("  Saving tabulation to session tempfile\n"))
    cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))
  }
  }
  
# FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

# RETURN
  invisible(x)

}
