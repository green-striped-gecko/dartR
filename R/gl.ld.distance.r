#' @name gl.ld.distance
#' @title Plots linkage disequilibrium against distance by population 
#' disequilibrium patterns
#' @description
#' The function creates a plot showing
#' the pairwise LD measure against distance in number of base pairs pooled over
#' all the chromosomes and a red line representing the threshold (R.squared = 
#' 0.2) that is commonly used to imply that two loci are unlinked (Delourme et
#' al., 2013; Li et al., 2014).
#' @param ld_report Output from function \code{\link{gl.report.ld.map}} 
#' [required].
#' @param ld_resolution Resolution at which LD should be reported in number of 
#' base pairs [default NULL].
#' @param pop_colors A color palette for box plots by population or a list
#'  with as many colors as there are populations in the dataset
#' [default NULL].
#' @param plot_theme User specified theme [default NULL].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @references
#' \itemize{
#' \item Delourme, R., Falentin, C., Fomeju, B. F., Boillot, M., Lassalle, G., 
#' André, I., . . . Marty, A. (2013). High-density SNP-based genetic map 
#' development and linkage disequilibrium assessment in Brassica napusL. BMC 
#' genomics, 14(1), 120.
#' \item Li, X., Han, Y., Wei, Y., Acharya, A., Farmer, A. D., Ho, J., . . . 
#' Brummer, E. C. (2014). Development of an alfalfa SNP array and its use to 
#' evaluate patterns of population structure and linkage disequilibrium. PLoS 
#' One, 9(1), e84329.
#'  }
#' @return A dataframe with information of LD against distance by population.  
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' if ((requireNamespace("snpStats", quietly = TRUE)) & (requireNamespace("fields", quietly = TRUE))) {
#' require("dartR.data")
#' x <- platypus.gl
#' x <- gl.filter.callrate(x,threshold = 1)
#' x <- gl.filter.monomorphs(x)
#' x$position <- x$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#' x$chromosome <- as.factor(x$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
#' ld_res <- gl.report.ld.map(x,ld_max_pairwise = 10000000)
#' ld_res_2 <- gl.ld.distance(ld_res,ld_resolution= 1000000)
#' }
#' @family ld functions
#' @export

gl.ld.distance <- function(ld_report,
                           ld_resolution = 100000,
                           pop_colors = NULL,
                           plot_theme = NULL,
                           plot.out = TRUE,
                           save2tmp = FALSE,
                           verbose = NULL){
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbose = verbose)

  # DO THE JOB
  
    ld_max_pairwise <- max(ld_report$distance)
    break_bins <- c(seq(1, ld_max_pairwise, ld_resolution), ld_max_pairwise)
    
    split_df <- split(ld_report, f = ld_report$pop)
    split_df <- lapply(split_df, function(x) {
      x[order(x$distance), ]
    })
    bins_ld_temp <-
      lapply(split_df, function(x) {
        fields::stats.bin(x$distance, x$ld_stat, breaks = break_bins)
      })
    bins_ld <-
      lapply(seq_along(bins_ld_temp), function(i) {
        as.data.frame(cbind(
          names(bins_ld_temp[i]),
          unname(bins_ld_temp[[i]]$breaks[2:length(bins_ld_temp[[i]]$breaks)]),
          unname(bins_ld_temp[[i]]$stats[2, ])
        ))
      })
    bins_ld <- rbindlist(bins_ld)
    colnames(bins_ld) <- c("pop", "distance", "ld_stat")
    bins_ld$pop <- as.factor(bins_ld$pop)
    bins_ld$distance <- as.numeric(bins_ld$distance)
    bins_ld$ld_stat <- as.numeric(bins_ld$ld_stat)
  
  # pairwise LD by population
    
    pops <- as.factor(unique(ld_report$pop))
      
    if(is.null(plot_theme)){
      plot_theme <- theme_dartR()
    }
    
    if(is.null(pop_colors)){
      pop_colors <- discrete_palette(length(levels(pops)))
    }
    
    if (is(pop_colors, "function")) {
      pop_colors <- pop_colors(length(levels(pops)))
    }
    
    if (!is(pop_colors,"function")) {
      pop_colors <- pop_colors
    }
    
  distance <- ld_stat <- NULL
  p3 <-
    ggplot(bins_ld, aes(x = distance, y = ld_stat, colour = pop)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_hline(aes(yintercept = 0.2, 
                   colour = "LD threshold for unlinked loci"),color="red",
               size = 1) +
    xlab("Base pairs") +
    ylab("Linkage disequilibrium") +
    labs(color = "") +
    scale_color_manual(values = pop_colors) +
    plot_theme +
    theme(legend.position = "bottom")
  
  
  # PRINTING OUTPUTS
  if(plot.out){
    print(p3)
  }
  print(bins_ld,row.names = FALSE)
  
  # SAVE INTERMEDIATES TO TEMPDIR creating temp file names
  if (save2tmp) {
    if (plot.out) {
      temp_plot <- tempfile(pattern = "Plot_")
      match_call <-
        paste0(names(match.call()),
               "_",
               as.character(match.call()),
               collapse = "_")
      # saving to tempdir
      saveRDS(list(match_call, p3), file = temp_plot)
      if (verbose >= 2) {
        cat(report("  Saving the ggplot to session tempfile\n"))
      }
    }
    temp_table <- tempfile(pattern = "Table_")
    saveRDS(list(match_call, bins_ld), file = temp_table)
    if (verbose >= 2) {
      cat(report("  Saving tabulation to session tempfile\n"))
      cat(
        report(
          "  NOTE: Retrieve output files from tempdir using
                    gl.list.reports() and gl.print.reports()\n"
        )
      )
    }
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  
  return(invisible(bins_ld))
  
}







