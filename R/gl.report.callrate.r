#' Report summary of Call Rate for loci or individuals
#'
#' SNP datasets generated by DArT have missing values primarily arising from failure to call a SNP because of a mutation
#' at one or both of the the restriction enzyme recognition sites. This function reports the number of missing values for each
#' of several quantiles. Quantiles are partitions of a finite set of values into q subsets of (nearly) equal sizes. In this function q = 20. Quantiles are useful measures because they are less susceptible to long-tailed distributions and outliers.
#' 
#' The function \code{\link{gl.filter.callrate}} will filter out the loci with call rates below a specified threshold.
#'
#' Tag Presence/Absence datasets (SilicoDArT) have missing values where it is not possible to determine reliably if the
#' sequence tag can be called at a particular locus.
#'
#' The minimum, maximum, mean and a tabulation of call rate quantiles against thresholds rate are provided. Output also includes a boxplot and a histogram.
#' 
#' \strong{Plots and table are saved to the temporal directory (tempdir) and can be accessed with the function \code{\link{gl.access.report}}. Note that they can be accessed only in the current R session because tempdir is cleared each time that an R session is closed.}
#' 
#' #' Examples of other themes that can be used can be consulted in 
#' \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'
#' @param x -- name of the genlight object containing the SNP or presence/absence (SilicoDArT) data [required]
#' @param method -- specify the type of report by locus (method="loc") or individual (method="ind") [default method="loc"]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @param plot_theme -- theme for the plot. See Details for options [default theme_dartR()]
#' @param plot_colours -- two colour names for borders and fill of the plots [default wes_palette("Zissou1")]
#' @return returns a genlight object with the file names of plots and table that were saved in the tempdir stored in the slot other$history
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' # SNP data
#'   out <- gl.report.callrate(testset.gl)
#' # Tag P/A data
#'   out <- gl.report.callrate(testset.gs)

gl.report.callrate <-
  function(x,
           method = "loc",
           plot_theme = theme_dartR(),
           plot_colours = wes_palette("Zissou1"),
           verbose = NULL) {
    
    # TRAP COMMAND, SET VERSION
    funname <- match.call()[[1]]
    build <- "Jacob"
    
    # ERROR CHECKING
    x <- utils.check.gl(x,verbose)
    verbose <- x@other$verbose
    
    # FLAG SCRIPT START
    if (verbose >= 1) {
      if (verbose == 5) {
        cat(report("Starting", funname, "[ Build =", build, "]\n\n"))
      } else {
        cat(report("Starting", funname, "\n\n"))
      }
    }

    # DO THE JOB
    
    # RECALCULATE CALL RATE, IF NOT PREVIOUSLY DONE
    # if (!x@other$loc.metrics.flags$monomorphs){
      x <- dartR:::utils.recalc.callrate(x, verbose=0)
     # }

    ########### FOR METHOD BASED ON LOCUS
      # get title for plots    
      if (method == "loc") {
      if (all(x@ploidy == 2)) {
        title1 <- "SNP data - Call Rate by Locus"
      } else {
        title1 <- "Fragment P/A data - Call Rate by Locus"
      }
      
      callrate <- data.frame(x@other$loc.metrics$CallRate)
      colnames(callrate) <- "callrate"
      
      p1 <- ggplot(callrate, aes(y = callrate)) +
        geom_boxplot(color=plot_colours[1],fill=plot_colours[2]) +
        coord_flip() +
        plot_theme +
        xlim(range = c(-1, 1)) +
        ylim(0, 1) + 
        ylab(" ") +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank())+
        ggtitle(title1)
      
      p2 <- ggplot(callrate, aes(x = callrate)) +
        geom_histogram(bins = 50,color=plot_colours[1],fill=plot_colours[2]) +
        coord_cartesian(xlim = c(0, 1)) +
        xlab("Call rate") +
        ylab("Count") +
        plot_theme
      
      # Print out some statistics
      cat("  Reporting Call Rate by Locus\n")
      cat("  No. of loci =", nLoc(x), "\n")
      cat("  No. of individuals =", nInd(x), "\n")
      cat("    Minimum Call Rate: ", round(
        min(x@other$loc.metrics$CallRate), 2
      ), "\n")
      cat("    Maximum Call Rate: ", round(
        max(x@other$loc.metrics$CallRate), 2
      ), "\n")
      cat("    Average Call Rate: ", round(
        mean(x@other$loc.metrics$CallRate), 3
      ), "\n")
      cat("    Missing Rate Overall: ", round(sum(
        is.na(as.matrix(x))
      ) / (
        nLoc(x) * nInd(x)
      ), 2), "\n\n")
      
      # Determine the loss of loci for a given threshold using quantiles
      quantile_res <- quantile(callrate$callrate,probs = seq(0,1,1/20))
      retained <- unlist(lapply(quantile_res, function(y){res <- length(callrate$callrate[callrate$callrate>=y]) }))
      pc.retained <- round(retained*100/nLoc(x),1)
      filtered <- nLoc(x) - retained
      pc.filtered <- 100 - pc.retained
      df <- data.frame(as.numeric(sub("%","",names(quantile_res))),quantile_res,retained,pc.retained,filtered,pc.filtered)
      colnames(df) <- c("Quantile","Threshold", "Retained", "Percent", "Filtered", "Percent")
      df <- df[order(-df$Quantile),]
      df$Quantile <- paste0(df$Quantile,"%")
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
      ind.call.rate <- data.frame(1 - rowSums(is.na(as.matrix(x))) / nLoc(x))
      colnames(ind.call.rate) <- "ind.call.rate"

      # Boxplot
      p1 <- ggplot(ind.call.rate, aes(y = ind.call.rate)) +
        geom_boxplot(color=plot_colours[1],fill=plot_colours[2]) +
        coord_flip() +
        plot_theme +
        xlim(range = c(-1, 1)) +
        ylim(0, 1) + 
        ylab(" ") +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank())+
        ggtitle(title1)
      
      # Histogram
      p2 <- ggplot(data.frame(ind.call.rate), aes(x = ind.call.rate)) +
        geom_histogram(bins = 50,color=plot_colours[1],fill=plot_colours[2]) +
        coord_cartesian(xlim = c(0, 1)) +
        xlab("Call rate") +
        ylab("Count") +
        plot_theme
      
      # Print out some statistics
      cat("  Reporting Call Rate by Individual\n")
      cat("  No. of loci =", nLoc(x), "\n")
      cat("  No. of individuals =", nInd(x), "\n")
      cat("    Minimum Call Rate: ", round(min(ind.call.rate$ind.call.rate), 2), "\n")
      cat("    Maximum Call Rate: ", round(max(ind.call.rate$ind.call.rate), 2), "\n")
      cat("    Average Call Rate: ", round(mean(ind.call.rate$ind.call.rate), 3), "\n")
      cat("    Missing Rate Overall: ", round(sum(
        is.na(as.matrix(x))
      ) / (
        nLoc(x) * nInd(x)
      ), 2), "\n\n")
      
      # Determine the loss of individuals for a given threshold using quantiles
      quantile_res <- quantile(ind.call.rate$ind.call.rate,probs = seq(0,1,1/20))
      retained <- unlist(lapply(quantile_res, function(y){res <- length(ind.call.rate$ind.call.rate[ind.call.rate$ind.call.rate>=y]) }))
      pc.retained <- round(retained*100/nInd(x),1)
      filtered <- nInd(x) - retained
      pc.filtered <- 100 - pc.retained
      df <- data.frame(as.numeric(sub("%","",names(quantile_res))),quantile_res,retained,pc.retained,filtered,pc.filtered)
      colnames(df) <- c("Quantile","Threshold", "Retained", "Percent", "Filtered", "Percent")
      df <- df[order(-df$Quantile),]
      df$Quantile <- paste0(df$Quantile,"%")
      rownames(df) <- NULL
    }
    
    # printing outputs
    p3 <- p1 / p2
    print(p3)
    print(df)
    # creating file names
    temp_plot <- tempfile(pattern = "plot_")
    temp_table <- tempfile(pattern = "table_")
    # saving to tempdir
    saveRDS(p3,file=temp_plot)
    saveRDS(df,file=temp_table)
  
    # ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- c(match.call(),temp_plot,temp_table)
    
    # FLAG SCRIPT END
    if (verbose >= 1) {
      cat(report("\n\nCompleted:", funname, "\n\n"))
    }
    
    cat(important("Plots and table were saved to the temporal directory (tempdir) and can be\n accesed with the function gl.access.report(). Note that they can be accessed\n only in the current R session because tempdir is cleared each time that the\n R session is closed.\n\n"))
  
     return(x)
  }
