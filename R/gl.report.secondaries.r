#' @name gl.report.secondaries
#'
#' @title Report loci containing secondary SNPs in sequence tags 
#'
#' @description 
#' SNP datasets generated by DArT include fragments with more than one SNP (that
#'  is, with secondaries) and record them separately with the same CloneID (=AlleleID).
#' These multiple SNP loci within a fragment are likely to be linked, and so you may 
#' wish to remove secondaries.
#' 
#' This function reports statistics associated with secondaries, and the consequences 
#' of filtering them out, and provides three plots. The first is a boxplot, the second is a barplot of the frequency 
#' of secondaries per sequence tag, and the third is the Poisson expectation for those
#' frequencies including an estimate of the zero class (no. of sequence tags with 
#' no SNP scored).
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param plot_theme Theme for the plot. See Details for options [default theme_dartR()].
#' @param plot_colours List of two color names for the borders and fill of the
#'  plots [default two_colors].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity]
#'
#' @details The function \code{\link{gl.filter.secondaries}} will filter out the
#'  loci with secondaries.
#'  
#' Heterozygosity as estimated by the function \code{\link{gl.report.heterozygosity}} 
#' is in a sense relative, because it is calculated
#' against a background of only those loci that are polymorphic somewhere in the dataset.
#' To allow intercomparability across studies and species, any measure of heterozygosity
#' needs to accommodate loci that are invariant. However, the number of invariant loci
#' are unknown given the SNPs are detected as single point mutational variants and 
#' invariant sequences are discarded, and because of
#' the particular additional filtering pre-analysis. Modelling the counts
#' of SNPs per sequence tag as a Poisson distribution in this script allows estimate of the zero class,
#' that is, the number of invariant loci. This is reported, and the veracity of the 
#' estimate can be assessed by the correspondence of the observed frequencies against
#' those under Poisson expectation in the associated graphs. The number of invariant loci can then be optionally
#' provided to the function \code{\link{gl.report.heterozygosity}} via the parameter n.invariants.
#'  
#'\strong{ Function's output }
#'
#'  Plots are saved to the temporal directory (tempdir) and can be accessed with the function \code{\link{gl.print.reports}} and listed with the function \code{\link{gl.list.reports}}. Note that they can be accessed only in the current R session because tempdir is cleared each time that the R session is closed.
#'   
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return A genlight object containing only those loci with secondaries
#'
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' out <- gl.report.secondaries(bandicoot.gl)
#'
#' @seealso \code{\link{gl.filter.secondaries}}, \code{\link{gl.list.reports}},
#'  \code{\link{gl.print.reports}},\code{\link{gl.report.heterozygosity}}
#'  
#' @family filters and filter reports
#' 
#' @importFrom stats dpois
#'
#' @import patchwork
#'
#' @export
#'  

gl.report.secondaries <- function(x, 
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

  # Extract the clone ID number
  a <- strsplit(as.character(x@other$loc.metrics$AlleleID),"\\|")
  b <- unlist(a)[ c(TRUE,FALSE,FALSE) ]
  if (verbose >= 2) {
    cat(report("Counting ....\n"))
  }
  x.secondaries <- x[,duplicated(b)]
   
  nloc.with.secondaries <- table(duplicated(b))[2]
  if (!is.na(nloc.with.secondaries)){
    
    freqs_1 <- c(0,as.numeric(table(b)))
    secondaries_plot <-  as.data.frame(freqs_1)
    colnames(secondaries_plot) <- "freqs"
    
    # Boxplot
    p1 <- ggplot(secondaries_plot, aes(y = freqs)) + 
      geom_boxplot(color = plot_colours[1], fill = plot_colours[2]) + 
      coord_flip() + 
      plot_theme + 
      xlim(range = c(-1, 1)) + 
      scale_y_discrete(limits=c(as.character(unique(freqs_1)))) +
      ylab("Frequency") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      ggtitle("Boxplot")
    
    # Barplot
    freqs_2 <- c(0,table(as.numeric(table(b))))
    secondaries_plot_2 <- as.table(freqs_2)
    names(secondaries_plot_2)<- seq(1:(length(secondaries_plot_2)))-1
    secondaries_plot_2 <- as.data.frame(secondaries_plot_2)
    colnames(secondaries_plot_2) <- c("freq","count")
    
    freq <- NULL
    
    p2 <- ggplot(secondaries_plot_2,aes(x=freq,y=count)) + 
      geom_col(color = plot_colours[1],fill = plot_colours[2]) + 
      xlab("Frequency") + 
      ylab("Count") + 
      ggtitle("Observed Frequency of SNPs per Sequence Tag") +
      plot_theme
    
    # Plot Histogram with estimate of the zero class
    if (verbose >= 2){
      cat(report("Estimating parameters (lambda) of the Poisson expectation\n"))
    }
    
      # Calculate the mean for the truncated distribution
        freqs <- c(0,table(as.numeric(table(b))))
        tmp <- NA
        for (i in 1:length(freqs)){
          tmp[i] <- freqs[i]*(i-1)
        }
        tmean <- sum(tmp)/sum(freqs)
        
      # Set a random seed, close to 1
        seed <- tmean
        
      # Set convergence criterion
        delta <- 0.00001
        
      # Use the mean of the truncated distribution to compute lambda for the untruncated distribution
        k <- seed
        for (i in 1:100){
          if (verbose >= 2){print(k)}
          k.new <- tmean*(1-exp(-k))
          if (abs(k.new - k) <= delta){
            if (verbose >= 2){cat("Converged on Lambda of",k.new,"\n")}
            fail <- FALSE
            break
          }
          if (i == 100){
            if(verbose >= 2){
              cat(important("Failed to converge: No reliable estimate of invariant loci\n"))
                                 }
            fail <- TRUE
            break
          }
          k <- k.new
        }
        
      # Size of the truncated distribution
        if (!fail) {
          n <- sum(freqs)  # Size of the truncated set 
          tp <- 1 - dpois( x=0, lambda=k ) # Fraction that is the truncated set
          rn <- round(n/tp,0) # Estimate of the whole set
          cat("Estimated size of the zero class",round(dpois(x=0,lambda=k)*rn,0),"\n")
          # Table for the reconstructed set  
            reconstructed <- dpois( x=0:(length(freqs)-1), lambda=k )*rn
            reconstructed <- as.table(reconstructed)
            names(reconstructed)<- seq(1:(length(reconstructed)))-1
            
            title <- paste0("Poisson Expectation (zero class ",round(dpois(x=0,lambda=k)*rn,0)," invariant loci)")
            
            reconstructed_plot <- as.data.frame(reconstructed)
            colnames(reconstructed_plot) <- c("freq","count")
            
             # Barplot
              p3 <- ggplot(reconstructed_plot,aes(x=freq,y=count)) + 
              geom_col(color = plot_colours[1],fill = plot_colours[2]) + 
              xlab("Frequency") + 
              ylab("Count") + 
              ggtitle(title) +
              plot_theme
              
              # PRINTING OUTPUTS
              # using package patchwork
              p4 <- p1/p2/p3
              print(p4)
        }
        
        if (fail) {
          p4 <- p1/p2
          print(p4)
        }
        
        # SAVE INTERMEDIATES TO TEMPDIR   
        # creating temp file names
        temp_plot <- tempfile(pattern = "dartR_plot_")
        match_call <- paste0(names(match.call()),"_",as.character(match.call()),collapse = "_")

        # saving to tempdir
        saveRDS(list(match_call,p4), file = temp_plot)
        if(verbose>=2){
          cat(report("  Saving the plot in ggplot format to the tempfile as",temp_plot,"using saveRDS\n"))
        }
        if(verbose>=2){
          cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))
        }
        if(verbose >= 2){
          cat(report("\nReturning a genlight object containing only those loci with secondaries (multiple entries per locus)\n"))
        }
        
  } else {
      cat(important("  Warning: No loci with secondaries, no plot produced\n"))
  }
  
# Identify secondaries in the genlight object
  cat("  Total number of SNP loci scored:",nLoc(x),"\n")
  if (is.na(table(duplicated(b))[2])) {
    cat("    Number of secondaries: 0 \n")
  } else {
      cat("   Number of sequence tags in total:",table(duplicated(b))[1],"\n")
      if (fail){
        cat("    Number of invariant sequence tags cannot be estimated\n")
      } else {
        cat("   Estimated number of invariant sequence tags:", round(dpois(x=0,lambda=k)*rn,0),"\n")
      }  
      cat("    Number of sequence tags with secondaries:",sum(table(as.numeric(table(b))))-table(as.numeric(table(b)))[1],"\n")
      cat("    Number of secondary SNP loci that would be removed on filtering:",table(duplicated(b))[2],"\n")
      cat("    Number of SNP loci that would be retained on filtering:",table(duplicated(b))[1],"\n")
      if(verbose >= 3){
        cat(" Tabular 1 to K secondaries (refer plot)\n",table(as.numeric(table(b))),"\n")
        }
  }  
  
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("\n\nCompleted:", funname, "\n\n"))
  }
  
  # RETURN
  invisible(x.secondaries)

}  







