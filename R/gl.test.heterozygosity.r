#' @name gl.test.heterozygosity
#'
#' @title Tests the difference in heterozygosity between populations taken pairwise
#'
#' @description 
#' Calculates the expected heterozygosities for each population in a genlight object,
#' and uses re-randomization to test the statistical significance of differences
#' in heterozygosity between populations taken pairwise.
#' 
#' @param x A genlight object containing the SNP genotypes [required]
#' @param nreps Number of replications of the re-randomization [default 1,000]
#' @param alpha1 First significance level for comparison with diff=0 on plot [default 0.05]
#' @param alpha2 Second significance level for comparison with diff=0 on plot [default 0.01]
#' @param plot If TRUE, plots a sampling distribution of the differences for each comparison [default TRUE]
#' @param max_plots Maximum number of plots to print per page [default 9]
#' @param plot_theme Theme for the plot. See Details for options [default theme_dartR()].
#' @param plot_colours List of two color names for the borders and fill of the
#'  plots [default two_colors].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity]
#'
#' @details 
#'  
#'\strong{ Function's output }
#'
#' If plot = TRUE, plots are created showing the sampling distribution for the difference between each pair
#' of heterozygosities, marked with the critical limits alpha1 and alpha2, the observed
#' heterozygosity, and the zero value (if in range).
#'  Plots and table are saved to the temporal directory (tempdir) and can be accessed with the function \code{\link{gl.print.reports}} and listed with the function \code{\link{gl.list.reports}}. Note that they can be accessed only in the current R session because tempdir is cleared each time that the R session is closed.
#'   
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return A dataframe containing population labels, heterozygosities and sample sizes
#'
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' out <- gl.test.heterozygosity(bandicoot.gl, nreps=10, verbose=3, plot=TRUE)
#'
#' @family Genetic variation within populations
#'
#' @import patchwork
#'
#' @export
#'  

gl.test.heterozygosity <- function(x, 
                                   nreps = 100,
                                   alpha1 = 0.05,
                                   alpha2 = 0.01,
                                   plot = TRUE,
                                   max_plots = 9,
                                   plot_theme = theme_dartR(), 
                                   plot_colours = two_colors, 
                                   verbose = NULL) { 
  
  
  # TRAP COMMAND
  
  funname <- match.call()[[1]]
  
  # SET VERBOSITY
  
  verbose <- gl.check.verbosity(verbose)
  
  # CHECKS DATATYPE 
  
  datatype <- utils.check.datatype(x)

  # SCRIPT SPECIFIC ERROR CHECKING
  
  # Upper and lower significance boundaries
  
  if(alpha1 < 0 || alpha1 > 1) {
    cat(warn("Warning: First alpha value should be between 0 and 1, set to 0.05\n"))
    alpha1 <- 0.05
  }
  if(alpha2 < 0 || alpha1 > 2) {
    cat(warn("Warning: Second alpha value should be between 0 and 1, set to 0.01\n"))
    alpha1 <- 0.01
  }
  
  upper1 <- 1-alpha1 # significance level #1
  upper2 <- 1-alpha2 # significance level #2
  if (upper1 > upper2){
    tmp <- upper2
    upper2 <- upper1
    upper1 <- tmp
  }
  lower1 <- 1-upper1
  lower2 <- 1-upper2
  
  #bon.critical <- round(lower1/((nPop(x)*nPop(x)-nPop(x))/2),4)
  
  # Set a population if none is specified (such as if the genlight object has been generated manually)
  if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
    if (verbose >= 2) {
      cat(warn("  No population assignments detected, 
                             individuals assigned to a single population labelled 'pop1'\n"))
    }
    pop(x) <- array("pop1", dim = nInd(x))
    pop(x) <- as.factor(pop(x))
  }
  
  # Check for monomorphic loci
  
  if (x@other$loc.metrics.flags$monomorphs==FALSE) {
    if(verbose >= 1){
      cat(warn("  Warning: genlight object contains monomorphic loci which will be factored into heterozygosity estimates\n"))
    }
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
  
  # Calculate a matrix of heterozygosities (He)
  out <- utils.het.pop(x)
  D <- outer(out,out,"-")
  
  # Simulate the distribution of differences in He by population, pairwise
  # Initialize the data matrix to hold simulation results
  if (verbose >= 2){
    cat(report("  Calculating the sampling distributions for pairwise differences between populations by re-randomization\n"))
    cat(important("  Please be patient .... go have a coffee\n"))
  }
  mat <- array(data=NA, dim=c(nPop(x),nPop(x),nreps))
  n <- nLoc(x)
  for (i in 1:nreps) {
    # Sample the dataset
    x2 <- x[,sample(1:n,n,replace=T)]
    # Calculate Heterozygosity
    out <- utils.het.pop(x2)
    # Save difference values away
    mat[,,i] <- outer(out,out,"-")
  }
  
  # Calculate the p values, significance for pairwise differences  
  # Dataframe to store output
  df <- data.frame(matrix(ncol = 5, nrow = (nPop(x)*nPop(x)-nPop(x))/2))
  colnames(df) <- c('pop1','pop2','diff','significance','pval')
  
  count <- 0
  
  if (plot){
    if (verbose >= 2){
      cat(report("  Cycling through the",(nPop(x)*nPop(x)-nPop(x))/2,"pairs of populations\n"))
    }
    if (verbose >= 2){
      cat(report("  Plotting sampling distributions\n"))
    }
    
    # set the limit values for the x axis in the plots
    x_axis_limits_lots <- c(min(mat),max(mat))
    # count how many plots are going to be created
    total_number_plots <- choose(nPop(x),2)
    # maximum number of plots per page
    max_plots <- max_plots
    # create list to contain plots
    p_list <- list()
    
    for (first in 1:(nPop(x)-1)){
      for (second in (first+1):nPop(x)) {
        count <- count + 1 
        
        # Prepare the parameters for the plot
        
        # Observed He
        obs <- D[first,second]
        
        # Upper and lower significance limits
        u1quantile <- as.numeric(quantile(mat[first,second,],upper1))
        u2quantile <- as.numeric(quantile(mat[first,second,],upper2))
        l1quantile <- as.numeric(quantile(mat[first,second,],lower1))
        l2quantile <- as.numeric(quantile(mat[first,second,],lower2))
        
        # Is zero within the range specified by the upper and lower limits
        if (0 < u2quantile && 0 > l2quantile){
          signif <- paste0("non-sig @",lower2)
        }
        if (0 < u1quantile && 0 > l1quantile){
          signif <- paste0("non-sig @",lower1)
        }
        if (0 > u1quantile || 0 < l1quantile){
          signif <- paste0("sig @",lower1)
        }
        if (0 > u2quantile || 0 < l2quantile){
          signif <- paste0("sig @",lower2)
        }
        
        # p value for zero
        values <- mat[first,second,]
        p <- min(length(values[values>0]),length(values[values<=0]))/length(abs(values))
        p <- round(p,4)
        
        # Store results in a df  
        df$pop1[count] <- popNames(x)[first]
        df$pop2[count] <- popNames(x)[second]
        df$diff[count] <- obs
        df$significance[count] <- signif
        df$pval[count] <- p
        
        # Plot the histogram of pairwise differences
        
        # Construct the label
        title <- paste0(popNames(x)[first]," vs ",popNames(x)[second])
        subtitle <- paste0(signif," (p=",p,")")
        
        plot_values <- as.data.frame(values)
        colnames(plot_values) <- "values"
        
        # these plots do not contain legends
        if(count < total_number_plots | count != max_plots ){
        
        p_temp <- ggplot(plot_values, aes(x = values)) + 
          geom_histogram(bins = 50, color = plot_colours[1],fill = plot_colours[2]) + 
           coord_cartesian(xlim = x_axis_limits_lots) + 
          xlab("Difference") + 
          ylab("Count") + 
          plot_theme +
          theme(plot.title = element_text(size = 12))+
          labs(title=title,subtitle = subtitle) + 
          # Add lines for the observed value of He, Zero, upper and lower levels of significance
          geom_vline(xintercept = 0,colour= "blue",size=1) +
          geom_vline(xintercept=D[first,second],colour="green",size=1)+
          geom_vline(xintercept=u1quantile,colour="firebrick4",size=1)+
          geom_vline(xintercept=u2quantile,colour="firebrick1",size=1)+
          geom_vline(xintercept=l1quantile,colour="firebrick4",size=1)+
          geom_vline(xintercept=l2quantile,colour="firebrick1",size=1)
        
        p_list[[count]] <- p_temp
        
        }
        
        # this plot contains the plot legends. it is just created every max_plots and in the last plot
        if(count == total_number_plots | (count %% max_plots) == 0 ){
          
          p_temp <- ggplot(plot_values, aes(x = values)) + 
            geom_histogram(bins = 50, color = plot_colours[1],fill = plot_colours[2]) + 
            coord_cartesian(xlim = x_axis_limits_lots) + 
            xlab("Difference") + 
            ylab("Count") + 
            plot_theme +
            theme(plot.title = element_text(size = 12))+
            labs(title=title,subtitle = subtitle)+
            geom_vline(aes(xintercept = 0,colour= "Zero_value"),size=1) +
            geom_vline(aes(xintercept=D[first,second],colour="Observed"),size=1)+
            geom_vline(aes(xintercept=u1quantile,colour="alpha1"),size=1)+
            geom_vline(aes(xintercept=u2quantile,colour="alpha2"),size=1)+
            geom_vline(xintercept=l1quantile,colour="firebrick4",size=1)+
            geom_vline(xintercept=l2quantile,colour="firebrick1",size=1)+
            scale_color_manual(name = "Values" ,values=c( Zero_value="blue",Observed="green",alpha1="firebrick4",alpha2="firebrick1"),labels=c(paste("Sig. ",alpha1),paste("Sig. ",alpha2),"Observed","Zero value"))+
            guides(color = guide_legend(override.aes = list(size = 5),ncol = 2))+
            theme(legend.position="bottom",legend.title = element_text(face = "bold"))
           
          p_list[[count]] <- p_temp
        }

      }
      
    }
    
    # PRINTING OUTPUTS
    match_call <- paste0(names(match.call()),"_",as.character(match.call()),collapse = "_")
    # using package patchwork
    for (i in 1:ceiling((length(p_list)/max_plots))) {
      seq_1 <- seq(1,length(p_list),max_plots)
      seq_2 <- seq(1,length(p_list),max_plots)-1
      seq_2 <- seq_2[-1]
      seq_2 <- c(seq_2,seq_1[length(seq_1)])
      
      p_final <- wrap_plots(p_list[seq_1[i]:seq_2[i]],ncol = 3)
  
      print(p_final)
      # SAVE INTERMEDIATES TO TEMPDIR             
      # creating temp file names
      temp_plot <- tempfile(pattern =paste0("dartR_plot_",seq_1[i],"_to_",seq_2[i]))
      # saving to tempdir
      saveRDS(list(match_call,p_final), file = temp_plot)
    }
    if(verbose>=2){
      cat(report("  Saving the ggplot to session tempfile\n"))
    }
  }

  print(df)
  
  # SAVE INTERMEDIATES TO TEMPDIR             
  # creating temp file names
  temp_table <- tempfile(pattern = "dartR_table_")
  match_call <- paste0(names(match.call()),"_",as.character(match.call()),collapse = "_")
  # saving to tempdir
  saveRDS(list(match_call,df), file = temp_table)
  if(verbose>=2){
    cat(report("  Saving tabulation to session tempfile\n"))
  }
  
  if(verbose>=2){
    cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))
  } 
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("\n\nCompleted:", funname, "\n\n"))
  }
  
  # RETURN
  
  invisible(df) 
}


