#' @name gl.pcoa.plot
#' @title  Plots of the results of an ordination generated using gl.pcoa()
#' @description 
#' This script takes the output from the ordination generated by \code{\link{gl.pcoa}}
#' and plots the individuals classified by population.
#'
#' @param glPca Name of the PCA or PCoA object containing the factor scores and 
#' eigenvalues [required].
#' @param x Name of the genlight object or fd object containing the SNP genotypes 
#' or a genlight object containing the Tag P/A (SilicoDArT) genotypes or the 
#' Distance Matrix used to generate the ordination [required].
#' @param scale Flag indicating whether or not to scale the x and y axes in 
#' proportion to \% variation explained [default FALSE].
#' @param ellipse Flag to indicate whether or not to display ellipses to 
#' encapsulate points for each population [default FALSE].
#' @param p Value of the percentile for the ellipse to encapsulate points for 
#' each population [default 0.95].
#' @param labels Flag to specify the labels to be added to the plot 
#' ["none"|"ind"|"pop", default = "pop"].
#' @param as.pop Assign another metric to represent populations for the plot 
#' from the slot $other$ind.metrics [default NULL].
#' @param overlaps Factor to exclude overlapping labels, higher values allow to 
#' overlap more labels [default 3].
#' @param interactive_plot Interactive plot of the results of a PCoA ordination 
#' using the package plotly [default FALSE].
#' @param three_D_plot 3D interactive plot of the results of a PCoA ordination
#' using the package plotly [default FALSE]. 
#' @param xaxis Identify the x axis from those available in the ordination
#' [default 1].
#' @param yaxis Identify the y axis from those available in the ordination
#' [default 2].
#' @param zaxis Identify the z axis from those available in the ordination, 
#' only used for three_D_plot = TRUE [default 3].
#' @param radius Size of the points [default 3].
#' @param plot.out If TRUE, returns a plot object compatible with ggplot, 
#' otherwise returns a dataframe [default TRUE].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param palette_discrete A discrete palette for the color of populations or a 
#' list with as many colors as there are populations in the dataset
#'  [default discrete_palette].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, 
#' progress log ; 3, progress and results summary; 5, full report 
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @details 
#'   
#' The factor scores are taken from the output of \code{\link{gl.pcoa}} and the population 
#' assignments are taken from from the original data file. The specimens are 
#' shown in a bivariate plot optionally with adjacent labels and enclosing 
#' ellipses.
#'
#' Any pair of axes can be specified from the ordination, provided they are 
#' within the range of the nfactors value provided to \code{\link{gl.pcoa}}, 
#' which can be defined with the parameters xaxis, yaxis and zaxis. 
#' 
#' Axes can be scaled to represent the proportion of variation explained (scale=TRUE).
#' In any case, the proportion of variation explained by each axis is provided 
#' in the axis label.
#'
#' Points displayed in the ordination can be identified if interactive_plot = TRUE,
#' in which case the resultant plot is a plotly plot. Identification of points 
#' is by moving the mouse over them. Refer to the plotly package for further information. 
#'
#' If plot.out=TRUE, returns an object of class ggplot so that layers can 
#' subsequently be added; if plot.out=FALSE, returns a dataframe with the 
#' individual labels, population labels and PCOA scores for subsequent plotting 
#' by the user with ggplot or other plotting software.
#' 
#' Plot themes can be obtained from \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' Resultant ggplots and the tabulation are saved to the session's temporary 
#' directory.
#' 
#' @return  A plot of the ordination [plot.out = TRUE] or a dataframe [plot.out = FALSE]
#' @author Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' gl_test <- bandicoot.gl
#' # five populations in gl_test
#' nPop(gl_test)
#' # color list for population colors 
#' pop_colours <- c("deepskyblue","green","gray","orange","deeppink")
#' # Ordination applied to genotypes
#' pca <- gl.pcoa(gl_test)
#' # interactive plot to examine labels
#' pca_plot_interactive <- gl.pcoa.plot(pca, gl_test, ellipse = TRUE, labels = "ind", interactive_plot = TRUE, palette_discrete = pop_colours)
#' # 3D interactive plot
#' pca_plot_3D_interactive <- gl.pcoa.plot(pca, gl_test, three_D_plot = TRUE, palette_discrete = pop_colours)
#' # using sex as population 
#' names(gl$other$ind.metrics)
#' pca_plot_sex <- gl.pcoa.plot(glPca=pca, gl_test, as.pop = "assigned.sex", palette_discrete = pop_colours)
#'
#' @seealso \code{\link{gl.pcoa}}
#' @family Exploration and visualisation functions
#' 
#' @export
#'  

gl.pcoa.plot <- function(glPca, 
                         x, 
                         scale = FALSE, 
                         ellipse = FALSE, 
                         p = 0.95, 
                         labels = "pop",
                         as.pop = NULL,
                         overlaps = 3,
                         interactive_plot = FALSE,
                         three_D_plot = FALSE,
                         xaxis = 1,  
                         yaxis = 2, 
                         zaxis = 3,
                         radius  = 3,
                         plot.out = TRUE, 
                         plot_theme = theme_dartR(),
                         palette_discrete = discrete_palette,
                         verbose = NULL) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(f=funname,build="Jackson",v=verbose)
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(x)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  if(interactive_plot == T | three_D_plot == T){
    pkg <- "plotly"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      stop(error("Package ",pkg," needed for this function to work. Please install it."))
    }
  }
  
  pkg <- "ggrepel"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error("Package ",pkg," needed for this function to work. Please install it."))
  } 
  
  pkg <- "ggthemes"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error("Package ",pkg," needed for this function to work. Please install it."))
  }
  
  pkg <- "dplyr"
  if (!(requireNamespace(pkg, quietly = TRUE))){
    stop("Package ",pkg," needed for this function to work. Please install it.")
  }

  if(class(glPca)!="glPca") {
    stop(error("Fatal Error: glPca object required as primary input (parameter glPca)!\n"))
  }
  
  if(class(x) != "genlight" && class(x) != "dist" && class(x)  != "fd") {
    stop(error("Fatal Error: genlight, fd or dist object required as second input (parameter x)!\n"))
  }
  
  if(class(x)=="fd"){
    x <- x$gl
  }
  
  if (labels != "none" && labels != "ind" && labels != "pop"){
    cat(warn("  Warning: Parameter 'labels' must be one of 'none'|'ind'|'pop', set to 'pop'\n"))
    labels <- "pop"
  }
  
  if (labels=="ind" && class(x)=="dist"){
    cat(warn("  Warning: Individual labels cannot be applied to a PCoA based on distances between populations\n  Using population labels\n"))
    labels="pop"
  }
  
  if (p < 0 | p > 1){
    cat(warn("  Warning: Parameter 'p' must fall between 0 and 1, set to 0.95\n"))
    p <- 0.95
  }
  
  if (xaxis < 1 | xaxis > ncol(glPca$scores)){
    cat(warn("  Warning: X-axis must be specified to lie between 1 and the number of retained dimensions of the ordination",ncol(glPca$scores),"; set to 1\n"))
    xaxis <- 1
  }
  
  if (xaxis < 1 | xaxis > ncol(glPca$scores)){
    cat(warn("  Warning: Y-axis must be specified to lie between 1 and the number of retained dimensions of the ordination",ncol(glPca$scores),"; set to 2\n"))
    yaxis <- 2
  }
  
  if (xaxis < 1 | xaxis > ncol(glPca$scores)){
    cat(warn("  Warning: Z-axis must be specified to lie between 1 and the number of retained dimensions of the ordination",ncol(glPca$scores),"; set to 3\n"))
    zaxis <- 3
  }
  
  if (class(glPca) == "dist" && !is.null(as.pop)){
    cat(warn("  Warning: Temporary reassignment of population assignment not available for distance matrices\n"))
    as.pop <- NULL
  }
  
  # Set a population if none is specified (such as if the genlight object has been generated manually)
  if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
    if (verbose >= 2){ 
      cat(important("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n"))
    }
    pop(x) <- array("pop1",dim = nLoc(x))
    pop(x) <- as.factor(pop(x))
  }
  
  # DO THE JOB 
  
  # Assign the new population list if as.pop is specified
  if (class(x) == "genlight"){
    if (!is.null(as.pop)){    
      if(as.pop %in% names(x@other$ind.metrics)){
        pop(x) <- as.matrix(x@other$ind.metrics[as.pop])
        if (verbose >= 2) {
          cat(report("  Setting population assignments to",as.pop,"as specified by the as.pop parameter\n"))
          }
      } else {
        stop(error("Fatal Error: individual metric assigned to 'pop' does not exist. Check names(gl@other$loc.metrics) and select again\n"))
      }
    }
  }
  
  # Create a dataframe to hold the required scores
  df <- as.data.frame(cbind(glPca$scores[,xaxis],glPca$scores[,yaxis],glPca$scores[,zaxis]))
  
  # Convert the eigenvalues to percentages
    s <- sum(glPca$eig[glPca$eig >= 0])
    e <- round(glPca$eig*100/s,1)
    
  # Labels for the axes and points
    if(class(x)=="genlight"){
      ind <- indNames(x)
      pop_name <- factor(pop(x))
      } else {
      ind <- rownames(as.matrix(x))
      pop_name <- ind
      }  
    
    PCoAx <- PCoAy <- PCoAz <- NA
    xlab <- paste("PCA Axis", xaxis, "(",e[xaxis],"%)")
    ylab <- paste("PCA Axis", yaxis, "(",e[yaxis],"%)")
    zlab <- paste("PCA Axis", zaxis, "(",e[zaxis],"%)")
    df <- cbind(df,ind,pop_name)
    colnames(df) <- c("PCoAx","PCoAy","PCoAz","ind","pop")
    
    # assigning colors to populations
    if(class(palette_discrete)=="function"){
      colours_pops <- palette_discrete(length(levels(pop(x))))
    }
    
    if(class(palette_discrete)!="function"){
      colours_pops <- palette_discrete
    }
    
    names(colours_pops) <- as.character(levels(x$pop))
    
    ########### for labels "pop"
    if (labels == "pop" & plot.out == T & three_D_plot == F) {
      
      if (class(x)=="genlight" & verbose>=2){
        cat(report("  Plotting individuals labelled by population in a PCA with loci as attributes\n"))
      }
      
      if(!class(x)=="genlight" & verbose>=2){
        cat(report("  Plotting populations in a PCoA with loci as attributes based on a distance matrix\n"))
      }
      
      # Plot
      suppressWarnings(
      p1 <- ggplot(df, aes(x = PCoAx, y = PCoAy, group = pop, colour = pop)) +
            geom_hline(yintercept = 0) +
            geom_vline(xintercept = 0) +
            geom_point(size = radius) +
            ggrepel::geom_label_repel(aes(label = pop),show.legend = FALSE, max.overlaps = overlaps) +
            labs(x = xlab, y = ylab) +
            plot_theme +
            theme(axis.title = element_text(size = 16)) +
            scale_color_manual(name = "Populations", values = colours_pops)
      )
          }
  
    ########### for labels "ind"
    if (labels == "ind" & plot.out == T & three_D_plot == F) {
      
      if (class(x)=="genlight" & verbose>=2){
        cat(report("  Plotting and labelling individuals in a PCA with loci as attributes\n"))
      }
      
      if(!class(x)=="genlight" & verbose>=2){
        cat(warn("  Warning: Data on individuals not available when applying PCoA to a distance matrix\n"))
        cat(report("    Displaying populations as entities in a space defined by allele frequencies\n"))
      }
      
      # Plot
      suppressWarnings(
        p1 <- ggplot(df, aes(x = PCoAx, y = PCoAy, group = pop, colour = pop, label = ind)) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          geom_point(size = radius) +
          ggrepel::geom_label_repel(aes(label = ind),show.legend = FALSE, max.overlaps = overlaps) +
          labs(x = xlab, y = ylab) +
          plot_theme +
          theme(axis.title = element_text(size = 16)) +
          scale_color_manual(name = "Populations", values = colours_pops)
      )
    }
    
    ########### for labels "none"
    if (labels == "none" | labels==FALSE & plot.out == T  & three_D_plot == F) {
      
      if (class(x)=="genlight" & verbose>=2){
        cat(report("  Plotting individuals without labels in a PCA with loci as attributes\n"))
      }
      
      if(!class(x)=="genlight" & verbose>=2){
        cat(report("  Plotting populations in a PCoA with loci as attributes based on a distance matrix\n"))
      }
      
      # Plot
      suppressWarnings(
        p1 <- ggplot(df, aes(x = PCoAx, y = PCoAy, group = pop, colour = pop)) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          geom_point(size = radius) +
          labs(x = xlab, y = ylab) +
          plot_theme +
          theme(axis.title = element_text(size = 16)) +
          scale_color_manual(name = "Populations", values = colours_pops)
      )
    }
    
    # Scale the axes in proportion to % explained, if requested
    if(scale==TRUE  & plot.out == T  & three_D_plot ==F) { 
      suppressWarnings(
      p1 <- p1 + coord_fixed(ratio=e[yaxis]/e[xaxis]) 
      )
    }
    
    # Add ellipses if requested
    if(ellipse==TRUE & plot.out == T  & three_D_plot ==F) {
      suppressWarnings(
      p1 <- p1 + stat_ellipse(aes(colour=pop), type="norm", level=p)
      )
    }
    
    ########### for interactive_plot
    if(interactive_plot == T  & plot.out == T  & three_D_plot ==F){
      if (verbose>=2){
      cat(report("  Displaying an interactive plot, mouse over for details for each point\n"))
      }
      suppressWarnings(
      p1 <- plotly::ggplotly(p1)
      )
    }
    
    ########### for three_D_plot
    if(three_D_plot == T & plot.out == T){
      if (verbose>=2){
        cat(report("  Displaying a three dimensional plot, mouse over for details for each point\n"))
      }
      suppressWarnings(
        p1 <- plotly::plot_ly(df,x=~PCoAx,y=~PCoAy,z=~PCoAz,
                              marker = list(size = radius*2),colors = colours_pops, text=ind) %>% 
          plotly::add_markers(color=~pop)%>% 
          plotly::layout(legend=list(title=list(text='Populations')),
                         scene = list(xaxis = list(title = xlab,titlefont = list(size = 16)),
                                      yaxis = list(title = ylab,titlefont = list(size = 16)),
                                      zaxis = list(title = zlab,titlefont = list(size = 16))))
      )
    }

    # PRINTING OUTPUTS
    if(plot.out == T){
      suppressWarnings(print(p1))
    }
    
    if(verbose>=2){
      cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))
    } 

    # SAVE INTERMEDIATES TO TEMPDIR    
    match_call <- paste0(names(match.call()),"_",as.character(match.call()),collapse = "_")
    # creating temp file names
    if(plot.out == T){
      temp_plot <- tempfile(pattern = "dartR_plot_")
    }
    temp_table <- tempfile(pattern = "dartR_table_")
    
    # saving to tempdir
    if(plot.out == T){
      saveRDS(list(match_call,p1), file = temp_plot)
      if(verbose>=2){
      cat(report("  Saving the ggplot to session tempfile\n"))
    }
    }
    
    saveRDS(list(match_call,df), file = temp_table)
    if(verbose>=2){
      cat(report("  Saving tabulation to session tempfile\n"))
    }
    
    if(class(x)=="dist"){
      df <- data.frame(id=labels(x),glPca$scores)
    } else {
      df <- data.frame(id=indNames(x), pop=as.character(pop(x)), glPca$scores)
      row.names(df) <- NULL
    }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
      cat(report("\nCompleted:", funname, "\n\n"))
    }
    
    # RETURN

  if(plot.out) {
    suppressWarnings(invisible(p1))
  } else {
    invisible(df)
  }
  
}

