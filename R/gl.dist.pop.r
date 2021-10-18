#' @name gl.dist.pop
#' @title Calculate a distance matrix for populations with SNP genotypes in a genlight object
#' @description
#' This script calculates various distances between populations based on allele frequencies. The distances are
#' calculated by scripts in the {stats} or {vegan} libraries, with the exception of the pcfixed (percent fixed
#' differences) distance.
#' @details
#' The distance measure can be one of "manhattan", "euclidean", "pcfixed", "pa", canberra", "bray", 
#' "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , 
#' "binomial", "chao", "cao", "mahalanobis", "maximum", "binary" or "minkowski". Refer to the documentation for
#' of functions \link[stats]{dist} (package stat) or \link[vegan]{vegdist} (package vegan) vegan for definitions. 
#' 
#' Distance pcfixed calculates the pair-wise count of fixed allelic differences between populations.
#'
#' @param x Name of the genlight containing the SNP genotypes [required]
#' @param method Specify distance measure [default euclidean]
#' @param plot.out If TRUE, display a histogram of the genetic distances, and a whisker plot [default TRUE]
#' @param binary Perform presence/absence standardization before analysis using decostand [default FALSE]
#' @param p The power of the Minkowski distance (typically a value ranging from 0.25 to infinity) [default 0.5]
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot_colors Vector with two color names for the borders and fill [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 
#' 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return An object of class 'dist' giving distances between populations
#' @importFrom stats dist
#' @importFrom vegan vegdist
#' @export
#' @author Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' D <- gl.dist.pop(testset.gl, method="euclidean")

gl.dist.pop <- function(x, 
                        method="euclidean", 
                        plot.out=TRUE, 
                        binary=FALSE, 
                        p=NULL, 
                        plot_theme = theme_dartR(), 
                        plot_colors = two_colors, 
                        save2tmp = FALSE,
                        verbose=NULL) {

# CHECK IF PACKAGES ARE INSTALLED
  pkg <- "reshape2"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package",pkg," needed for this function to work. Please install it.")
  }
  
# SET VERBOSITY
      verbose <- gl.check.verbosity(verbose)
      
# FLAG SCRIPT START
      funname <- match.call()[[1]]
      utils.flag.start(func=funname,build="Jody",verbosity=verbose)
      
# CHECK DATATYPE 
      datatype <- utils.check.datatype(x,accept="SNP",verbose=verbose)
      
# FUNCTION SPECIFIC ERROR CHECKING

# DO THE JOB

  veganmethod <- c("bray", 
                   "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , 
                   "binomial", "chao", "cao", "mahalanobis")
  distmethod <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "Chebyshev")

  if (!(method%in%veganmethod || method%in%distmethod || method=="pcfixed")){
    stop(error("Fatal Error: Specified distance method is not among those available.\n  Specify one of",paste(paste(veganmethod,collapse = ", "), paste(distmethod,collapse = ", "),collapse = ", "),"or pcfixed.\n"))
  }
  hard.min.p <- 0.25

  m <- method
  b <- binary
  pr <- p
  u <- TRUE
  d <- TRUE
  
  # Calculate allele frequencies for each population and locus
    f <- gl.percent.freq(x,verbose=0)
  # Select only pop, locus, frequency columns  
    f <- f[,c(1,2,6)]
  # Convert to a pop x locus matrix
    f <- reshape2::dcast(f, popn ~ locus, value.var="frequency")
  # Reassign names to the populations, and convert from percentages to proportions 
    row.names(f)=f[,1]
    f <- f[,-c(1)]
    f <- f/100
    
  # Calculate distance using dist {stat}
    if (m %in% distmethod) {
      if (verbose >= 2) {
        cat(paste("  Calculating distances: ",m,"\n"))
        cat("  Refer to dist {stats} documentation for algorithm\n")
      }  
      if (method == "minkowski"){
        if (pr < 0.25) {
          if (verbose >= 2){cat("  Warning:",hard.min.p,"is the practical minimum for Minkowski distance, set to,",hard.min.p,"\n\n")}
          pr <- hard.min.p
        }
        if (pr == 1) {
          if (verbose >= 2) {cat("  Note: for p = 1, Minkowski distance is equivalent to Manhattan distance\n\n")}
        }
        if (pr == 2) {
          if (verbose >= 2) {cat("  Note: for p = 2, Minkowski distance is equivalent to Euclidean distance\n\n")}
        }
        if (pr >= 30) {
          if (verbose >= 2) {cat("  Note: for large p, Minkowski distance is equivalent to the Maxiumum Metric distance\n\n")}
        }
        if (pr < 1) {
          if (verbose >= 2) {cat("  Note: for p < 1, Minkowski distance is not a metric distance, and so should be considered a measure of dissimilarity\n\n")}
        }
      }
      dd <- stats::dist(f, method=m, diag=d, upper=u, p=pr)
    }
    
    # Calculate distance using vegdist {vegan}
    if (m %in% veganmethod) {
      dd <- vegan::vegdist(f, method=m, binary=b, diag=d, upper=u, na.rm=TRUE)
      if (verbose >= 2) {
        cat(report(paste("  Calculating distances: ",m,"\n")))
        cat(report("    Refer to vegdist {vegan} documentation for algorithm\n"))
      }
      if (method == "bray"){
        if (verbose >= 2) {cat("  Note: the Bray-Curtis distance is non-metric, and so should be considered a dissimilarity measure. A metric alternative is the Jaccard distance.\n\n")}
      }
    }
    
    if (m == "pcfixed"){
      dd <- gl.fixed.diff(x,verbose=0)[[3]]
      if (verbose >= 2) {
        cat(report("  Calculating percent fixed differences\n"))
        cat(warn("Note: this distance may be non-metric, and so should be considered a dissimilarity measure\n"))
      }  
    }
    # if (m == "pa"){
    #   dd <- gl.report.pa.pop(x)
    #   dd <- dd[,c(3,4,10)]
    #   tmp <- dd
    #   tmp2 <- dd$pop1
    #   dd$pop1 <- dd$pop2
    #   dd$pop2 <- tmp2
    #   dd <- rbind(dd,tmp)
    #   dd <- dcast(dd, pop1 ~ pop2, value.var="totalpriv")
    #   row.names(dd) <- dd[,1]
    #   dd <- dd[,2:length(dd[1,])]
    #   dd <- dd[order(row.names(dd)),]
    #   dd <- dd[,order(colnames(dd))]
    #   if (verbose >= 2) {
    #     cat("  Calculating total private alleles\n")
    #     cat("  Note: this distance may be non-metric, and so should be considered a dissimilarity measure\n")
    #   }
    # }
    
    dd <- as.dist(dd) 
    
  # # Revert to original order  
  #   ord <- rank(popNames(x))
  #   mat <- as.matrix(dd)[ord, ord]
  #   dd <- as.dist(mat)
    mat <- as.matrix(dd)
    
# PLOT
    # Plot Box-Whisker plot

  if (plot.out){
      if (datatype=="SNP"){
        title_plot <-  paste0("SNP data\nUsing ",method," distance")
      } else {
        title_plot <- paste0("Tag P/A data (SilicoDArT)\nUsing ",method," distance")
      }  
      values <- NULL
      df_plot <- data.frame(values =as.vector(mat))
      
      # Boxplot
      p1 <- ggplot(df_plot, aes(y = values)) + 
        geom_boxplot(color = plot_colors[1], fill = plot_colors[2]) + 
        coord_flip() + 
        plot_theme + 
        xlim(range = c(-1, 1)) + 
        ylim(min(df_plot$values,na.rm=TRUE),max(df_plot$values,na.rm=TRUE)) +
        ylab(" ") + 
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
        ggtitle(title_plot)
      
      # Histogram
      p2 <- ggplot(df_plot, aes(x = values)) + 
        geom_histogram(bins = 20, color = plot_colors[1], fill = plot_colors[2]) +
        xlim(min(df_plot$values,na.rm=TRUE),max(df_plot$values,na.rm=TRUE)) +
        xlab("Distance") + 
        ylab("Count") + 
        plot_theme
    }
    
  
# SUMMARY 
    # Print out some statistics
  if(verbose >= 3){
    cat("  Reporting inter-population distances\n")
    cat("  Distance measure:",method,"\n")
    cat("    No. of populations =", nPop(x), "\n")
    cat("    Average no. of individuals per population =", nInd(x)/nPop(x), "\n")
    cat("    No. of loci =", nLoc(x), "\n")
    cat("    Minimum Distance: ",round(min(dd),2),"\n")
    cat("    Maximum Distance: ",round(max(dd),2),"\n")
    cat("    Average Distance: ",round(mean(dd),3),"\n")
  }  
    
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
      saveRDS(list(match_call,dd), file = temp_table)
      if(verbose>=2){
        cat(report("  Saving tabulation to session tempfile\n"))
        cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))
      }
    }    
    
# PRINTING OUTPUTS
    if(plot.out){
      # using package patchwork
      p3 <- (p1/p2) + plot_layout(heights = c(1, 4))
      suppressWarnings(print(p3))
    }
    
# FLAG SCRIPT END
  if (verbose > 0) {
    cat(report("Completed:",funname,"\n"))
  }
    
  return(dd)
}