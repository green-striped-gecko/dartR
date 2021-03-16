#' Report summary of the slot $other$loc.metrics
#'
#' This script uses any field with numeric values stored in $other$loc.metrics to produce summary statistics (mean, minimum, average, percentiles), histograms and boxplots to assist the decision of choosing thresholds for the filter function gl.filter.locmetric().
#' The fields that are included in dartR, and a short description, are found below. Optionally, the user can also set his/her own field by adding a vector into $other$loc.metrics as shown in the example. You can check the names of all available loc.metrics via: names(gl$other$loc.metrics).
#' 
#' - SnpPosition - position (zero is position 1) in the sequence tag of the defined SNP variant base 
#' - CallRate - proportion of samples for which the genotype call is non‐ missing (that is, not “‐” ) 
#' - OneRatioRef - proportion of samples for which the genotype score is 0 
#' - OneRatioSnp - proportion of samples for which the genotype score is 2 
#' - FreqHomRef - proportion of samples homozygous for the Reference allele 
#' - FreqHomSnp - proportion of samples homozygous for the Alternate (SNP) allele 
#' - FreqHets - proportion of samples which score as heterozygous, that is, scored as 1 
#' - PICRef - polymorphism information content (PIC) for the Reference allele 
#' - PICSnp - polymorphism information content (PIC) for the SNP 
#' - AvgPIC - average of the polymorphism information content (PIC) of the Reference and SNP alleles 
#' - AvgCountRef - sum of the tag read counts for all samples, divided by the number of samples with non‐zero tag read counts, for the Reference allele row 
#' - AvgCountSnp - sum of the tag read counts for all samples, divided by the number of samples with non‐zero tag read counts, for the Alternate (SNP) allele row 
#' - RepAvg - proportion of technical replicate assay pairs for which the marker score is consistent 
#' - rdepth - read depth
#' 
#' @param x -- name of the genlight object containing the SNP or presence/absence (SilicoDArT) data [required]
#' @param metric -- name of the metric to be used for filtering [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return returns a tabulation of locmetrics against different thresholds
#' @export
#' @author Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' # adding dummy data
#' test <- testset.gl
#' test$other$loc.metrics$test <- 1:nLoc(test)
#' # SNP data
#' out <- gl.report.locmetric(test,metric="test")
#' 
#' # adding dummy data
#' test.gs <- testset.gs
#' test.gs$other$loc.metrics$test <- 1:nLoc(test.gs)
#' # Tag P/A data
#' out <- gl.report.locmetric(test.gs,metric="test")

gl.report.locmetric <- function(x, metric, verbose=NULL) {
  
# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
# STANDARD ERROR CHECKING

  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
    if (all(x@ploidy == 1)){
      cat("  Processing Presence/Absence (SilicoDArT) data\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)!\n")
    }

  # Check for monomorphic loci

  if (!x@other$loc.metrics.flags$monomorphs) {
    if(verbose >= 1){cat("  Warning: genlight object contains monomorphic loci\n")}
  }

  # FUNCTION SPECIFIC ERROR CHECKING
# check whether the field exists in the genlight object
  if (!(metric %in% colnames(x$other$loc.metrics)) ) {
    stop("  Fatal Error: name of the metric not found\n")
  }
    if (!is.numeric(unlist(x$other$loc.metrics[metric]))) {
    stop("  Fatal Error: metric is not numeric\n")
  }

# DO THE JOB
  
  field <- which(colnames(x@other$loc.metrics) == metric)
  
  p1 <- p2 <- NULL
  
    # Plot Box-Whisker plot
    if (all(x@ploidy==2)){
      title1 <- paste0("SNP data - ",metric, " by Locus")
    } else {
      title1 <- paste0("Fragment P/A data - ", metric, " by Locus")
    }
  
  metric_df <- data.frame(x$other$loc.metrics[field])
  colnames(metric_df) <- "field"
  
    p1 <- ggplot(metric_df, aes(y=field)) +
      geom_boxplot() +
      theme()+
      coord_flip() +
      ylab(metric) +
      ggtitle(title1)

  p2 <-  ggplot(metric_df, aes(x=field)) +
    geom_histogram(bins = 50) +
    xlab(metric) +
    ggtitle(title1)

  
  gridExtra::grid.arrange(p1,p2)
  
  # Print out some statistics
    cat("  Reporting", metric, "Call Rate by Locus\n")
    cat("  No. of loci =", nLoc(x), "\n")
    cat("  No. of individuals =", nInd(x), "\n")
    cat("    Minimum ", metric, ": ",round(min(metric_df$field),2),"\n")
    cat("    Maximum ", metric, ": ",round(max(metric_df$field),2),"\n")
    cat("    Average ", metric, ": ",round(mean(metric_df$field),3),"\n")
    cat("    Missing ", metric, "Overall : ",round(sum(is.na(as.matrix(x)))/(nLoc(x)*nInd(x)),2),"\n")

  # Determine the loss of loci for a given filter cut-off
    percentile <- quantile(metric_df$field,probs = seq(0,1,1/20))
    retained <- unlist(lapply(percentile, function(y){res <- length(metric_df$field[metric_df$field>=y]) }))
    pc.retained <- round(retained*100/nLoc(x),1)
    filtered <- nLoc(x) - retained
    pc.filtered <- 100 - pc.retained
    df <- cbind(percentile,retained,pc.retained,filtered,pc.filtered)
    df <- data.frame(df)
    colnames(df) <- c("Threshold", "Retained", "Percent", "Filtered", "Percent")
    df <- df[order(-df$Threshold),]
    rownames(df) <- NULL

 # FLAG SCRIPT END
    
    if (verbose >= 1) {
      cat("Completed:",funname,"\n")
    }
    
    return(list(metric=df, boxplot=p1, hist=p2))

}
