#' Tests the difference in heterozyosity between populations taken pairwise. 
#' 
#' Calculates the expected heterozygosities for each population in a genlight object,
#' and uses re-randomization to test the statistical significance of differences
#' in heterozyosity between populations taken pairwise.
#' 
#' Optionally plots the sampling distribution for the difference between each pair
#' of heterozygosities, marked with the critical limits alpha1 and alpha2, the observed
#' heterozygosity, and the zero value (if in range). Can output the plots to jpeg.
#' 
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param nreps -- number of replications of the re-randomization [10,000]
#' @param alpha1 -- significance level for comparison with diff=0 on plot [0.05]
#' @param alpha2 -- second significance level for comparison with diff=0 on plot [0.01]
#' @param plot -- if TRUE, plots a sampling distribution of the differences for each comparison [FALSE]
#' @param plot.out -- if TRUE, outputs the plots as jpeg [FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return returns a dataframe containing population labels, heterozygosities and sample sizes
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom grDevices rainbow
#' @importFrom graphics par
#' @importFrom plyr join
#' @importFrom robustbase adjbox
#' @examples
#' gl.report.heterozygosity(testset.gl)
#' out <- gl.test.heterozygosity(testset.gl, nreps=100, verbose=3)

gl.test.heterozygosity <- function(x, 
                                   nreps=10000,
                                   alpha1=0.05,
                                   alpha2=0.01,
                                   plot=FALSE,
                                   plot.out=FALSE,
                                   verbose=NULL) { 
 
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
    stop("Fatal Error: genlight object required!")
  }
  
  if (all(x@ploidy == 1)){
    stop("  Processing  Presence/Absence (SilicoDArT) data, heterozygosity can only be calculated for SNP data\n")
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
  
# SCRIPT SPECIFIC ERROR CHECKING
  
  # Upper and lower significance boundaries
  
  if(alpha1 < 0 || alpha1 > 1) {
    cat("Warning: First alpha value should be between 0 and 1, set to 0.05\n")
    alpha1 <- 0.05
  }
  if(alpha2 < 0 || alpha1 > 2) {
    cat("Warning: Second alpha value should be between 0 and 1, set to 0.01\n")
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
  
  # Check for monomorphic loci
  
  if (x@other$loc.metrics.flags$monomorphs==FALSE) {
    if(verbose >= 1){
      cat("  Warning: genlight object contains monomorphic loci which will be factored into heterozygosity estimates\n")
    }
  }
# DO THE JOB
  
  # Calculate a matrix of heterozyosities (He)
  out <- gl.report.heterozygosity(x, v=0)
  D <- outer(out$He,out$He,"-")

  # Simulate the distribution of differences in He by population, pairwise
    # Initialize the data matrix to hold simulation results
    if (verbose >= 2){
      cat("  Calculating the sampling distributions for pairwise differences between populations by re-randomization\n")
      cat("  Please be patient .... go have a coffee\n")
    }
    mat <- array(data=NA, dim=c(nPop(x),nPop(x),nreps))
    n <- nLoc(x)
    for (i in 1:nreps) {
      # Sample the dataset
        x2 <- x[,sample(1:n,n,replace=T)]
      # Calculate Heterozyosity
        out <- gl.report.heterozygosity(x2, verbose=0)
      # Save difference values away
        mat[,,i] <- outer(out$He,out$He,"-")
    }

  # Calculate the p values, significance for pairwise differences  
    
    # Dataframe to store output
    df <- data.frame(matrix(ncol = 5, nrow = (nPop(x)*nPop(x)-nPop(x))/2))
    colnames(df) <- c('pop1','pop2','diff','significance','pval')
    
    count <- 0
    if (verbose >= 2){cat("  Cycling through the",(nPop(x)*nPop(x)-nPop(x))/2,"pairs of populations\n")}
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
      if (plot){
      if (verbose >= 2){cat("  Plotting sampling distributions\n")}
      # Contruct the label
      label <- paste0(popNames(x)[first]," x ",popNames(x)[second]," -- ",signif," (p = ",p,")")
            a <- hist(values,
              breaks=100,
              main=label,
              #          xlim=c(min(0,min(mat[first,second,])),max(0,max(mat[first,second,]))),
              xlab="difference",
              freq=F,
              col="gray"
            )
      
      # Add lines for the observed value of He, Zero, upper and lower levels of significance
      abline(v=0,col="blue",lty=2)
      abline(v=D[first,second],col="green")
      abline(v=u1quantile,col="red")
      abline(v=u2quantile,col="red")
      abline(v=l1quantile,col="red")
      abline(v=l2quantile,col="red")
      
      # Label the lines  
      text(u1quantile,max(a$density),upper1,col='blue',cex=0.6,adj=0)
      text(u2quantile,max(a$density),upper2,col='blue',cex=0.6,adj=0)
      text(l1quantile,max(a$density),lower1,col='blue',cex=0.6,adj=1)
      text(l2quantile,max(a$density),lower2,col='blue',cex=0.6,adj=1)
      text(obs,max(a$density),"OBS",col='blue',cex=0.6, adj=0.5)
      text(0,max(a$density),"0",col='blue',cex=1, adj=1)
      
      }
      
      # Write the plots out
      if(plot.out){
        if (verbose >= 2){cat("  Outputting plots as jpeg files\n")}
        plot.name <- paste0(popNames(x)[first]," x ",popNames(x)[second],".jpg")
        dev.copy(jpeg,plot.name)
        dev.off()
      }

    }
    }
    
# FLAG SCRIPT END
    
    if (verbose > 0) {
      cat("Completed:",funname,"\n")
    }
    
    return(df) 
}
    

