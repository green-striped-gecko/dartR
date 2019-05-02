#' Reports observed and expected heterozygosity by population or by individual. 
#' 
#' Calculates the observed and expected heterozygosities for each population (method="pop") 
#' or the observed heterozyosity for each individual (method="ind") in a genlight object.
#' 
#' Observed heterozygosity for populations is the observed proportion of heterozygotes
#' averaged over individuals.
#' 
#' Output for method="pop" is an ordered barchart of observed heterozygosity across 
#' populations and a table of expected heterozygosity and estimated variance for 
#' expected heterozygosity by population.
#' 
#' Expected heterozygosity for a population is calculated by calculating the expected
#' proportion of heterozygotes, under Hardy-Weinberg equilibrium, for each locus, then
#' averaging across the loci. The calculations are according to Nei, M. (1987) Molecular 
#' evolutionary genetics. New York: Columbia University Press.
#' 
#' Output for method="ind" is a histogram of heterozygosity calculated across individuals.
#' The histogram is accompanied by a box and whisker plot presented either in standard 
#' (boxplot="standard") or adjusted form (boxplot=adjusted). 
#' 
#' Refer to Tukey (1977, Exploratory Data Analysis. Addison-Wesley) for standard
#' Box and Whisker Plots and Hubert & Vandervieren (2008, An Adjusted Boxplot for Skewed
#' Distributions, Computational Statistics & Data Analysis 52:5186-5201) for adjusted
#' Box and Whisker Plots.
#' 
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param boxplot -- if 'standard', plots a standard box and whisker plot; if 'adjusted',
#' plots a boxplot adjusted for skewed distributions [default 'adjusted']
#' @param range -- specifies the range for delimiting outliers [default = 1.5 interquartile ranges]
#' @param cex.labels -- sets the size of the population labels [default 0.7]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return a dataframe containing population labels, heterozygosities and sample sizes
#' @export
#' @author Bernd Gruber, Arthur Georges and Renee Catullo (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom grDevices rainbow
#' @importFrom graphics par
#' @importFrom plyr join
#' @importFrom robustbase adjbox
#' @examples
#' tmp <- gl.report.heterozygosity(testset.gl,verbose=3)
#' tmp <- gl.report.heterozygosity(testset.gl,method='ind',verbose=3)

# Last amended 1-May-19

gl.report.heterozygosity <- function(x, 
                                     method="pop", 
                                     boxplot="adjusted",
                                     range=1.5,
                                     cex.labels=0.7, 
                                     verbose=2) {
  
# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  
  if (!(method=="pop" | method == "ind")) {
    cat("Warning: Method must either be by population or by individual, set to method='pop'\n")
    method <- "pop"   
  }
  
  if (!(boxplot=="standard" | boxplot == "adjusted")) {
    cat("Warning: Box-whisker plots must either standard or adjusted for skewness, set to boxplot='adjusted'\n")
    boxplot <- 'adjusted'   
  }

  # Work around a bug in adegenet if genlight object is created by subsetting
    x@other$loc.metrics <- x@other$loc.metrics[1:nLoc(x),]

  # Set a population if none is specified (such as if the genlight object has been 
  # generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, 
                             individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nLoc(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
      cat("  Warning: genlight object contains monomorphic loci\n")
    }

# DO THE JOB FOR POPULATIONS
    
  if (method=="pop"){

# Split the genlight object into a list of populations
  sgl <- seppop(x)
  
# OBSERVED HETEROZYGOSITY  
  # if (verbose >=2){cat("  Calculating Observed Heterozygosities by population\n")}
  # 
  # # Calculate heterozygosity for each population in the list
  #   Ho <- data.frame(lapply(sgl, function(x) mean(colMeans(as.matrix(x, na.rm=TRUE)==1), na.rm=TRUE) ))
  #   
  # # Calculate sample sizes
  #   sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x, na.rm=TRUE)==1),
  #                                                   na.rm=TRUE) ))
  #   n <- t(sums/Ho)
  #   n <- cbind(row.names(n),n)
  #   
  # # Join the sample sizes with the heteozygosities
  #   df1 <- data.frame(pop=names(Ho), Ho=as.numeric(Ho[1,]))
  #   df2 <- data.frame(n)
  #   names(df2)<- c("pop","nInd")
  #   df <- plyr::join(df1,df2, by="pop")
  #   
  # # Plot the results
  #   op <- par(mfrow=c(1,1),mai=c(1.2,0.5,0.2,0),oma=c(2,2,2,0), pty="m")
  #   df.ordered <- df[order(df$Ho),]
  #   barplot(df.ordered$Ho, names.arg=paste(df.ordered$pop, df.ordered$nInd, sep=" | "), las=2, cex.names=cex.labels, space=0, border=F, col=rainbow(nrow(df.ordered)), 
  #           main="Observed Heterozygosity by Population")

# EXPECTED HETEROZYGOSITY
    if (verbose >=2){cat("  Calculating Expected Heterozygosities, averaged across loci, for each population\n")}
    
    Hexp <- array(NA,length(sgl))
    # For each population
    for (i in 1:length(sgl)){
      gl <- sgl[[i]]
      
      p <- array(NA,nLoc(gl))  
      q <- array(NA,nLoc(gl))
      He <- array(NA,nLoc(gl))
      
      # For each locus
      for (j in 1:nLoc(gl)){
        
        # Calculate proportion individuals for which the focal locus is heterozygous
        # First sum up the number of homozygous and heterozygous loci
        tbl <- table(as.matrix(gl)[,j],useNA = "always")
        if (is.na(tbl["0"])) {tbl["0"] <- 0}
        if (is.na(tbl["1"])) {tbl["1"] <- 0}
        if (is.na(tbl["2"])) {tbl["2"] <- 0}
        
        # Back-calculate the allele frequencies
        p[j] <- (2*tbl["0"] + tbl["1"]) # Raw
        q[j] <- (tbl["1"] + 2*tbl["2"])
        
        # Total alleles
        n <- p[j] + q[j]
        
        # Compute proportions
        p[j] <- p[j]/n # Proportions
        q[j] <- q[j]/n
        
        He[j] <- 1 - (p[j]*p[j] + q[j]*q[j])
      }
      Hexp[i] <- mean(as.numeric(He),na.rm=TRUE)
    }  
      
# OBSERVED HETEROZYGOSITY
    if (verbose >=2){cat("  Calculating Observed Heterozygosities, averaged across individuals, for each population\n")}
  
  Hobs <- array(NA,length(sgl))
  # For each population
  for (i in 1:length(sgl)){
    gl <- sgl[[i]]
    
    p <- array(NA,nInd(gl))  
    q <- array(NA,nInd(gl))
    Ho <- array(NA,nInd(gl))

    # For each individual
    for (j in 1:nInd(gl)){
      
      # Calculate heterozygosity as the proportion of heterozygous loci
      # First sum up the number of homozygous and heterozygous loci
      tbl <- table(as.matrix(gl)[j,],useNA = "always")
      if (is.na(tbl["0"])) {tbl["0"] <- 0}
      if (is.na(tbl["1"])) {tbl["1"] <- 0}
      if (is.na(tbl["2"])) {tbl["2"] <- 0}
    
      # Back-calculate the allele frequencies
      p[j] <- (2*tbl["0"] + tbl["1"]) # Raw
      q[j] <- (tbl["1"] + 2*tbl["2"])
      
      # Total alleles
      n <- p[j] + q[j]
      
      # Compute proportions
      p[j] <- p[j]/n # Proportions
      q[j] <- q[j]/n
    
      # Calculate the observed heterozygosity = proportion of heterozygous loci
      Ho[j] <- tbl["1"]/(n/2)
    }
    p <- as.numeric(p)
    q <- as.numeric(q)
    Hobs[i] <- mean(as.numeric(Ho),na.rm=TRUE)
  }
    
  df <- data.frame(popNames(x),as.numeric(table(pop(x))),Hobs,Hexp)
  names(df) <- c("pop","nInd","Ho","He")
  
  op <- par(mfrow=c(1,1),mai=c(1.2,0.5,0.2,0),oma=c(2,2,2,0), pty="m")
  df.ordered <- df[order(df$Ho),]
  barplot(df.ordered$Ho, names.arg=paste(df.ordered$pop, df.ordered$nInd, sep=" | "), las=2, cex.names=cex.labels, space=0, border=F, col=rainbow(nrow(df.ordered)), 
          main="Observed Heterozygosity by Population")
  
# PLOT OBSERVED AGAINST EXPECTED 
  
  # if (plot){
  #   par(mai=c(1,1,1,1),pty="s")
  #   p <- ggplot(df, aes(x=Hobs, y=Hexp))
  #   p <- p + geom_point(cex=2,col="red") 
  #   p <- p + geom_abline(intercept = 0, slope = 1)
  #   p <- p + xlab("Observed Heterozygosity") + ylab("Expected Heterozygosity")
  #   p <- p + theme(panel.background = element_rect(fill = "white"),
  #                plot.margin = margin(2, 2, 2, 2, "cm"),
  #                plot.background = element_rect(fill = "grey90",colour = "black",size = 1),
  #                axis.title=element_text(size=15),
  #                axis.text=element_text(size=15)
  #   )
  #   plot(p)
  # }  

# OUTPUT REPORT
  if (verbose >= 3){
    cat("Reporting Heterozygosity by Population\n")
    cat("No. of loci =", nLoc(x), "\n")
    cat("No. of individuals =", nInd(x), "\n")
    cat("No. of populations =", nPop(x), "\n")
  
    cat("  Miniumum Observed Heterozygosity: ",round(min(df$Ho,na.rm=TRUE),4),"\n")
    cat("  Maximum Observed Heterozygosity: ",round(max(df$Ho,na.rm=TRUE),4),"\n")
    cat("  Average Observed Heterozygosity: ",round(mean(df$Ho,na.rm=TRUE),4),"\n")
    cat("  Average Expected Heterozygosity: ",round(mean(df$He,na.rm=TRUE),4),"\n\n")
  
    print(df)
  } 
  
  }
  
  # DO THE JOB FOR INDIVIDUALS
  
  if (method=="ind"){
    # Convert to matrix
    m <- as.matrix(x)
    
    # For each individual determine counts of hets, homs and NAs
    c.na <- array(NA, nInd(x))
    c.hets <- array(NA, nInd(x))
    c.hom0 <- array(NA, nInd(x))
    c.hom2 <- array(NA, nInd(x))
    for (i in 1:nInd(x)){
      c.na[i] <- sum(is.na(m[i,]))
      c.hets[i] <- sum(m[i,]==1,na.rm=TRUE)/(nLoc(x)-c.na[i])
      c.hom0[i] <- sum(m[i,]==0,na.rm=TRUE)/(nLoc(x)-c.na[i])
      c.hom2[i] <- sum(m[i,]==2,na.rm=TRUE)/(nLoc(x)-c.na[i])
    }
    
    # Join the sample sizes with the heteozygosities
    df <- cbind.data.frame(x@ind.names, c.hets, c.hom0, c.hom2)
    names(df)<- c("ind.name", "c.hets", "c.hom0", "c.hom2")
    
    # Prepare for plotting
    # Save the prior settings for mfrow, oma, mai and pty, and reassign
      op <- par(mfrow = c(2, 1), oma=c(1,1,1,1), mai=c(0.5,0.5,0.5,0.5),pty="m")
    # Set margins for first plot
    par(mai=c(1,0.5,0.5,0.5))
    # Plot Box-Whisker plot
    adjbox(c.hets,
           horizontal = TRUE,
           col='red',
           range=range,
           main = "Heterozygosity by Individual")
    # Set margins for second plot
    par(mai=c(0.5,0.5,0,0.5))  
    # Plot Histogram
      hist(c.hets, col='red', main=NULL)
      
    # OUTPUT REPORT
      if (verbose >= 3){
        cat("Reporting Heterozygosity by Individual\n")
        cat("No. of loci =", nLoc(x), "\n")
        cat("No. of individuals =", nInd(x), "\n")

        cat("  Miniumum Observed Heterozygosity: ",round(min(df$c.hets),4),"\n")
        cat("  Maximum Observed Heterozygosity: ",round(max(df$c.hets),4),"\n")
        cat("  Average Observed Heterozygosity: ",round(mean(df$c.hets),4),"\n\n")
        
        #print(df)
      }   
  }  
    
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  # Reset the par options    
    par(op)
  
  # Return the result
  return(df) 
}
