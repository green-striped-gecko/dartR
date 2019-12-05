#' Reports observed and expected heterozygosity by population or by individual. 
#' 
#' Calculates the observed and expected heterozygosities for each population (method="pop") 
#' or the observed heterozyosity for each individual (method="ind") in a genlight object.
#' 
#' Observed heterozygosity for a population takes the proportion of heterozygous
#' loci for each individual then averages over the individuals in that population. 
#' The calculations take into account missing values.
#' 
#' Expected heterozygosity for a population takes the expected proportion of
#' heterozygotes, that is, expected under Hardy-Weinberg equilibrium, for each locus, then
#' averages this across the loci for an average estimate for the population.
#' The calculations of expected heterozygosity use the unbiassed estimates of Nei, M. (1987) 
#' Molecular evolutionary genetics. New York: Columbia University Press.
#'
#' Output for method="pop" is an ordered barchart of observed heterozygosity across 
#' populations together with a table of observed and expected heterozygosity by population.
#' 
#' Observed heterozygosity for individuals is calculated as the proportion of loci that
#' are heterozygous for that individual.
#' 
#' Output for method="ind" is a histogram of heterozygosity across individuals.
#' The histogram is accompanied by a box and whisker plot presented either in standard 
#' (boxplot="standard") or adjusted for skewness (boxplot=adjusted). 
#' 
#' Refer to Tukey (1977, Exploratory Data Analysis. Addison-Wesley) for standard
#' Box and Whisker Plots and Hubert & Vandervieren (2008, An Adjusted Boxplot for Skewed
#' Distributions, Computational Statistics & Data Analysis 52:5186-5201) for adjusted
#' Box and Whisker Plots.
#' 
#' Finally, the loci that are invariant across all individuals in the dataset (that is,
#' across populations), is typically unknown. This can render estimates of heterozygosity
#' analysis specific, and so it is not valid to compare such estimates across species
#' or even across different analyses. This is a similar problem faced by microsatellites.
#' If you have an estimate of the number of invariant sequence tags (loci) in your data,
#' such as provided by gl.report.secondaries, you can specify it with the n.invariant
#' parameter to standardize your estimates of heterozygosity.
#' 
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param method -- calculate heterozygosity by population (method='pop') or by individual (method='ind') [default 'pop']
#' @param n.invariant -- an estimate of the number of invariant sequence tags used to adjust the heterozygosity rate [default 0]
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

# Last amended 4-May-19

gl.report.heterozygosity <- function(x, 
                                     method="pop", 
                                     n.invariant=0,
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
  
  if(!is(x, "genlight")) {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  
  if (!(method=="pop" | method == "ind")) {
    cat("Warning: Method must either be by population or by individual, set to method='pop'\n")
    method <- "pop"   
  }
  
  if (n.invariant < 0){
    cat("Warning: Number of invariant loci must be non-negative, set to zero\n")
    n.invariant <- 0
  }
  
  if (!(boxplot=="standard" | boxplot == "adjusted")) {
    cat("Warning: Box-whisker plots must either standard or adjusted for skewness, set to boxplot='adjusted'\n")
    boxplot <- 'adjusted'   
  }


  # Set a population if none is specified (such as if the genlight object has been 
  # generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, 
                             individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
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
    if (verbose >=2){cat("  Calculating Observed Heterozygosities, averaged across loci, for each population\n")}
    
    # Parse across populations
  
    Hobs <- array(NA,length(sgl))
    Hobs.adj <- array(NA,length(sgl))
    obs.ind <- array(NA,length(sgl))
    obs.loc <- array(NA,length(sgl))
    for (i in 1:length(sgl)){
      gl <- sgl[[i]]
      
      # Parse across individuals
      
      Hobs.ind <- array(NA,nInd(gl))
      loc.obs.ind <- array(NA,nInd(gl))
      Hobs.ind.adj <- array(NA,nInd(gl))
      for (j in 1:nInd(gl)){
        
        # Parse across loci
      
        het.flag <- array(NA,nLoc(gl))
        for (k in 1:nLoc(gl)){
        
          # flag if the locus is heterozygous
          if (!is.na(as.matrix(gl)[j,k])){
            if (as.matrix(gl)[j,k] == 1){
              het.flag[k] <- 1
            } else {
              het.flag[k] <- 0
            } 
          }  else {
            het.flag[k] <- NA
          }
        }
        
      # Calculate the observed heterozygosity for each individual
      loc.obs.ind[j] <- sum(!is.na(het.flag))
      Hobs.ind[j] <- sum(het.flag, na.rm=TRUE)/loc.obs.ind[j]
      Hobs.ind.adj[j] <- sum(het.flag, na.rm=TRUE)/(loc.obs.ind[j] + n.invariant)
      
      }
      
      # Average the observed heterozygosity across individuals
      obs.ind[i] <- sum(!is.na(Hobs.ind))
      obs.loc[i] <- round(mean(loc.obs.ind,na.rm=TRUE),1)
      Hobs[i] <- sum(Hobs.ind, na.rm=TRUE)/obs.ind[i]
      Hobs.adj[i] <- Hobs[i]*obs.loc[i]/(obs.loc[i] + n.invariant)
      
    }  

    # EXPECTED HETEROZYGOSITY
    if (verbose >=2){cat("  Calculating Expected Heterozygosities, averaged across loci, for each population\n")}
    
    Hexp <- array(NA,length(sgl))
    Hexp.adj <- array(NA,length(sgl))
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
      Hexp.adj[i] <- Hexp[i]*obs.loc[i]/(obs.loc[i] + n.invariant)
    }  
    
  df <- data.frame(popNames(x),
                   as.numeric(table(pop(x))),
                   obs.loc,
                   round(Hobs,6),
                   round(Hobs.adj,6),
                   round(Hexp,6),
                   round(Hexp.adj,6)
                   )
  names(df) <- c("pop","nInd","nLoc","Ho","Ho.adj","He","He.adj")
  
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
  
    cat("  Miniumum Observed Heterozygosity: ",round(min(df$Ho,na.rm=TRUE),6),"\n")
    cat("  Maximum Observed Heterozygosity: ",round(max(df$Ho,na.rm=TRUE),6),"\n")
    cat("  Average Observed Heterozygosity: ",round(mean(df$Ho,na.rm=TRUE),6),"\n\n")
    
    cat("  Miniumum Observed Heterozygosity: ",round(min(df$He,na.rm=TRUE),6),"\n")
    cat("  Maximum Observed Heterozygosity: ",round(max(df$He,na.rm=TRUE),6),"\n")
    cat("  Average Observed Heterozygosity: ",round(mean(df$He,na.rm=TRUE),6),"\n\n")
    
    if (n.invariant > 0){
      cat("  Average correction factor for invariant loci =",mean(obs.loc/(obs.loc + n.invariant),rn.na=TRUE),"\n")
    }
  
    if (n.invariant > 0 ) {
      print(df)
    } else {
      print(df[,c("pop","nInd","nLoc","Ho","He")])
    }
 
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

        cat("  Miniumum Observed Heterozygosity: ",round(min(df$c.hets),6),"\n")
        cat("  Maximum Observed Heterozygosity: ",round(max(df$c.hets),6),"\n")
        cat("  Average Observed Heterozygosity: ",round(mean(df$c.hets),6),"\n\n")
        
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
